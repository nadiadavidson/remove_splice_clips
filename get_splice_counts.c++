/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <experimental/string_view>
#include <algorithm>

#include <sam.h>
#include <bam.h>

#include <gperftools/profiler.h>
using namespace std;
using namespace std::experimental;

static const int FLANK_SIZE=30;
static const bool ALLOW_MISMATCH=false;
static const int MIN_GAP=200000;
static const string END_LABEL="END";
static const string START_LABEL="START";
static const char EXON_ID_DELIM=':';
static const int MIN_COUNTS=2;

static int n_first_match=0;
static int n_perfect_match=0;

//complimenting code taken from stackoverflow
char compliment(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}


static class Counts {
  unordered_map< string, int > _counts;
  void print_if_interesting_junction(string name, int & read_support){
    //get positions and genes from the exon junction ids
    stringstream ss(name);    
    string field;
    vector<string> gene_info;
    while(getline(ss,field, EXON_ID_DELIM)) {
      gene_info.push_back(field);
    } //in case the format looks wrong.
    if(gene_info.size()!=8){
      cerr << "Issue with exon sequence IDs." 
	   << "Format should be Gene:Chrom:Start:END/START" << endl;
      exit(1);
    }

    //otherwise fill in the gene info
    vector<string> gene{gene_info.at(0),gene_info.at(4)};
    vector<string> chrom{gene_info.at(1),gene_info.at(5)};
    vector<int>pos{atoi(gene_info.at(2).c_str()),atoi(gene_info.at(6).c_str())};
    
    //now check if the junction looks interesting
    bool different_chrom = chrom[0]!=chrom[1];
    bool non_linear_order = pos[1] < pos[0];
    bool distal = ((pos[1]-pos[0])>MIN_GAP) & (gene[0]!=gene[1]);
    bool enough_support = read_support >= MIN_COUNTS;
    if( (different_chrom | non_linear_order | distal ) & enough_support ){
      cout << gene[0] << "\t" << chrom[0] << "\t" << pos[0] << "\t" 
	   << gene[1] << "\t" << chrom[1] << "\t" << pos[1] << "\t" 
	   << read_support << endl;
    }
  };
  
public:
  void increment(string& end,string& start){
    string pair = end + EXON_ID_DELIM + start;
    _counts[pair]++;
  };
  void print_table(){
    unordered_map< string , int >::iterator counts_itr=_counts.begin();
    for(;counts_itr!=_counts.end(); counts_itr++){
      print_if_interesting_junction(counts_itr->first,counts_itr->second);
    }
  };

} counts;

class JunctionSeq { //read the fasta
  unordered_map<string_view,string > junc_seq;
  unordered_map<string,string> _seqs;
public:
  void read_fasta( string & flank_fasta,const string type){
    ifstream file;
    file.open(flank_fasta);
    if(!(file.good())){
      cout << "Unable to open file " << flank_fasta << endl;
      exit(1);
    } //read the fasta files
    string id="";
    string line;
    while ( getline (file,line) ){
      int start=line.find(">")+1;
      if(start==1){ //if this is the ID line...
        int end=line.find_first_of("\t\n ")-1;
	id=line.substr(start,end);
      } else {
        _seqs[id]=_seqs[id]+line;
      }
    }
    //loop through the sequences and sort into start and end flanking sequence
    //mark any duplicate sequences for later removal
    vector<string_view> to_erase; //list of junction sequences that aren't unique.
    unordered_map<string,string>::iterator seq_itr=_seqs.begin();
    for(; seq_itr!=_seqs.end(); seq_itr++){
      //if more than one junction with this sequence
      //will need to remove later.
      string id=seq_itr->first;
      string_view seq=seq_itr->second;
      int second_last_char=id.size()-type.size()-1;
      bool is_type = id.find(type,second_last_char)!=string::npos;
      bool right_size = seq.size()==FLANK_SIZE;
      if(is_type & right_size){
	if(junc_seq.find(seq)!=junc_seq.end())
	  to_erase.push_back(seq);
	junc_seq[seq]=id;
      }
    }
    if(junc_seq.size()==0){
      cerr << "Found no compatible sequences in "<< flank_fasta << endl;
      exit(1);
    }
    sort(to_erase.begin(),to_erase.end());
    to_erase.erase(unique(to_erase.begin(),to_erase.end()),to_erase.end());
    //now loop again and remove all the black listed junctions
    for(int i=0; i<to_erase.size(); i++)
      junc_seq.erase(to_erase.at(i));
    cout << "Done reading fasta" << endl;
  };
  inline unordered_map<string_view,string>::iterator find(string_view & key){
    return junc_seq.find(key);
  };
  inline unordered_map<string_view,string>::iterator end(){
    return junc_seq.end();
  };
  
};

static JunctionSeq junc_seq_start;
static JunctionSeq junc_seq_end;

unordered_map<string_view,string>::iterator 
find_non_exact_match(string_view & orig_kmer,JunctionSeq& junc_seq){
  string kmer(orig_kmer);
  for(int base=0; base < FLANK_SIZE ; base++){
    vector<string> nuc{"A","G","C","T"};
    for(int n=0; n< nuc.size(); n++){
      kmer.replace(base,1,nuc.at(n));
      string_view sv_kmer=kmer;
      unordered_map<string_view,string>::iterator match=junc_seq.find(sv_kmer);
      if(match!=junc_seq.end())
	return match;
    }
  }
  return junc_seq.end();
}

//loop through the sequence to find a match
bool get_match(string & seq){
  if(seq.size()<(2*FLANK_SIZE)) return false;
  unordered_map<string_view,string>::iterator end; //end of exon1
  unordered_map<string_view,string>::iterator start; //joins to start of exon2
  string_view sv_seq = seq;
  string_view kmer1,kmer2;
  //search in the forward direction
  for(int pos=0; pos < (seq.size()-FLANK_SIZE) ; pos++){
    kmer1=sv_seq.substr(pos,FLANK_SIZE);
    end=junc_seq_end.find(kmer1);
    //if a match is found. look for other side of the junction
    if(end!=junc_seq_end.end()){
      n_first_match++;
      kmer2=sv_seq.substr(pos+FLANK_SIZE,FLANK_SIZE);
      start=junc_seq_start.find(kmer2);
      //if other end is found
      if(start!=junc_seq_start.end()){
	n_perfect_match++;
	//cout << "FOUND perfect match" << endl;
	counts.increment(end->second,start->second);
	return true;
      }
      //start not found. Try permutating the bases to account for 1 mismatch
      if(ALLOW_MISMATCH){
	start=find_non_exact_match(kmer2,junc_seq_start);
	if(start!=junc_seq_start.end()){
	  counts.increment(end->second,start->second);
	  return true;
	}
      }
    }
  }
  if(ALLOW_MISMATCH){
    //check again in reverse, permutating the end bases:
    for(int pos=seq.size()-FLANK_SIZE-1; pos >= FLANK_SIZE ; pos--){
      kmer1=sv_seq.substr(pos,FLANK_SIZE);
      start=junc_seq_start.find(kmer1);
      //if a match is found. look for other side of the junction
      if(start!=junc_seq_start.end()){
	kmer2=sv_seq.substr(pos-FLANK_SIZE,FLANK_SIZE);
	end=find_non_exact_match(kmer2,junc_seq_end);
	if(end!=junc_seq_end.end()){
	  counts.increment(end->second,start->second);
	  return true;
	}
      }
    }
  }
  return false;
}

int main(int argc, char *argv[]){

  if(argc!=3){
    cerr << "Usage: get_non_linear_region <exon_flanking_seq.fasta> <in.bam>" << endl;
    exit(1);
  }
  std::string flank_fasta=argv[1];
  std::string in_filename=argv[2];
  
  cerr << "Reading fasta file of junction sequences: " << flank_fasta << endl;
  junc_seq_start.read_fasta(flank_fasta,START_LABEL);
  junc_seq_end.read_fasta(flank_fasta,END_LABEL);

  ProfilerStart("prof.out");

  //Bam file reader
  samfile_t *in = 0 ;
  in = samopen(in_filename.c_str(), "br", NULL);
  if ((in==0) | (in->header == 0)) {
    cerr << "fail to open "<< in_filename << " for reading." << endl;
    exit(1);
  }
  bam1_t *b = bam_init1();
  int r;
  int i=0;
  int nread_processed=0;
  int f_count=0;
  int r_count=0;
  int nread=0;
  int unmapped=0;
  while ((r = samread(in, b)) >= 0) {
    nread++;
    if( nread % 1000000 == 0 ) cerr << nread/1000000 << " million reads processed" << endl;
    int seq_length=b->core.l_qseq;
    uint32_t *cigar = bam1_cigar(b);
    int matched=0;
    int rest=0;
    for(int k=0; k < b->core.n_cigar; k++){
      int c_oper=cigar[k]&BAM_CIGAR_MASK;
      int c_size=cigar[k]>>BAM_CIGAR_SHIFT;
      if(c_oper==BAM_CMATCH)
	matched+=c_size;
      else
	rest+=c_size;
    }
    if((seq_length-matched)<FLANK_SIZE && rest < MIN_GAP) continue ;
    nread_processed++;
    //get the read sequence
    char qseq[seq_length];
    uint8_t * s = bam1_seq(b);
    for(int n=0; n<seq_length; n++){
      char v = bam1_seqi(s,n);
      qseq[n] = bam_nt16_rev_table[v];
    }
    string seq(qseq,seq_length);
    //loop through the sequence to find the first match
    if(get_match(seq)){ 
      f_count++;
      //also check the reverse compliment if the read is unmapped.
    } else if (matched==0)  {
      unmapped++;
      reverse(seq.begin(),seq.end());
      transform(seq.begin(),seq.end(),seq.begin(),compliment);
      if(get_match(seq)){
	r_count++;
      }
    }
    //find all matches
  }
  bam_destroy1(b);
  samclose(in);


  ProfilerStop();

  cerr << "Reads Total=" << nread << endl;
  cerr << "Reads Processed=" << nread_processed << endl;
  cerr << "One match=" << n_first_match << endl;
  cerr << "Junction counts=" << f_count << "   " << r_count << endl;
  cerr << "Perfect matches=" << n_perfect_match << endl;
  cerr << "Unmapped=" << unmapped << endl;

  //print out the table of counts
  counts.print_table();

}

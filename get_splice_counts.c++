/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <numeric>  
#include <fstream>
#include <unordered_map>

#include <SeqLib/BamReader.h>
#include <SeqLib/GenomicRegionCollection.h>
#include <SeqLib/GenomicRegion.h>
#include <SeqLib/FastqReader.h>
#include <SeqLib/UnalignedSequence.h>

using namespace std;

static const int FLANK_SIZE=30;
static const int MIN_GAP=200000;
static const string END_LABEL="END";
static const string START_LABEL="START";
static const string EXON_ID_DELIM=":";

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
  const int MIN=1;
 public:
  void increment(string& end,string& start){
    string pair = end + "\t" + start;
    _counts[pair]++;
  };
  void print_table(){
    unordered_map< string , int >::iterator counts_itr=_counts.begin();
    for(;counts_itr!=_counts.end(); counts_itr++){
      if(counts_itr->second>MIN){
	cout << counts_itr->first << "\t"
	     << counts_itr->second << endl;
      }
    }
  };

} counts;

class JunctionSeq { //read the fasta
  unordered_map<string,string> junc_seq;
public:
  void read_fasta( string & flank_fasta,const string type){
    SeqLib::FastqReader fr;
    if(!fr.Open(flank_fasta)){
      cerr << "Trouble opening "<<flank_fasta << endl; 
      exit(1);
    }
    SeqLib::UnalignedSequence s;
    vector<string> to_erase; //list of junction sequences that aren't unique.
    while(fr.GetNextSequence(s)){
      //if more than one junction with this sequence
      //will need to remove later.
      bool is_type = s.Name.find(type,s.Name.size()-type.size()-1)!=string::npos; //at the end?
      bool right_size = s.Seq.size()==FLANK_SIZE;
      if(is_type & right_size & (junc_seq.find(s.Seq)!=junc_seq.end()))
	to_erase.push_back(s.Seq);
      junc_seq[s.Seq]=s.Name;
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
  };
  inline unordered_map<string,string>::iterator find(string & key){
    return junc_seq.find(key);
  };
  inline unordered_map<string,string>::iterator end(){
    return junc_seq.end();
  };
  
};

static JunctionSeq junc_seq_start;
static JunctionSeq junc_seq_end;

unordered_map<string,string>::iterator 
find_non_exact_match(string& kmer,JunctionSeq& junc_seq){
  for(int base=0; base < FLANK_SIZE ; base++){
    vector<string> nuc{"A","G","C","T"};
    for(int n=0; n< nuc.size(); n++){
      kmer.replace(base,1,nuc.at(n));
      unordered_map<string,string>::iterator match=junc_seq.find(kmer);
      if(match!=junc_seq.end())
	return match;
    }
  }
  return junc_seq.end();
}

//loop through the sequence to find a match
bool get_match(string & seq){
  if(seq.size()<(2*FLANK_SIZE)) return false;
  unordered_map<string,string>::iterator end; //end of exon1
  unordered_map<string,string>::iterator start; //joins to start of exon2
  //search in the forward direction
  for(int pos=0; pos < (seq.size()-FLANK_SIZE) ; pos++){
    string kmer1=seq.substr(pos,FLANK_SIZE);
    end=junc_seq_end.find(kmer1);
    //if a match is found. look for other side of the junction
    if(end!=junc_seq_end.end()){
      n_first_match++;
      string kmer2=seq.substr(pos+FLANK_SIZE,FLANK_SIZE);
      start=junc_seq_start.find(kmer2);
      //if other end is found
      if(start!=junc_seq_start.end()){
	n_perfect_match++;
	//cout << "FOUND perfect match" << endl;
	counts.increment(end->second,start->second);
	return true;
      }
      //start not found. Try permutating the bases to account for 1 mismatch
      start=find_non_exact_match(kmer2,junc_seq_start);
      if(start!=junc_seq_start.end()){
	counts.increment(end->second,start->second);
	return true;
      }
    }
  }
  //check again in reverse, permutating the end bases:
  for(int pos=seq.size()-FLANK_SIZE-1; pos >= FLANK_SIZE ; pos--){
    string kmer1=seq.substr(pos,FLANK_SIZE);
    start=junc_seq_start.find(kmer1);
    //if a match is found. look for other side of the junction
    if(start!=junc_seq_start.end()){
      string kmer2=seq.substr(pos-FLANK_SIZE,FLANK_SIZE);
      end=find_non_exact_match(kmer2,junc_seq_end);
      if(end!=junc_seq_end.end()){
	counts.increment(end->second,start->second);
	return true;
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
  
  //Bam file reader
  SeqLib::BamReader bw;
  if(!bw.Open(in_filename)){
    cerr << "Trouble opening "<< in_filename << endl;
    exit(1);
  }
  SeqLib::BamRecord r;
  //loop through bam records
  cerr << "Reading bam file:" << in_filename << endl;
  int nread=0;
  int nread_processed=0;
  int f_count=0;
  int r_count=0;
  while(bw.GetNextRecord(r)){
    nread++;
    if( nread % 100000 == 0 ) cerr << nread << endl;
    if((r.Length()-r.NumAlignedBases())<FLANK_SIZE && (r.PositionEnd() - r.Position()) < MIN_GAP){
      //      cout << r.CigarString() << endl;
      continue;
    }
    nread_processed++;
    string seq=r.Sequence();
    //loop through the sequence to find the first match
    if(get_match(seq)){ 
      f_count++;
      //no need to reverse compliment because the bam sequence has already been done.
      /**   } else {
      reverse(seq.begin(),seq.end());
      transform(seq.begin(),seq.end(),seq.begin(),compliment);
      if(get_match(seq,junc_seq,counts)){
	cout << "Break found reverse - " << r.ReverseFlag() << endl;
	r_count++;
	}**/
    }
    //find all matches
  }
  bw.Close();

  cerr << "Reads Total=" << nread << endl;
  cerr << "Reads Processed=" << nread_processed << endl;
  cerr << "One match=" << n_first_match << endl;
  cerr << "Junction counts=" << f_count << "   " << r_count << endl;
  cerr << "Perfect matches=" << n_perfect_match << endl;

  //print out the table of counts
  counts.print_table();

}

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

static int n_first_match=0;
static int n_perfect_match=0;

//complimenting code taken from stackoverflow
char compliment(char c){
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
  bool read_fasta( string & flank_fasta){
    SeqLib::FastqReader fr;
    if(!fr.Open(flank_fasta)) return false;
    SeqLib::UnalignedSequence s;
    vector<string> to_erase; //list of junction sequences that aren't unique.
    while(fr.GetNextSequence(s)){
      //if more than one junction with this sequence
      //will need to remove later.
      if(junc_seq.find(s.Seq)!=junc_seq.end())
	to_erase.push_back(s.Seq);
      junc_seq[s.Seq]=s.Name;
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

//loop through the sequence to find a match
bool get_match(string & seq){
  if(seq.size()<(2*FLANK_SIZE)) return false;
  unordered_map<string,string>::iterator end; //end of exon1
  unordered_map<string,string>::iterator start; //joins to start of exon2
  //search in the forward direction
  for(int pos=0; pos < (seq.size()-FLANK_SIZE) ; pos++){
    string kmer=seq.substr(pos,FLANK_SIZE);
    end=junc_seq_end.find(kmer);
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
      for(int base=0; base < FLANK_SIZE ; base++){
	vector<string> nuc{"A","G","C","T"};
	for(int n=0; n< nuc.size(); n++){
	  kmer2.replace(base,1,nuc.at(n));
	  start=junc_seq_start.find(kmer2);
	  if(start!=junc_seq_start.end()){
	    counts.increment(end->second,start->second);
	    return true;
	  }
	}
      }
    }
  }
  //check again in reverse, permutating the end bases:
  for(int pos=seq.size()-FLANK_SIZE-1; pos >= FLANK_SIZE ; pos--){
    string kmer=seq.substr(pos,FLANK_SIZE);
    start=junc_seq_start.find(kmer);
    //if a match is found. look for other side of the junction
    if(start!=junc_seq_start.end()){
      string kmer2=seq.substr(pos-FLANK_SIZE,FLANK_SIZE);
      for(int base=0; base < FLANK_SIZE ; base++){
        vector<string> nuc{"A","G","C","T"};
        for(int n=0; n< nuc.size(); n++){
          kmer2.replace(base,1,nuc.at(n));
          end=junc_seq_end.find(kmer2);
          if(end!=junc_seq_end.end()){
            counts.increment(end->second,start->second);
            return true;
          }
        }
      }
    }
  }
  return false;
}

int main(int argc, char *argv[]){

  if(argc!=4){
    cerr << "Usage: get_non_linear_region <exon_starts.fasta> <exon_ends.fasta> <in.bam>" << endl;
    exit(1);
  }
  std::string flank_start_fasta=argv[1];
  std::string flank_end_fasta=argv[2];
  std::string in_filename=argv[3];
  
  cerr << "Reading fasta files of junction sequences: " << flank_start_fasta
       << " and " << flank_end_fasta << endl;
  if(!junc_seq_start.read_fasta(flank_start_fasta)){
    cerr << "Trouble opening "<<flank_start_fasta << endl; exit(1);
  }
  if(!junc_seq_end.read_fasta(flank_end_fasta)){
    cerr << "Trouble opening "<<flank_end_fasta << endl; exit(1);
  }
  
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

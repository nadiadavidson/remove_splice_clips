/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <numeric>  
#include <fstream>

#include <SeqLib/BamReader.h>
#include <SeqLib/GenomicRegionCollection.h>
#include <SeqLib/GenomicRegion.h>
#include <SeqLib/FastqReader.h>
#include <SeqLib/UnalignedSequence.h>

using namespace std;

//code taken from stackoverflow
char compliment(char c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

//loop through the sequence to find a match
bool get_match(string & seq, map<string,string> & junc_seq, map< pair<string,string > , int > & counts ){
  if(seq.size()<60) return false;
  for(int pos=0; pos < (seq.size()-30) ; pos++){
    string k_mer=seq.substr(pos,30);
    map<string,string>::iterator match=junc_seq.find(k_mer);
    //a match is found. look for other side of the junction
    if(match!=junc_seq.end()){
      string k_mer2=seq.substr(pos+30,30);
      map<string,string>::iterator match2=junc_seq.find(k_mer2);
      //other end if found
      if(match2!=junc_seq.end()){
	counts[make_pair(match->second,match2->second)]++;
	//	cout << "Found pair: " << match->second << "   " << match2->second << endl;
	return true;
      }
    }
  }
  return false;
}

int main(int argc, char *argv[]){

  if(argc!=3){
    cerr << "Usage: get_non_linear_region <exon_flanking_sequence.fasta> <in.bam>" << endl;
    exit(1);
  }
  std::string flank_fasta=argv[1];
  std::string in_filename=argv[2];

  //Bam file reader
  SeqLib::BamReader bw;
  if(!bw.Open(in_filename)){
    cerr << "Trouble opening "<< in_filename << endl;
    exit(1);
  }
  SeqLib::BamRecord r;
  
  //read the fasta
  cerr << "Reading fasta file of junction sequence: " << flank_fasta << endl;
  SeqLib::FastqReader fr;
  if(!fr.Open(flank_fasta)){ 
    cerr << "Trouble opening "<<flank_fasta << endl;
    exit(1);
  }
  SeqLib::UnalignedSequence s;
  map<string,string> junc_seq;
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

  //loop through bam records
  cerr << "Reading bam file:" << in_filename << endl;
  int nread=0;
  int nread_processed=0;
  int f_count=0;
  int r_count=0;
  map< pair<string,string > , int > counts;
  while(bw.GetNextRecord(r)){
    nread++;
    if( nread % 100000 == 0 ) cerr << nread << endl;
    if((r.Length()-r.NumAlignedBases())<30 && (r.PositionEnd() - r.Position()) < 100000) continue;
    nread_processed++;
    string seq=r.Sequence();
    //loop through the sequence to find the first match
    if(get_match(seq,junc_seq,counts)){ 
      f_count++;
    } else {
      reverse(seq.begin(),seq.end());
      transform(seq.begin(),seq.end(),seq.begin(),compliment);
      if(get_match(seq,junc_seq,counts))
	r_count++;
    }
    //find all matches
  }
  cerr << "Junction counts=" << f_count << "   " << r_count << endl;
  cerr << "Reads Total=" << nread << endl;
  cerr << "Reads Processed=" << nread_processed << endl;

  //print out the table of counts
  map< pair<string,string > , int >::iterator counts_itr=counts.begin();
  for(;counts_itr!=counts.end(); counts_itr++){
    if(counts_itr->second>1){
      cout << counts_itr->first.first << "\t"
	   << counts_itr->first.second << "\t"
	   << counts_itr->second << endl;
    }
  }
  bw.Close();
}


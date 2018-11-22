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
bool get_match(string & seq, map<string,string> & junc_seq){
  for(int pos=0; pos < (seq.size()-30) ; pos++){
    string k_mer=seq.substr(pos,30);
    map<string,string>::iterator match=junc_seq.find(k_mer);
    //a match is found. look for other side of the junction
    if(match!=junc_seq.end()){
      string k_mer2=seq.substr(pos+30,30);
      map<string,string>::iterator match2=junc_seq.find(k_mer2);
      //other end if found
      if(match2!=junc_seq.end()){
	cout << "Found pair: " << match->second << "   " << match2->second << endl;
	return true;
      }
    }
  }
  return false;
}

int main(int argc, char *argv[]){

  if(argc!=3){
    cout << "Usage: get_non_linear_region <exon_flanking_sequence.fasta> <in.bam>" << endl;
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
  SeqLib::FastqReader fr;
  if(!fr.Open(flank_fasta)){ 
    cerr << "Trouble opening "<<flank_fasta << endl;
    exit(1);
  }
  SeqLib::UnalignedSequence s;
  map<string,string> junc_seq;
  while(fr.GetNextSequence(s)){
    junc_seq[s.Seq]=s.Name;
  }

  //loop through bam records
  cout << "Reading bam file" << endl;
  while(bw.GetNextRecord(r)){
    string seq=r.Sequence();
    //loop through the sequence to find the first match
    if(!get_match(seq,junc_seq)){
      cout << seq << " ";
      reverse(seq.begin(),seq.end());
      transform(seq.begin(),seq.end(),seq.begin(),compliment);
      cout << seq << endl;
      if(get_match(seq,junc_seq))
	cout << "Rev match found" << endl;
    }
    //find all matches
  }

  bw.Close();
}


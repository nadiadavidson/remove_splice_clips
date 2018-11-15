/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <numeric>

#include <SeqLib/BamReader.h>

using namespace std;


static const int BIN_SIZE=100000;

int main(int argc, char *argv[]){

  if(argc!=2){
    cout << "Usage: get_non_linear_region <in.bam>" << endl;
    exit(1);
  }
  std::string in_filename=argv[1];

  //Bam file reader
  SeqLib::BamReader bw;
  bw.Open(in_filename);
  SeqLib::BamRecord r;

  //loop through bam records
  cout << "Reading bam file" << endl;
  //first iteration is to check for barcode contamination
  map<string,float> density;
  float bad_map=0;
  float all_map=0;
  while(bw.GetNextRecord(r) && r.MappedFlag()){
    all_map++;
    if(!r.ProperPair()) bad_map++;
    if(all_map==BIN_SIZE){
      stringstream position;
      position << r.ChrName() << ":" << r.Position() << endl;
      density[position.str()]=bad_map/all_map;
      if(bad_map>10000)
	cout << position.str() << " " << bad_map << "/" << all_map << endl;
      all_map=0;
      bad_map=0;
    }
  }
  bw.Close();

}


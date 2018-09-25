/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <numeric>

#include <SeqLib/BamReader.h>
#include <SeqLib/BamWriter.h>
#include <SeqLib/BWAWrapper.h>
using namespace std;

static const int THREADS=8;
static const int MIN_TRANS_MAP_SC_LENGTH=12;

int main(int argc, char *argv[]){

  std::string in_filename=argv[1];
  std::string out_filename=argv[2];
  std::string bwa_index=argv[3];

  //Bam file reader
  SeqLib::BamReader bw;
  SeqLib::BamRecord read;
  bw.Open(in_filename);

  //open the output BAM
  SeqLib::BamWriter writer; // or writer(SeqLib::SAM) or writer(SeqLib::CRAM) 
  writer.SetHeader(bw.Header());
  writer.Open(out_filename);
  SeqLib::ThreadPool t(THREADS);
  writer.SetThreadPool(t);
  writer.WriteHeader();

  //open BWA Wrapper
  SeqLib::BWAWrapper bwa;
  bwa.LoadIndex(bwa_index); 

  while(bw.GetNextRecord(read)){
    bool keep_read=true;
    if(read.NumSoftClip()>MIN_TRANS_MAP_SC_LENGTH){
      SeqLib::BamRecordVector results;
      bwa.AlignSequence(read.Sequence(), read.Qname(), results,false,1.0,0);
      //cout << "Soft clipped:"<< r.NumSoftClip() << " "; 
      for(int a=0; a< results.size(); a++){
       	//cout << "Realigned: " << results[a].NumSoftClip() << " , " ;
	if(results[a].NumSoftClip()<MIN_TRANS_MAP_SC_LENGTH) keep_read=false; //remove is alignment without softclip found
      }
    }
    if(keep_read) writer.WriteRecord(read);
    //    cout << endl;
  }
  bw.Close();
}


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

//static const char SOFT_CLIP_TAG='S';
static const int MAX_CLIP_LENGTH=150;
static const int CONTAM_SAMPLE_SIZE=50000;
static const int THREADS=8;

int check_contam(vector<int> clip_table){
  int max=*max_element(clip_table.begin()+1,clip_table.end());
  int total=accumulate(clip_table.begin(),clip_table.end(),0);
  int length=1;
  while(clip_table.at(length)!=max) length++;
  if(max>0.2*total){
    cout << length << "bp contamination suspected (max/total:"<<max<<"/"<<total<<")" << endl;
    return length ;
  }
  cout << "no contamination found (max/total:"<<max<<"/"<<total<<"),"<<length<< endl;
  return 0;
}

SeqLib::Cigar fix_cigar( SeqLib::Cigar& cigar){
  SeqLib::Cigar new_cigar;
  new_cigar.add(SeqLib::CigarField(cigar[1].Type(),cigar[1].Length()+cigar[0].Length()));
  for(int s=2; s<cigar.size() ; s++)
    new_cigar.add(SeqLib::CigarField(cigar[s].Type(),cigar[s].Length()));
  return new_cigar;
}

int main(int argc, char *argv[]){

  if(argc!=3){
    cout << "Usage: fix_contam_clips <in.bam> <out.bam>" << endl;
    exit(1);
  }
  std::string in_filename=argv[1];
  std::string out_filename=argv[2];

  //Bam file reader
  SeqLib::BamReader bw;
  bw.Open(in_filename);
  SeqLib::BamRecord r;
  //loop through bam records
  //first iteration is to check for barcode contamination
  vector<int> SC_length_R1(MAX_CLIP_LENGTH);
  vector<int> SC_length_R2(MAX_CLIP_LENGTH);
  cout << "Checking for sequence contamination at start of reads" << endl; 
  int contam_sample_reads=0;
  while (bw.GetNextRecord(r) && contam_sample_reads <= CONTAM_SAMPLE_SIZE){
    int start = r.AlignmentPosition();
    if(r.ReverseFlag())
      start = r.AlignmentPositionReverse();
    if(r.FirstFlag()){
      SC_length_R1[start]++;
    } else 
      SC_length_R2[start]++;
    contam_sample_reads++;
  }
  cout << "Read end 1:" << endl;
  int r1_snip=check_contam(SC_length_R1);
  cout << "Read end 2:" << endl;
  int r2_snip=check_contam(SC_length_R2);
  bw.Close();

  //Now loop through the reads again
  //remove soft clipping flag from reads with contamination
  //remap soft clipped reads to transcriptome.
  //open the output BAM
  SeqLib::BamWriter writer; // or writer(SeqLib::SAM) or writer(SeqLib::CRAM) 
  writer.SetHeader(bw.Header());
  writer.Open(out_filename);
  SeqLib::ThreadPool t(THREADS);
  writer.SetThreadPool(t);
  writer.WriteHeader();

  //open BWA Wrapper
  //SeqLib::BWAWrapper bwa;
  //bwa.LoadIndex(bwa_index); 

  //Bam file reader again
  SeqLib::BamReader bw2;
  SeqLib::BamRecord read;
  bw2.Open(in_filename);
  while (bw2.GetNextRecord(read)){
    if((read.NumSoftClip())>0){
      //If this is a case of contamination, remove the tag (make them mismatches)
      int start = read.AlignmentPosition();
      if(read.ReverseFlag())
	start = read.AlignmentPositionReverse();
      if((r1_snip!=0 && read.FirstFlag() && start==r1_snip) ||
	 (r2_snip!=0 && !read.FirstFlag() && start==r2_snip)){
	//do something...
	SeqLib::Cigar cigar=read.GetCigar();
	SeqLib::Cigar new_cigar;
	if(!read.ReverseFlag()){
	  new_cigar = fix_cigar(cigar);
	} else {
	  reverse(cigar.begin(),cigar.end());
	  new_cigar=fix_cigar(cigar);
	  reverse(new_cigar.begin(),new_cigar.end());
	}
	read.SetCigar(new_cigar);
      }
    }
    writer.WriteRecord(read);
    //    cout << endl;
  }
  bw2.Close();
}


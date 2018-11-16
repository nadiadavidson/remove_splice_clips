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

using namespace std;

int main(int argc, char *argv[]){

  if(argc!=3){
    cout << "Usage: get_non_linear_region <genes.bed> <in.bam>" << endl;
    exit(1);
  }
  std::string gene_bed=argv[1];
  std::string in_filename=argv[2];

  //Bam file reader
  SeqLib::BamReader bw;
  bw.Open(in_filename);
  SeqLib::BamRecord r;

  //read the bed file and make gene regions for counting.
  //read the bed file
  ifstream infile(gene_bed);
  map<string,SeqLib::GenomicRegionCollection<SeqLib::GenomicRegion>> gene_regions;
  string line;
  while (getline(infile, line)){
    if(!line.find("#") == 0){
      stringstream ss(line);
      string chrom; ss >> chrom;
      string start; ss >> start;
      string end; ss >> end;
      string name; ss >> name;
      SeqLib::GenomicRegion exon_gr(chrom,start,end,bw.Header());
      if(start!=end)
	gene_regions[name].add(exon_gr);  
    }
  }
  // now loop through the genes and merge overlapping intervals
  map<string,SeqLib::GenomicRegionCollection<SeqLib::GenomicRegion>>::iterator gene_regions_itr;
  SeqLib::GenomicRegionCollection<SeqLib::GenomicRegion> gene_gr;
  for(gene_regions_itr=gene_regions.begin(); gene_regions_itr!=gene_regions.end(); gene_regions_itr++){
    gene_regions_itr->second.MergeOverlappingIntervals();
    gene_gr.Concat(gene_regions_itr->second);
  }
  gene_gr.CreateTreeMap();
  //now loop through the gene_ranges again to get the index
  map<string,vector<int>> gene_to_range_indx;
  for(gene_regions_itr=gene_regions.begin(); gene_regions_itr!=gene_regions.end(); gene_regions_itr++){
    for(int i=0; i<gene_regions_itr->second.size(); i++){
      for(int j=0; j<gene_gr.size(); j++){
	if(gene_gr[j]==gene_regions_itr->second[i]){
	  gene_to_range_indx[gene_regions_itr->first].push_back(j);
	  break;
	}
      }
    }
    /**gene_regions_itr->second.CreateTreeMap();
    vector<int> query_id;
    vector<int> subject_id;
    gene_regions_itr->second.FindOverlaps(gene_gr,query_id,subject_id,true);
    string name=gene_regions_itr->first;
    gene_to_range_indx[name]=subject_id;**/
  }

  //loop through bam records
  cout << "Reading bam file" << endl;
  //first iteration is to check for barcode contamination
  map<int,int> all_reads;
  map<int,int> pair_unmapped;
  map<int,int> pair_other_chrom;
  map<int,int> pair_proper;
  map<int,int> pair_bad_order;

  while(bw.GetNextRecord(r) && r.MappedFlag()){
    //find the gene interval
    SeqLib::GenomicRegion read_gr=r.AsGenomicRegion();
    vector<int> gene_indx=gene_gr.FindOverlappedIntervals(read_gr,true);
    for(int i=0; i < gene_indx.size(); i++){
      int indx=gene_indx[i];
      all_reads[indx]++;
      if(!r.PairMappedFlag()) pair_unmapped[indx]++;
      if(r.Interchromosomal()) pair_other_chrom[indx]++;
      if(!r.ProperPair()) pair_proper[indx]++;
      if(!r.ProperOrientation()) pair_bad_order[indx]++;
    }
  }

  //output the counts
  map<string,vector<int>>::iterator gene_to_range_indx_itr;
  for(gene_to_range_indx_itr=gene_to_range_indx.begin();
      gene_to_range_indx_itr!=gene_to_range_indx.end();
      gene_to_range_indx_itr++){
    cout << gene_to_range_indx_itr->first << "\t" ;
    vector<int> gene_indx=gene_to_range_indx_itr->second;
    //add up the counts if the count is found in mutiple places on the genome
    int all=0, unmapped=0, other_chrom=0, proper=0, bad_order=0;
    for(int i=0; i<gene_indx.size(); i++){
      int index=gene_indx.at(i);
      all+=all_reads[index]; unmapped+=pair_unmapped[index];
      other_chrom+=pair_other_chrom[index]; proper+=pair_proper[index];
      bad_order+=pair_bad_order[index];
    }
    cout << all << "\t" << unmapped << "\t" << other_chrom;
    cout << "\t" << proper << "\t" << bad_order << endl;
  }
  bw.Close();
}


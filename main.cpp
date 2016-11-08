#include <algorithm>
#include <iostream>
#include <getopt.h>
#include <stdlib.h>

#include "lca.h"
#include "kmer_counter.h"
#include "read_stitcher.h"

int    max_read_len;
int    max_k;
int    min_bp_overlap;
double min_frac_correct;
int    check_entropy;
double min_entropy;

void printStitching(std::string& s1, std::string& s2, int index){
  std::cout << s1 << std::endl;
  std::string spacing = "";
  for (int i = 0; i < index; i++)
    spacing += " ";
  std::cout << spacing << s2 << std::endl;
}

void generateSequences(int num_seqs, int seq_length, std::vector<std::string>& seqs){
  srand(1187238845);
  char bases[4] = {'A', 'C', 'G', 'T'};
  seqs.clear();
  for(int i = 0; i < num_seqs; i++){
    std::string res = "";
    for(int j = 0; j < seq_length; j++){
      res += bases[rand()%4];
    }
    seqs.push_back(res);
  }
}

void mutate(std::vector<std::string>& reads, double mutation_rate){
  mutation_rate *= 4.0/3.0; // 1 out of 4 times we mutate the base to itself
  char bases[4]  = {'A', 'C', 'G', 'T'};
  for(int i = 0; i < reads.size(); i++){
    for(int j = 0; j < reads[i].size(); j++){
      double f = ((double)rand())/RAND_MAX;
      if (f <= mutation_rate)
	reads[i][j] = bases[rand()%4];
    }
  }
}

void generateReads(int num_reads, int read_length, std::vector<std::string>& reads_a, std::vector<std::string>& reads_b, 
		   std::vector<int>& pos_a, std::vector<int>& pos_b, double mutation_rate){
  std::vector<std::string> seqs;
  int seq_length = 2*read_length;
  generateSequences(num_reads, seq_length, seqs);
  reads_a.clear();
  reads_b.clear();
  pos_a.clear();
  pos_b.clear();
  for(int i = 0; i < num_reads; i++){
    int loc_a = rand()%(seq_length-read_length+1);
    int loc_b = rand()%(seq_length-read_length-loc_a+1) + loc_a;
    reads_a.push_back(seqs[i].substr(loc_a, read_length));
    reads_b.push_back(seqs[i].substr(loc_b, read_length));
    pos_a.push_back(loc_a);
    pos_b.push_back(loc_b);
  }
  
  mutate(reads_a, mutation_rate);
  mutate(reads_b, mutation_rate);
}


int main(int argc, char **argv){
  // Default parameter settings
  max_read_len     = 100;
  max_k            = 10;
  min_bp_overlap   = 10;
  min_frac_correct = 0.9;
  check_entropy    = 1;
  min_entropy      = 0.0;
	
  std::string f1 = "";
  std::string f2 = "";
  
  if (argc == 1){
    std::cout 
      << "Usage: ReadStitcher --f1 <fq_1> --f2 <fq_2> [options]" << std::endl
      << "\t" << "--check-entropy"              << std::endl
      << "\t" << "--min-entropy      <FLOAT>"   << " Minimum entropy" << std::endl
      << "\t" << "--min-frac-correct <FLOAT>"   << " Minimum fraction of overlapping bases that must match"  << std::endl
      << "\t" << "--max-read-length  <INT>  "   << " Maximum read length to be considered"                   << std::endl
      << "\t" << "--max-mismatches   <INT>  "   << " Maximum number of overlapping bases that can not match" << std::endl
      << "\t" << "--min-overlap      <INT>  "   << " Minimum number of overlapping bases required "          << std::endl;
    exit(0);
  }

  static struct option long_options[] = {
    {"check-entropy",    no_argument,  &check_entropy, 1},
    {"f1",               required_argument, 0, 'a'},
    {"f2",               required_argument, 0, 'b'},
    {"min-entropy",      required_argument, 0, 'e'},
    {"min-frac-correct", required_argument, 0, 'f'},
    {"max-read-length",  required_argument, 0, 'l'},
    {"max-mismatches",   required_argument, 0, 'm'},
    {"min-overlap",      required_argument, 0, 'o'}
  };

  int c;
  while (1){
    int option_index = 0;
    c = getopt_long(argc, argv, "a:b:e:f:l:m:o:", long_options, &option_index);
    if (c == -1)
      break;
    switch (c){
    case 0:
      break;
    case 'a':
      f1 = std::string(optarg);
      break;
    case 'b':
      f2 = std::string(optarg);
      break;
    case 'e':
      min_entropy = atof(optarg);
      break;
    case 'f':
      min_frac_correct = atof(optarg);
      break;
    case 'l':
      max_read_len = atoi(optarg);
      break;
    case 'm':
      max_k = atoi(optarg);
      break;
    case 'o':
      min_bp_overlap = atoi(optarg);
      break;
    case '?':
      break;
    default:
      abort();
      break;
    }
  }
  

  if (f1.empty() || f2.empty()){
    std::cerr << "ERROR: Missing input file(s). Exiting..." << std::endl;
    exit(1);
  }
  
  ReadStitcher stitcher(max_read_len, max_k, min_bp_overlap, min_frac_correct);
  
  std::vector<std::string> l_reads;
  std::vector<std::string> r_reads;
  std::vector<int>         l_start;
  std::vector<int>         r_start;
  
  //100,100
  
  generateReads(10000, max_read_len, l_reads, r_reads, l_start, r_start, 0.01);
  for(int i = 0; i < l_reads.size(); i++){
    int best_frac_idx = stitcher.stitch_reads(l_reads[i], r_reads[i]);
    if (best_frac_idx != -1){
      printStitching(l_reads[i], r_reads[i], best_frac_idx);
      std::cout << std::endl;
    }
  }
}

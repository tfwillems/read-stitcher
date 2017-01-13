#include <algorithm>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <stdlib.h>

#include "error.h"
#include "kmer_counter.h"
#include "lca.h"
#include "read_stitcher.h"
#include "stringops.h"
#include "version.h"

int    max_read_len;
int    max_k;
int    min_bp_overlap;
double min_frac_correct;

bool file_exists(std::string path){
  return (access(path.c_str(), F_OK) != -1);
}

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

void print_usage(){
   std::cout
      << "Usage: ReadStitcher --f1 <fq_1.gz> --f2 <fq_2.gz> --out <prefix> --log <log_file.txt> [options]"                        << "\n"
      << "\t" << "--f1               <fq_1.gz>      " << "\t" << " Bgzipped FASTQ containing first  set of reads"                 << "\n"
      << "\t" << "--f2               <fq_2.gz>      " << "\t" << " Bgzipped FASTQ containing second set of reads"                 << "\n"
      << "\t" << "--out              <prefix>       " << "\t" << " Prefix for output files for stitched and unstitched reads"     << "\n"
      << "\t" << "--log              <log_file.txt> " << "\t" << " Path for log file output"                                      << "\n"
      << "\t" << "--min-frac-correct <FLOAT>        " << "\t" << " Minimum fraction of overlapping bases that must match (Default = "  << min_frac_correct << ")" << "\n"
      << "\t" << "--max-read-length  <INT>          " << "\t" << " Maximum read length to be considered (Default = "                   << max_read_len     << ")" << "\n"
      << "\t" << "--max-mismatches   <INT>          " << "\t" << " Maximum number of overlapping bases that can not match (Default = " << max_k            << ")" << "\n"
      << "\t" << "--min-overlap      <INT>          " << "\t" << " Minimum number of overlapping bases required (Default = "           << min_bp_overlap   << ")" << "\n"
      << "\t" << "--help                            " << "\t" << " Print this help message and exit"                                                              << "\n"
      << "\t" << "--version                         " << "\t" << " Print ReadStitcher version and exit"                                                           << "\n" << std::endl;
    exit(0);
}


int main(int argc, char **argv){
  // Default parameter settings
  max_read_len      = 100;
  max_k             = 10;
  min_bp_overlap    = 10;
  min_frac_correct  = 0.9;
  std::string f1    = "";
  std::string f2    = "";
  std::string out   = "";
  std::string log   = "";
  int print_version = 0, print_help = 0;
  
  if (argc == 1)
    print_usage();

  static struct option long_options[] = {
    {"f1",               required_argument, 0, 'a'},
    {"f2",               required_argument, 0, 'b'},
    {"min-frac-correct", required_argument, 0, 'f'},
    {"max-read-length",  required_argument, 0, 'l'},
    {"max-mismatches",   required_argument, 0, 'm'},
    {"min-overlap",      required_argument, 0, 'o'},
    {"out",              required_argument, 0, 'p'},
    {"log",              required_argument, 0, 'r'},
    {"help",        no_argument, &print_help,    1},
    {"version",     no_argument, &print_version, 1},
    {0, 0, 0, 0}
  };

  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "a:b:f:l:m:o:", long_options, &option_index);
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
    case 'p':
      out = std::string(optarg);
      break;
    case 'r':
      log = std::string(optarg);
      break;
    case '?':
      printErrorAndDie("Unrecognized command line option");
      break;
    default:
      abort();
      break;
    }
  }

  if (optind < argc) {
    std::stringstream msg;
    msg << "Did not recognize the following command line arguments:" << "\n";
    while (optind < argc)
      msg << "\t" << argv[optind++] << "\n";
    msg << "Please check your command line syntax or type ./ReadStitcher --help for additional information" << "\n";
    printErrorAndDie(msg.str());
  }

  if (print_version == 1){
    std::cerr << "ReadStitcher version " << VERSION << std::endl;
    exit(0);
  }

  if (print_help == 1)
    print_usage();

  if (f1.empty())
    printErrorAndDie("--f1 argument required");
  if (f2.empty())
    printErrorAndDie("--f2 argument required");
  if (out.empty())
    printErrorAndDie("--out argument required");
  if (log.empty())
    printErrorAndDie("--log argument required");
  if (!string_ends_with(f1, ".gz"))
    printErrorAndDie("Argument to --f1 must be a bgzipped FASTQ file (and end in .gz)");
  if (!string_ends_with(f2, ".gz"))
    printErrorAndDie("Argument to --f2 must be a bgzipped FASTQ file (and end in .gz)");
  if (!file_exists(f1))
    printErrorAndDie("Argument to --f1 is not a valid file path");
  if (!file_exists(f2))
    printErrorAndDie("Argument to --f2 is not a valid file path");

  std::ofstream log_stream;
  log_stream.open(log, std::ofstream::out);
  if (!log_stream.is_open())
    printErrorAndDie("Failed to open the log file: " + log);
  
  ReadStitcher stitcher(max_read_len, max_k, min_bp_overlap, min_frac_correct);
  std::vector<std::string> l_reads;
  std::vector<std::string> r_reads;
  std::vector<int>         l_start;
  std::vector<int>         r_start;
  
  //100,100

  /*  
  generateReads(10000, max_read_len, l_reads, r_reads, l_start, r_start, 0.01);
  for(int i = 0; i < l_reads.size(); i++){
    int best_frac_idx = stitcher.stitch_reads(l_reads[i], r_reads[i]);
    if (best_frac_idx != -1){
      printStitching(l_reads[i], r_reads[i], best_frac_idx);
      std::cout << std::endl;
    }
  }
  */

  stitcher.stitch_fastq(f1, f2, out, log_stream);
  stitcher.print_base_qual_stats(log_stream);
  log_stream.close();
}

#ifndef READ_STITCHER_H
#define READ_STITCHER_H

#include <iostream>
#include <map>
#include <string>

#include "read_info.h"
#include "lca.h"

class ReadStitcher {
private:
  int    max_read_len;
  int    max_k;
  int    min_bp_overlap;
  double min_frac_correct;
  LCA *  lca;

  std::map<char, int64_t> match_base_quals_;
  std::map<char, int64_t> mismatch_base_quals_;

  void printStitching(const std::string& s1, const std::string& s2, int index);
  ReadInfo merge_read_information(ReadInfo& r1, ReadInfo& r2, int stitch_index, int num_bp_overlap, int num_mismatches);

public:
  ReadStitcher(int max_read_len, int max_k, int min_bp_overlap, double min_frac_correct);
  ~ReadStitcher();

  int stitch_reads(const std::string& s1, const std::string& s2, int& num_bp_overlap, int& num_mismatches);
  void stitch_fastq(std::string fastq_f1, std::string fastq_f2, std::string output_prefix);
  void kMismatch(const std::string& s1, const std::string& s2, int* best_frac_idx, double* best_frac, int& num_bp_overlap, int& num_mismatches);
  void print_base_qual_stats(std::ostream& out);
};

#endif

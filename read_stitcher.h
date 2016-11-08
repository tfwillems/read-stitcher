#ifndef READ_STITCHER_H
#define READ_STITCHER_H

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

  void printStitching(const std::string& s1, const std::string& s2, int index);
  ReadInfo merge_read_information(ReadInfo& r1, ReadInfo& r2, int stitch_index);

public:
  ReadStitcher(int max_read_len, int max_k, int min_bp_overlap, double min_frac_correct);
  ~ReadStitcher();

  int stitch_reads(const std::string& s1, const std::string& s2);
  void stitch_fastq(std::string fastq_f1, std::string fastq_f2, std::string output_prefix);
  void kMismatch(const std::string& s1, const std::string& s2, int* max_match_idx, int* max_matches, int* best_frac_idx, double* best_frac);
};

#endif

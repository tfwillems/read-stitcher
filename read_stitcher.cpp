#include <iostream>
#include <fstream>
#include <sstream>

#include "error.h"
#include "fastq_reader.h"
#include "fastq_writer.h"
#include "read_stitcher.h"
#include "stringops.h"

ReadStitcher::ReadStitcher(int max_read_len, int max_k, int min_bp_overlap, double min_frac_correct){
  this->max_read_len      = max_read_len;
  this->max_k             = max_k;
  this->min_bp_overlap    = min_bp_overlap;
  this->min_frac_correct  = min_frac_correct;
  lca = new LCA(2*(2*max_read_len+2)); // +2 due to separator character and terminating character
}

ReadStitcher::~ReadStitcher(){
  delete lca;
}

void ReadStitcher::printStitching(const std::string& s1, const std::string& s2, int index){
  std::cout << s1 << std::endl;
  std::string spacing = "";
  for (int i = 0; i < index; i++)
    spacing += " ";
  std::cout << spacing << s2 << std::endl;
}

void ReadStitcher::kMismatch(const std::string& s1, const std::string& s2, int* best_frac_idx, double* best_frac, int& num_bp_overlap, int& num_mismatches){
  std::string exp = s1 + '#' + s2;
  char separator  = '#';
  SuffixTree tree(exp);

  lca->processTree(tree);
  *best_frac      = 0;
  *best_frac_idx  = -1;
  num_bp_overlap  = -1;
  num_mismatches  = -1;
  //*max_matches    = -1;
  //*max_match_idx  = -1;

  for (int i = 0; i < ((int)s1.size()) - min_bp_overlap; i++){
    int sfx_offset_1 = 0;
    int sfx_offset_2 = 0;

    int k;
    for (k = 0; k < max_k; k++){
      // +1 due the way in which we increment the offsets (as we assume the match is followed by a mismatch)
      if (sfx_offset_2 == 1+s2.size())
	break;

      // +1 due to separator character
      int nmatch = lca->longestPrefix(tree, i+sfx_offset_1, s1.size()+1+sfx_offset_2);

      if (exp[i+sfx_offset_1+nmatch] == separator){
	sfx_offset_1 += nmatch;
	sfx_offset_2 += nmatch;
	break;
      }

      sfx_offset_1 += nmatch+1;
      sfx_offset_2 += nmatch+1;
    }


    // Check if stitching was successful
    if (i+sfx_offset_1 == s1.size() || (s2.size() >= min_bp_overlap && sfx_offset_2 == 1+s2.size())){
      double frac = 1.0*(sfx_offset_1-k)/sfx_offset_1;

      // Check if stitching satisfies minimum fraction requirement
      if (frac > min_frac_correct){
	if (frac > *best_frac){
	  *best_frac     = frac;
	  *best_frac_idx = i;
	  num_bp_overlap = sfx_offset_1;
	  num_mismatches = k;

	  // We've found a perfect match that meets all of the requirements. Because we're scanning from left to right,
	  // this is also the maximal match, so we're certain this is the best result and can abort the search
	  if (k == 0)
	    return;
	}

	/*
	if (sfx_offset_1-k > *max_matches){
	  *max_matches   = sfx_offset_1-k;
	  *max_match_idx = i;
	}
	*/
      }
    }
  }
}

ReadInfo ReadStitcher::merge_read_information(ReadInfo& r1, ReadInfo& r2, int stitch_index, int num_bp_overlap, int num_mismatches){
  std::string sequence = r1.get_sequence().substr(0, stitch_index);
  std::string quality  = r1.get_quality().substr(0, stitch_index);

  // Leading portion of stitched read
  std::string s1 = r1.get_sequence().substr(stitch_index);
  std::string s2 = r2.get_sequence();
  std::string q1 = r1.get_quality().substr(stitch_index);
  std::string q2 = r2.get_quality();

  // For the overlapping portion of the reads, select the base
  // with the highest quality score
  int i;
  for (i = 0; i < std::min(q1.length(), q2.length()); i++){
    if (q1[i] >= q2[i]){
      sequence += s1[i];
      quality  += q1[i];
    }
    else {
      sequence += s2[i];
      quality  += q2[i];
    }

    if (s1[i] != s2[i])
      mismatch_base_quals_[std::min(q1[i], q2[i])]++;
    else
      match_base_quals_[std::min(q1[i], q2[i])]++;
  }

  // Trailing portion of stitched read
  if (i < q1.length()){
    sequence += s1.substr(i);
    quality  += q1.substr(i);
  }
  else {
    sequence += s2.substr(i);
    quality  += q2.substr(i);
  }

  return ReadInfo("STITCHED_" + std::to_string(num_bp_overlap) + "_" + std::to_string(num_mismatches) + "_" + r1.get_identifier() , sequence, quality, false);
}

int ReadStitcher::stitch_reads(const std::string& s1, const std::string& s2, int& num_bp_overlap, int& num_mismatches){
  int best_frac_idx;
  double best_frac;

  kMismatch(s1, s2, &best_frac_idx, &best_frac, num_bp_overlap, num_mismatches);
  if (best_frac_idx != -1)
    return best_frac_idx;
  else
    return -1;
}

void ReadStitcher::stitch_fastq(std::string fastq_f1, std::string fastq_f2, std::string output_prefix, std::ostream& log){
  FASTQReader f1_reader(fastq_f1, true, false);
  FASTQReader f2_reader(fastq_f2, true, true);
  ReadInfo f1_read       = f1_reader.next_read();
  ReadInfo f2_read       = f2_reader.next_read();
  std::string prev_f1_id = f1_read.get_identifier();
  std::string prev_f2_id = f2_read.get_identifier();

  FASTQWriter f1_writer(output_prefix + "_1.fq.gz");
  FASTQWriter f2_writer(output_prefix + "_2.fq.gz");
  FASTQWriter stitched(output_prefix  + "_stitched.fq.gz");
  int32_t N_skip_count = 0, length_skip_count = 0, success_count = 0, fail_count = 0;

  while (true){
    if (f1_reader.is_empty())
      break;
    if (f2_reader.is_empty())
      break;

    f1_read = f1_reader.next_read();
    f2_read = f2_reader.next_read();
    if (f1_read.get_identifier().compare(f2_read.get_identifier()) != 0){
      std::stringstream error;
      error << "Mismatched read ids in FASTQ files:" << "\n"
	    << "\t" << f1_read.get_identifier() << " and " << f2_read.get_identifier();
      printErrorAndDie(error.str());
    }

    // Remove N's on ends of reads and low quality flanks
    f1_read.trimNTails();
    f2_read.trimNTails();
    char min_qual = '5';
    f1_read.trimLowQualityEnds(min_qual);
    f2_read.trimLowQualityEnds(min_qual);
    if (f1_read.empty() || f2_read.empty()){
      fail_count++;
      continue;
    }

    // Skip reads that exceed the max length, as they'll break the LCA computation
    if (f1_read.get_sequence().size() > max_read_len || f2_read.get_sequence().size() > max_read_len){
      f1_writer.write_read(f1_read);
      f2_writer.write_read(f2_read);
      length_skip_count++;
      continue;
    }

    // Attempt to stitch the reads together
    int best_frac_idx;
    double best_frac;
    int num_bp_overlap, num_mismatches;

    // Skip reads with N's, as the suffix tree doesn't accommodate it
    if (f1_read.get_sequence().find("N") != std::string::npos ||
	f2_read.get_sequence().find("N") != std::string::npos){
      N_skip_count++;
      continue;
    }

    kMismatch(f1_read.get_sequence(), f2_read.get_sequence(), &best_frac_idx, &best_frac, num_bp_overlap, num_mismatches);
    if (best_frac_idx != -1){
      // Stitching met requirements
      //std::cout << num_bp_overlap << " " << num_mismatches << std::endl;
      //printStitching(f1_read.get_sequence(), f2_read.get_sequence(), best_frac_idx);
      ReadInfo stitched_read = merge_read_information(f1_read, f2_read, best_frac_idx, num_bp_overlap, num_mismatches);
      stitched.write_read(stitched_read);
      success_count++;
    }
    else {
      // Retry stitching, reversing which read we assume comes upstream
      kMismatch(f2_read.get_sequence(), f1_read.get_sequence(), &best_frac_idx, &best_frac, num_bp_overlap, num_mismatches);
      if (best_frac_idx != -1){
	// Stitching met requirements
	//printStitching(f2_read.get_sequence(), f1_read.get_sequence(), best_frac_idx);
	ReadInfo stitched_read = merge_read_information(f2_read, f1_read, best_frac_idx, num_bp_overlap, num_mismatches);
	stitched.write_read(stitched_read);
	success_count++;
      }
      else {
	// Stitching did not meet requirements
	f1_writer.write_read(f1_read);
	f2_writer.write_read(f2_read);
	fail_count++;
      }
    }
  }

  if (length_skip_count != 0)
    log << "Skipped " << length_skip_count << " reads whose length was greater than " << max_read_len << "\n"
	      << "\t" << "If this is a significant fraction of your dataset, consider increasing --max-read-length" << std::endl;
  if (N_skip_count != 0)
    log << "Skipped " << N_skip_count << " reads with N bases" << std::endl;
  log << "Stitching succeeded for " << success_count << " out of " << (success_count+fail_count) << " remaining pairs of reads (" << (100.0*success_count/(success_count+fail_count)) << "%)" << std::endl;

  f1_reader.close();
  f2_reader.close();
  f1_writer.close();
  f2_writer.close();
  stitched.close();
}


void ReadStitcher::print_base_qual_stats(std::ostream& out){
  int64_t match_total = 0;
  for (auto count_iter = match_base_quals_.begin(); count_iter != match_base_quals_.end(); count_iter++)
      match_total += count_iter->second;
  for (auto count_iter = match_base_quals_.begin(); count_iter != match_base_quals_.end(); count_iter++)
    out << count_iter->first << "\t" << count_iter->second << "\t" << 100.0*count_iter->second/match_total << "\n";
  out << "\n";

  int64_t mismatch_total = 0;
  for (auto count_iter = mismatch_base_quals_.begin(); count_iter != mismatch_base_quals_.end(); count_iter++)
    mismatch_total += count_iter->second;
  for (auto count_iter = mismatch_base_quals_.begin(); count_iter != mismatch_base_quals_.end(); count_iter++)
    out << count_iter->first << "\t" << count_iter->second << "\t" << 100.0*count_iter->second/mismatch_total << "\n";
  out << "\n";
}

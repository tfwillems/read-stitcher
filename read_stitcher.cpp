#include <iostream>
#include <fstream>

#include "fastq_reader.h"
#include "fastq_writer.h"
#include "read_stitcher.h"

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

void ReadStitcher::kMismatch(const std::string& s1, const std::string& s2, int* max_match_idx, int* max_matches, int* best_frac_idx, double* best_frac){
  std::string  exp = s1 + '#' + s2;
  char separator = '#';
  SuffixTree tree(exp);

  lca->processTree(tree);
  *best_frac      = 0;
  *best_frac_idx  = -1;
  *max_matches    = -1;
  *max_match_idx  = -1;

  for(int i = std::max(0, ((int)s1.size())-((int)s2.size())); i < ((int)s1.size()) - min_bp_overlap; i++){
    int sfx_offset_1 = 0;
    int sfx_offset_2 = 0;

    int k;
    for(k = 0; k < max_k; k++){
      //+1 due to separator character
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
    if (i+sfx_offset_1 == s1.size()){
      double frac = 1.0*(sfx_offset_1-k)/sfx_offset_1;

      // Check if stiching satisfies minimum fraction requirement
      if (frac > min_frac_correct){
	if (frac > *best_frac){
	  *best_frac     = frac;
	  *best_frac_idx = i;
	}

	if (sfx_offset_1-k > *max_matches){
	  *max_matches   = sfx_offset_1-k;
	  *max_match_idx = i;
	}
      }
    }
  }
}

ReadInfo ReadStitcher::merge_read_information(ReadInfo& r1, ReadInfo& r2, int stitch_index){
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

  return ReadInfo(r1.get_identifier(), sequence, quality);
}

int ReadStitcher::stitch_reads(const std::string& s1, const std::string& s2){
  int max_match_idx;
  int max_matches;
  int best_frac_idx;
  double best_frac;

  kMismatch(s1, s2, &max_match_idx, &max_matches, &best_frac_idx, &best_frac);
  if (best_frac_idx != -1)
    return best_frac_idx;
  else
    return -1;
}

void ReadStitcher::stitch_fastq(std::string fastq_f1, std::string fastq_f2, std::string output_prefix){
  FASTQReader f1_reader(fastq_f1, true);
  FASTQReader f2_reader(fastq_f2, true);
  ReadInfo f1_read      = f1_reader.next_read();
  ReadInfo f2_read      = f2_reader.next_read();
  std::string prev_f1_id = f1_read.get_identifier();
  std::string prev_f2_id = f2_read.get_identifier();

  FASTQWriter f1_writer(output_prefix + "_1.fastq");
  FASTQWriter f2_writer(output_prefix + "_2.fastq");
  FASTQWriter stitched(output_prefix  + "_stitched.fastq");
  bool read_f1 = false;
  bool read_f2 = false;

  while (true){
    if (read_f1){
      if (f1_reader.is_empty())
	break;
      f1_read = f1_reader.next_read();
      if (f1_read.get_identifier().compare(prev_f1_id) <= 0){
	std::cerr << "ERROR: FIle provided to read stitcher must be sorted by read identifier. Exiting..." << std::endl;
	exit(1);
      }
      prev_f1_id = f1_read.get_identifier();
      read_f1    = false;
    }

    if (read_f2){
      if (f2_reader.is_empty())
	break;
      f2_read = f2_reader.next_read();
      if (f2_read.get_identifier().compare(prev_f2_id) <= 0){
	std::cerr << "ERROR: FIle provided to read stitcher must be sorted by read identifier. Exiting..." << std::endl;
	exit(1);
      }
      prev_f2_id = f2_read.get_identifier();
      read_f2    = false;
    }

    int comp = f1_read.get_identifier().compare(f2_read.get_identifier());
    if (comp == 0){
      // Attempt to stitch the reads together
      int max_match_idx;
      int max_matches;
      int best_frac_idx;
      double best_frac;

      kMismatch(f1_read.get_sequence(), f2_read.get_sequence(), &max_match_idx, &max_matches, &best_frac_idx, &best_frac);
      if (best_frac_idx != -1){
	// Stitching met requirements
	printStitching(f1_read.get_sequence(), f2_read.get_sequence(), best_frac_idx);
	ReadInfo stitched_read = merge_read_information(f1_read, f2_read, best_frac_idx);
	stitched.write_read(stitched_read);
      }
      else {
	// Stitching did not meet requirements
	f1_writer.write_read(f1_read);
	f2_writer.write_read(f2_read);
      }
      read_f1 = true;
      read_f2 = true;
    }
    else if (comp < 0){
      read_f1 = true;
      f1_writer.write_read(f1_read);
    }
    else {
      read_f2 = true;
      f2_writer.write_read(f2_read);
    }	
  }

  // Output any remaining reads in f1 without mates
  if (!read_f1)
    f1_writer.write_read(f1_read);
  while (!f1_reader.is_empty()){
    f1_read = f1_reader.next_read();
    f1_writer.write_read(f1_read);
  }

  // Output any remaining reads in f2 without mates
  if (!read_f2)
    f2_writer.write_read(f2_read);
  while (!f2_reader.is_empty()){
    f2_read = f2_reader.next_read();
    f2_writer.write_read(f2_read);
  }

  f1_reader.close();
  f2_reader.close();
  f1_writer.close();
  f2_writer.close();
  stitched.close();
}

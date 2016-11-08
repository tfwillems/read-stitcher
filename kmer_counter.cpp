#include <algorithm>
#include <iostream>
#include <vector>

#include <math.h>

#include "kmer_counter.h"

int KmerCounter::string_to_int(std::string s){
  int val = 0;
  for(int i = 0; i < s.length(); i++){
    val <<= num_bits;
    val  |= indices[s[i]];
  }
  return val;
}

std::string KmerCounter::int_to_string(int n){
  std::string res = "";
  for(int i = 0; i < k; i++){
    res += bases[(first_digit_mask & n) >> (k-1)*num_bits];
    n  <<= num_bits;
  }
  return res;
}

int KmerCounter::reverse_complement(int n){
  int comp = (~n) & all_digits_mask;
  int rc   = 0;
  for(int i = 0; i < k; i++){
    rc   <<= num_bits;
    rc    |= (comp & last_digit_mask);
    comp >>= num_bits; 
  }
  return rc;
}

int KmerCounter::min_permutation(int n){
  int min_val = n;
  for(int i = 0; i < k; i++){
    min_val = min_val < n ? min_val : n;
    n = (n & last_digit_mask) << (num_bits*(k-1)) | (n >> num_bits);
  }
  return min_val;
}

void KmerCounter::create_periodic_mappings(){
  min_mapping = new int[1 << (num_bits*k)];
  for (int i = 0; i < (1 << (num_bits*k)); i++)
    min_mapping[i] = std::min(min_permutation(i), min_permutation(reverse_complement(i)));
}

KmerCounter::KmerCounter(int k){
  this->k = k;
  if (k*num_bits > 32){
    std::cerr << "ERROR: Kmer counter unable to accomodate requested value of k as it is too large. Exiting... " << std::endl;
    exit(1);
  }

  // Construct the masks 
  last_digit_mask = 0;
  for(int i = 0; i < num_bits; i++){
    last_digit_mask <<= 1;
    last_digit_mask  |= 1;
  }
  all_digits_mask = 0;
  for(int i = 0; i < num_bits*k; i++){
    all_digits_mask <<= 1;
    all_digits_mask  |= 1;
  }
  first_digit_mask = (last_digit_mask << (k-1)*num_bits);

  // Map each DNA base to an integer
  // Indexing preserves base complement relationships
  for(int i = 0; i < 256; i++)
    indices[i] = -1;

  indices['A'] = 0;
  indices['C'] = 1;
  indices['G'] = 2;
  indices['T'] = 3;
  bases[0]     = 'A';
  bases[1]     = 'C';
  bases[2]     = 'G';
  bases[3]     = 'T';

  // Determine the minimum transformation for each kmer integer representation
  create_periodic_mappings();
}

std::unordered_map<int,int> KmerCounter::count_kmer_indexes(std::string s){
  if (s.length() < k)
    return std::unordered_map<int,int>();

  std::unordered_map<int,int> counts;
  int val = 0;
  for(int i = 0; i < k; i++){
    val <<= num_bits;
    if(indices[s[i]] == -1){
      std::cerr << "ERROR: Invalid character encountered in count_kmers. Exiting..." << std::endl;
      exit(1);
    }
    val  |= indices[s[i]]; 
  }
  counts[min_mapping[val]] = 1;

  for(int i = k; i < s.length(); i++){
    val = (val << num_bits) & all_digits_mask;
    if(indices[s[i]] == -1){
      std::cerr << "ERROR: Invalid character encountered in count_kmers. Exiting..." << std::endl;
      exit(1);
    }
    val  |= indices[s[i]];

    if (counts.find(min_mapping[val]) == counts.end())
      counts[min_mapping[val]]  = 1;
    else
      counts[min_mapping[val]] += 1;
  }
  return counts;
}

std::unordered_map<std::string,int> KmerCounter::count_kmer_words(std::string s){
  std::unordered_map<int,int> counts = count_kmer_indexes(s);
  std::unordered_map<std::string,int> res;

  for (std::unordered_map<int,int>::iterator it = counts.begin(); it != counts.end(); it++)
    res[int_to_string(it->first)] = it->second;
  return res;
}

double KmerCounter::calc_entropy(std::string s){
  std::unordered_map<int,int> counts = count_kmer_indexes(s);
  double entropy   = 0.0;
  double num_kmers = s.length() - k + 1;

  if(num_kmers < 1)
    return -1.0;

  for (std::unordered_map<int,int>::iterator it = counts.begin(); it != counts.end(); it++){
    double frac = it->second/num_kmers;
    entropy    += -frac*log2(frac);
  }
  return entropy;
}

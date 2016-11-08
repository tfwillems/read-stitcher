#ifndef KMER_COUNTER_H
#define KMER_COUNTER_H

#include <unordered_map>
#include <string>

class KmerCounter {
protected:
public:
  int indices[256];
  char bases[4]; 
  const static int num_bits  = 2;
  int first_digit_mask, last_digit_mask, all_digits_mask;

  int k;             // Length of kmer
  int* min_mapping;  // Mapping from a kmer's index to the index for its transformed value

  /* Returns the index associated with the kmer */
  int string_to_int(std::string s);

  /* Returns the string associated with a kmer's integer representation */
  std::string int_to_string(int n);

  /* Returns the index for the reverse complement of the kmer associated with the provided index */
  int reverse_complement(int n);

  /* Returns the index for the minimum permutation of the kmer associated with the provided index */
  int min_permutation(int n);

  /* Determines the transformed index associated with each possible kmer index */
  void create_periodic_mappings();

 public:
  KmerCounter(int k);

  ~KmerCounter(){ delete [] min_mapping;}

  /* 
   *  Returns an unordered_map, where the keys are the indexes of each kmer and the values 
   *  are the number of kmer occurrences within the provided string 
   */	 
  std::unordered_map<int,int> count_kmer_indexes(std::string s);

  /* 
   *  Returns an unordered_map, where the keys are the kmers and the values 
   *  are the number of kmer occurrences within the provided string 
   */	 
  std::unordered_map<std::string,int> count_kmer_words(std::string s);

  /* Returns the entropy of the kmers contained in s */
  double calc_entropy(std::string s);
};

#endif

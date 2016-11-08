#include <algorithm>

#include "error.h"
#include "stringops.h"

void reverse_complement(std::string& sequence){
  for (unsigned int i = 0; i < sequence.size(); i++){
    switch(sequence[i]){
    case 'A':
      sequence[i] = 'T';
      break;
    case 'C':
      sequence[i] = 'G';
      break;
    case 'G':
      sequence[i] = 'C';
      break;
    case 'T':
      sequence[i] = 'A';
      break;
    case 'N':
      sequence[i] = 'N';
      break;
    default:
      printErrorAndDie("Invalid character encountered in reverse_complement function");
      break;
    }
  }
  std::reverse(sequence.begin(), sequence.end());
}


#include "read_info.h"

void ReadInfo::trimNTails(){
  int start = 0, end = sequence.size()-1;
  while (start < sequence.size()){
    if (sequence[start] != 'N')
      break;
    start++;
  }

  while (end >= 0){
    if (sequence[end] != 'N')
      break;
    end--;
  }

  if (start <= end){
    std::string new_seq  = sequence.substr(start, end-start+1);
    std::string new_qual = quality.substr(start,  end-start+1);
    sequence = new_seq;
    quality  = new_qual;
  }
  else {
    sequence = "";
    quality  = "";
  }
}


void ReadInfo::trimLowQualityEnds(char min_qual){
  int start = 0, end = sequence.size()-1;
  while (start < quality.size()){
    if (quality[start] >= min_qual)
      break;
    start++;
  }

  while (end >= 0){
    if (quality[end] >= min_qual)
      break;
    end--;
  }

  if (start <= end){
    std::string new_seq  = sequence.substr(start, end-start+1);
    std::string new_qual = quality.substr(start,  end-start+1);
    sequence = new_seq;
    quality  = new_qual;
  }
  else {
    sequence = "";
    quality  = "";
  }
}

#include "read_info.h"

void ReadInfo::trimNTails(){
  int start = 0, end = sequence_.size()-1;
  while (start < sequence_.size()){
    if (sequence_[start] != 'N')
      break;
    start++;
  }

  while (end >= 0){
    if (sequence_[end] != 'N')
      break;
    end--;
  }

  if (start <= end){
    std::string new_seq  = sequence_.substr(start, end-start+1);
    std::string new_qual = quality_.substr(start,  end-start+1);
    ltrim_   += start;
    rtrim_   += (sequence_.size()-1-end);
    sequence_ = new_seq;
    quality_  = new_qual;
  }
  else {
    sequence_ = "";
    quality_  = "";
    ltrim_   += start;
  }
}

void ReadInfo::trimLowQualityEnds(char min_qual){
  int start = 0, end = sequence_.size()-1;
  while (start < quality_.size()){
    if (quality_[start] >= min_qual)
      break;
    start++;
  }

  while (end >= 0){
    if (quality_[end] >= min_qual)
      break;
    end--;
  }

  if (start <= end){
    std::string new_seq  = sequence_.substr(start, end-start+1);
    std::string new_qual = quality_.substr(start,  end-start+1);
    ltrim_   += start;
    rtrim_   += (sequence_.size()-1-end);
    sequence_ = new_seq;
    quality_  = new_qual;
  }
  else {
    sequence_ = "";
    quality_  = "";
    ltrim_   += start;
  }
}

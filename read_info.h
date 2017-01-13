#ifndef READ_INFO_H
#define READ_INFO_H

#include <assert.h>

#include <string>
#include <fstream>

class ReadInfo {
private:
  std::string identifier_;
  std::string sequence_;
  std::string quality_;
  bool rev_comp_;
  int ltrim_, rtrim_; // Amount of the sequence and quality scores that's been trimmed

public:
  ReadInfo(std::string identifier, std::string sequence, std::string quality, bool rev_complement){
    assert(sequence.size() == quality.size());
    identifier_ = identifier;
    sequence_   = sequence;
    quality_    = quality;
    rev_comp_   = rev_complement;
    ltrim_      = 0;
    rtrim_      = 0;
  }	

  const std::string& get_identifier(){ return identifier_; }
  const std::string& get_sequence()  { return sequence_;   }
  const std::string& get_quality()   { return quality_;    }
  bool reverse_complement()          { return rev_comp_;   }

  void trimNTails();

  void trimLowQualityEnds(char min_qual);

  bool empty(){
    return sequence_.size() == 0;
  }
};

#endif

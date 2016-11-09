#ifndef READ_INFO_H
#define READ_INFO_H

#include <assert.h>

#include <string>
#include <fstream>

class ReadInfo {
private:
  std::string identifier;
  std::string sequence;
  std::string quality;

public:
  ReadInfo(std::string identifier, std::string sequence, std::string quality){
    assert(sequence.size() == quality.size());
    this->identifier = identifier;
    this->sequence   = sequence;
    this->quality    = quality;
  }	

  const std::string& get_identifier(){ return this->identifier;}
  const std::string& get_sequence()  { return this->sequence;  }
  const std::string& get_quality()   { return this->quality;   }

  void trimNTails();

  void trimLowQualityEnds(char min_qual);

  bool empty(){
    return sequence.size() == 0;
  }
};

#endif

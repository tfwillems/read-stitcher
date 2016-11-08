#ifndef READ_INFO_H
#define READ_INFO_H

#include <string>
#include <fstream>

class ReadInfo {
private:
  std::string identifier;
  std::string sequence;
  std::string quality;

public:
  ReadInfo(std::string identifier, std::string sequence, std::string quality){
    this->identifier = identifier;
    this->sequence  = sequence;
    this->quality   = quality;
  }	

  const std::string& get_identifier(){ return this->identifier;}
  const std::string& get_sequence()  { return this->sequence;  }
  const std::string& get_quality()   { return this->quality;   }
};

#endif

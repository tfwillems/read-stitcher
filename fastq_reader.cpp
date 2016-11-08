#include <iostream>

#include "fastq_reader.h"

FASTQReader::FASTQReader(std::string filename, bool paired_end){
  this->filename   = filename;
  this->paired_end = paired_end;
  input.open(filename, std::ifstream::in);
}

FASTQReader::~FASTQReader(){ close(); }

bool FASTQReader::is_empty(){ return !input; }

ReadInfo FASTQReader::next_read(){
  if (!input){
    std::cerr << "ERROR: Attempt to read line in FASTQ_READER when stream is empty. Exiting..." << std::endl;
    exit(1);
  }
  std::string identifier, sequence, quality; 

  std::getline(input, identifier);
  if (!input){
    std::cerr << "ERROR: Attempt to read line in FASTQ_READER when stream is empty. Exiting..." << std::endl;
    exit(1);
  }

  std::getline(input, sequence);
  if (!input){
    std::cerr << "ERROR: Attempt to read line in FASTQ_READER when stream is empty. Exiting..." << std::endl;
    exit(1);
  }

  std::getline(input, quality);
  if (!input){
    std::cerr << "ERROR: Attempt to read line in FASTQ_READER when stream is empty. Exiting..." << std::endl;
    exit(1);
  }

  std::getline(input, quality);  
  if (identifier.at(0) != '@'){
    std::cerr << "ERROR: Read identifier is fastq file must begin with @ character. Exiting..." << std::endl;
    exit(1);
  }
	
  if (paired_end){
    if (identifier.length() > 2 && (identifier.back() == '1' || identifier.back() == '2') && identifier[identifier.length()-2] == '/') {
      identifier.pop_back();
      identifier.pop_back();
    }
    else {
      std::cerr << "ERROR: Read identifiers for paired-end files must end in /1 or /2. Exiting..." << std::endl;
      exit(1);
    }
  }
  return ReadInfo(identifier.substr(1), sequence, quality);
}

void FASTQReader::close(){
  if (input.is_open())
    input.close();
}

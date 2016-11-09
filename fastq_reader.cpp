#include <algorithm>
#include <iostream>

#include "error.h"
#include "fastq_reader.h"
#include "stringops.h"

FASTQReader::FASTQReader(std::string filename, bool paired_end, bool reverse_complement){
  this->filename   = filename;
  this->paired_end = paired_end;
  this->rev_complement = reverse_complement;
  input.open(filename.c_str());
  std::getline(input, next_line);
}

FASTQReader::~FASTQReader(){ close(); }

bool FASTQReader::is_empty(){ return !input; }

ReadInfo FASTQReader::next_read(){
  std::string identifier, sequence, quality;
  identifier = next_line;
  size_t space = identifier.find(" ");
  if (space != std::string::npos)
    identifier = identifier.substr(0, space);

  if (!input)
    printErrorAndDie("Attempt to read line in FASTQ_READER when stream is empty");
  std::getline(input, sequence);

  if (!input)
    printErrorAndDie("Attempt to read line in FASTQ_READER when stream is empty");
  std::getline(input, quality);

  if (!input)
    printErrorAndDie("Attempt to read line in FASTQ_READER when stream is empty");
  std::getline(input, quality);  

 if (identifier.at(0) != '@')
   printErrorAndDie("Read identifier is FASTQ file must begin with @ character");
	
  if (paired_end){
    if (identifier.length() > 2 && identifier[identifier.length()-2] == '/') {
      if (identifier.back() == '1' || identifier.back() == '2'){
	identifier.pop_back();
	identifier.pop_back();
      }
      else
	printErrorAndDie("Read identifiers for paired-end files must end in /1 or /2");
    }
  }

  if (rev_complement){
    // Reverse complement the sequence and reverse the quality scores
    reverse_complement(sequence);
    std::reverse(quality.begin(), quality.end());
  }

  std::getline(input, next_line);
  return ReadInfo(identifier.substr(1), sequence, quality, rev_complement);
}

void FASTQReader::close(){
  input.close();
}

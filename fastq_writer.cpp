#include "fastq_writer.h"

FASTQWriter::FASTQWriter(std::string filename){
  this->filename = filename;
  output.open(filename, std::ofstream::out);
}

FASTQWriter::~FASTQWriter(){
  close();
}

void FASTQWriter::close(){
  if (output.is_open())
    output.close();
}

void FASTQWriter::write_read(ReadInfo& read){
  std::string info = "@" + read.get_identifier() + "\n" + read.get_sequence() + "\n" + "+" + "\n" + read.get_quality() + "\n";
  output << info;
}

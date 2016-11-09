#include "fastq_writer.h"
#include "stringops.h"

FASTQWriter::FASTQWriter(std::string filename){
  this->filename = filename;
  output.open(filename.c_str());
}

FASTQWriter::~FASTQWriter(){
  close();
}

void FASTQWriter::close(){
  output.close();
}

void FASTQWriter::write_read(ReadInfo& read){
  std::string bases = read.get_sequence();
  std::string quals = read.get_quality();

  if (read.reverse_complement()){
    reverse_complement(bases);
    std::reverse(quals.begin(), quals.end());
  }

  output << "@"   << read.get_identifier() << "\n"
	 << bases << "\n"
	 << "+"   << "\n"
	 << quals << "\n";
}

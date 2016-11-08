#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <fstream>
#include <string>

#include "read_info.h"

class FASTQReader {
private:
  std::string filename;
  std::ifstream input;
  bool paired_end;

public:
  FASTQReader(std::string file, bool paired_end);
  ~FASTQReader();

  bool is_empty();
  ReadInfo next_read();
  void close();
};

#endif

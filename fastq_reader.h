#ifndef FASTQ_READER_H
#define FASTQ_READER_H

#include <fstream>
#include <string>

#include "read_info.h"
#include "bgzf_streams.h"

class FASTQReader {
private:
  std::string filename;
  bgzfistream input;
  bool paired_end;
  bool rev_complement;
  std::string next_line;

public:
  FASTQReader(std::string file, bool paired_end, bool reverse_complement);
  ~FASTQReader();

  bool is_empty();
  ReadInfo next_read();
  void close();
};

#endif

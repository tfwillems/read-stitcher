#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <fstream>

#include "bgzf_streams.h"
#include "read_info.h"

class FASTQWriter {
 private:
  std::string filename;
  bgzfostream output;

 public:
  FASTQWriter(std::string filename);
  ~FASTQWriter();

  void close();
  void write_read(ReadInfo& read);
};

#endif

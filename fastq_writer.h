#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include "read_info.h"

class FASTQWriter {
 private:
  std::string filename;
  std::ofstream output;

 public:
  FASTQWriter(std::string filename);
  ~FASTQWriter();

  void close();
  void write_read(ReadInfo& read);
};

#endif

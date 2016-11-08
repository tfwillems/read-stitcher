CXXFLAGS=-O3
CXX=gcc

SRCS=fastq_reader.cpp fastq_writer.cpp kmer_counter.cpp lca.cpp read_info.cpp read_stitcher.cpp suffix_tree.cpp main.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: ReadStitcher

ReadStitcher: $(OBJS)
	g++ -o ReadStitcher $(OBJS)

main.o: kmer_counter.h lca.h suffix_tree.h

read_stitcher.o: read_stitcher.h  lca.h
kmer_counter.o: kmer_counter.h
lca.o: lca.h suffix_tree.h
suffix_tree.o: suffix_tree.h 
read_info.o: read_info.h
fastq_reader.o: fastq_reader.h read_info.h
fastq_writer.o: fastq_writer.h read_info.h

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) main

##
## Makefile for all executables
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D__STDC_LIMIT_MACROS -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX -pthread -Wno-c++11-narrowing #-pedantic -Wunreachable-code -Weverything

## To build static executables, run:
##   make STATIC=1
## verify:
##   ReadStitcher
##
## To Create a static distribution file, run:
##   make static-dist
ifeq ($(STATIC),1)
LDFLAGS=-static
else
LDFLAGS=
endif

## Source code files, add new files to this list
SRC_COMMON  = error.cpp fastq_reader.cpp fastq_writer.cpp kmer_counter.cpp lca.cpp read_info.cpp read_stitcher.cpp stringops.cpp suffix_tree.cpp version.cpp
SRC_MAIN    = main.cpp

# For each CPP file, generate an object file
OBJ_COMMON  := $(SRC_COMMON:.cpp=.o)
OBJ_MAIN    := $(SRC_MAIN:.cpp=.o)

HTSLIB_ROOT=htslib
LIBS              = -L./ -lm -L$(HTSLIB_ROOT)/ -lz
HTSLIB_LIB        = $(HTSLIB_ROOT)/libhts.a

.PHONY: all
all: version ReadStitcher
	rm version.cpp
	touch version.cpp

# Create a tarball with static binaries
.PHONY: static-dist
static-dist:
	rm -f ReadStitcher
	$(MAKE) STATIC=1
	( VER="$$(git describe --abbrev=7 --dirty --always --tags)" ;\
	  DST="ReadStitcher-$${VER}-static-$$(uname -s)-$$(uname -m)" ; \
	  mkdir "$${DST}" && \
            cp ReadStitcher README.md "$${DST}" && \
            tar -czvf "$${DST}.tar.gz" "$${DST}" && \
            rm -r "$${DST}/" \
        )

version:
	git describe --abbrev=7 --dirty --always --tags | awk '{print "#include \"version.h\""; print "const std::string VERSION = \""$$0"\";"}' > version.cpp

# Clean the generated files of the main project only (leave Bamtools/vcflib alone)
.PHONY: clean
clean:
	rm -f *.o *.d ReadStitcher

# Clean all compiled files, including bamtools
.PHONY: clean-all
clean-all: clean

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

# The resulting binary executable

ReadStitcher: $(OBJ_COMMON) $(OBJ_MAIN) $(HTSLIB_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@

# Rebuild htslib library if needed
$(HTSLIB_LIB):
	git submodule update --init --recursive htslib
	git submodule update --recursive htslib
	cd htslib && $(MAKE)

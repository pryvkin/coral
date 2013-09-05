CXXFLAGS = -O3 -Wall

CXX = g++

# samtools
SAMTOOLS_DIR = samtools-0.1.18
SAMTOOLS_CFLAGS = -I$(SAMTOOLS_DIR)
SAMTOOLS_LDFLAGS = -L$(SAMTOOLS_DIR) -lbam -lz

# RSEQtools
RSEQTOOLS_DIR = $(shell pwd)/RSEQtools
BIOINFOCONFDIR = $(RSEQTOOLS_DIR)/bios/conf
BGRSEGMENTER = $(RSEQTOOLS_DIR)/mrf/bgrSegmenter
RSEQTOOLSBIOS = $(RSEQTOOLS_DIR)/bios/lib/libbios.a

# RSEQtools/bios depends on GSL
# bios: make && make prod
# bgrSegmenter: make bgrSegmenter

CXXFLAGS = -Wall -O3 $(SAMTOOLS_CFLAGS)
#CXXFLAGS = -Wall -g $(SAMTOOLS_CFLAGS)
LDFLAGS = $(SAMTOOLS_LDFLAGS)

PROGS = bin/annotate_smrna_loci bin/compute_genomic_lenvectors bin/index_genomic_lenvectors \
        bin/compute_locus_lenvectors
SRCS = $(wildcard *.cpp)
PROGS = $(patsubst %.cpp,bin/%,$(SRCS))
SCRIPTS = $(patsubst %,bin/%,$(wildcard *.sh))
SCRIPTS += $(patsubst %,bin/%,$(wildcard *.awk))
SCRIPTS += $(patsubst %,bin/%,$(wildcard *.rb))
SCRIPTS += $(patsubst %,bin/%,$(wildcard *.R))

all: Makefile bin $(RSEQTOOLSBIOS) $(BGRSEGMENTER) $(SAMTOOLS_DIR)/libbam.a $(PROGS) $(SCRIPTS)

bin:
	mkdir bin

bin/%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

bin/%.sh: %.sh
	cd bin && ln -s ../$< $<
bin/%.awk: %.awk
	cd bin && ln -s ../$< $<
bin/%.rb: %.rb
	cd bin && ln -s ../$< $<
bin/%.R: %.R
	cd bin && ln -s ../$< $<

$(SAMTOOLS_DIR)/libbam.a: 
	@echo "Building samtools..."
	cd $(SAMTOOLS_DIR) && make

$(RSEQTOOLSBIOS):
	@echo "Building RSEQtools/bios..."
	@export BIOINFOCONFDIR=$(BIOINFOCONFDIR); cd $(RSEQTOOLS_DIR)/bios && make && make prod

$(BGRSEGMENTER):
	@echo "Building RSEQtools/bgrSegmenter..."
	@export BIOINFOCONFDIR=$(BIOINFOCONFDIR); cd $(RSEQTOOLS_DIR)/mrf && make bgrSegmenter
	@ln -s $(RSEQTOOLS_DIR)/mrf/bgrSegmenter bin/bgrSegmenter

clean:
	rm -f bin/*
	cd $(SAMTOOLS_DIR) && make clean && cd -
	cd RSEQtools/bios && make clean && cd -
	cd RSEQtools/mrf && make clean && rm -f *.o && cd -



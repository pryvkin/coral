//  Copyright (c) 2013 University of Pennsylvania
//
//  Permission is hereby granted, free of charge, to any person obtaining a 
//  copy of this software and associated documentation files (the "Software"), 
//  to deal in the Software without restriction, including without limitation 
//  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//  and/or sell copies of the Software, and to permit persons to whom the 
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in 
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
//  DEALINGS IN THE SOFTWARE.


// compute length vectors (counts of reads of various lengths)
// at each position across the genome

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <deque>
#include "sam.h"

using namespace std;

struct SeqRead {
  float weight;
  int remaining_length;
  int length;
  SeqRead(float w, int l) : weight(w), remaining_length(l), length(l) { }
};

typedef deque< SeqRead > read_deque;

class WeightedReadLengthCoverageComputer {
  read_deque q;
  string chr;
  int front_pos;
  ostream &os;
  int min_length, max_length, length_range;
  char strand;

public:
  WeightedReadLengthCoverageComputer(const string &chr, ostream &o, int minlen, int maxlen, char str) :
    q(), chr(chr), front_pos(0), os(o), min_length(minlen), max_length(maxlen), length_range(1+maxlen-minlen),
    strand(str) { }

  void process_read(int start, int length, float weight);
  void finish();
};

void WeightedReadLengthCoverageComputer::process_read(int start, int length,
						float weight) {
  if (q.empty())
    front_pos = start;
  
  //cout << "\tfront pos: " << front_pos << "\n";
  //cout << "\tread: " << Read(start,length,weight) << "\n";
  
  // read at pos > front
  while(front_pos < start) {
    vector<float> length_sums(length_range, 0.0);
    float sum(0);

    for(read_deque::iterator p=q.begin(); p != q.end(); ++p) {
      if(p->length >= min_length && p->length <= max_length) {
	length_sums[p->length - min_length] += p->weight;   // accumulate wgt at this pos and length
	sum += p->weight;
      }

      --(p->remaining_length);
      if (p->remaining_length == 0)
	p->weight = 0;    // nullify weight when no length remaining
    }
    // output length vector for this pos
    if (sum > 0) {
      os << chr << "\t" << strand << "\t" << front_pos;
      for(int i=0; i < length_range; ++i)
	os << "\t" << length_sums[i];
      os << "\n";
    }
    
    ++front_pos;
  }
  
  // clean out dead entries
  while((!q.empty()) && q.front().weight == 0) {
    q.pop_front();
    //cout << "pop_front\n";
  }
  
  // add new entry

  q.push_back( SeqRead(weight, length) );
  //cout << "push_back (" << weight << "," << length << ")\n";  
}

void WeightedReadLengthCoverageComputer::finish() {
  if (!q.empty())
    process_read(front_pos + q.back().remaining_length, 1, 1);
}


int main(int argc, char **argv) {

  if (argc < 6) {
    cerr << "USAGE: " << argv[0] << " in.bam outplus outminus min_len max_len\n";
    return 1;
  }

  int min_length(atoi(argv[4]));
  int max_length(atoi(argv[5]));

  // open BAM file
  bamFile fp;
  bam_header_t *hdr;
  if ((fp = bam_open(argv[1], "r")) == 0) {
    cerr << "Failed to open BAM file " << argv[1] << "\n";
    return 1;
  }

  // open output files
  ofstream out_files[2];   // [plus, minus]
  for(int i=0; i < 2; ++i ) {
    out_files[i].open(argv[2+i]);
    if (!out_files[i].is_open() ) {
      cerr << "Failed to open file for writing: " << argv[2+i] << "\n";
      return 1;
    }
  }

  hdr = bam_header_read(fp);

  bam1_t *b = bam_init1();

  WeightedReadLengthCoverageComputer *pwc[2];  // [plus, minus]
  pwc[0] = NULL;
  pwc[1] = NULL;

  int prev_ref(-1);   // BAM ID of reference sequence (chromosome)

  string chr_name;
  char strand_name[2] = {'+', '-'};

  while(bam_read1(fp, b) > 0) {

    int curr_ref = b->core.tid;
    int strand = ((b->core.flag & 0x0010) > 0);
    int read_start = b->core.pos;
    int read_len = b->core.l_qseq;

    //uint8_t *aux_data = bam_aux_get(b, "XW");
    //float read_weight = bam_aux2f( aux_data );
    float read_weight = 1.0/float( bam_aux2i(bam_aux_get(b, "NH")));

    if (curr_ref != prev_ref) {
      // finished chromosome, so dump the queue and close the files
      if (prev_ref >= 0) {
	for(int i=0; i < 2; ++i) {
	  pwc[i]->finish();
	  delete pwc[i];
	}
      }

      chr_name = hdr->target_name[curr_ref];
      cout << chr_name << "... ";
      cout.flush();
      for(int i=0; i < 2; ++i)
	pwc[i] = new WeightedReadLengthCoverageComputer(chr_name, out_files[i],
							min_length, max_length,
							strand_name[i]);
    }

    pwc[strand]->process_read( read_start, read_len, read_weight );

    prev_ref = curr_ref;
  }

  // finish up
  if (prev_ref >= 0) {
    for(int i=0; i < 2; ++i) {
      pwc[i]->finish();
      delete pwc[i];
      out_files[i].close();
    }
  }

  cout << "\n";

  bam_destroy1(b);
  bam_header_destroy(hdr);
  bam_close(fp);

  return 0;
}

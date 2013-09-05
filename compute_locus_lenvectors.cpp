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
// for each locus in the genome

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include <cmath>

using namespace std;

bool verbose = true;

// a vector of read length counts at a genomic bp
struct LengthVector {
  string chr, strand;
  size_t bp;
  vector<int> data;
};

// parse a line from a genomic lenvector file
//   load_data=false to only read chr,strand (for indexing)
void parse_lenvec_line(const string &line, LengthVector& result, 
		      bool load_data=true) {
  result.data.clear();
  istringstream line_str(line);
  getline(line_str, result.chr, '\t');
  getline(line_str, result.strand, '\t');
  string bp_s;
  getline(line_str, bp_s, '\t');
  result.bp = atol(bp_s.c_str());
  if (load_data) {
    string s;
    while(getline(line_str, s, '\t'))
      result.data.push_back(atol(s.c_str()));
  }
}

struct BEDEntry {
  string chr;
  size_t start, end;
  string name;
  int score;
  string strand;
};
void parse_bed_line(const string &line, BEDEntry &result) {
  istringstream line_str(line);
  getline(line_str, result.chr, '\t');
  string s;
  getline(line_str, s, '\t');
  result.start = atol(s.c_str());
  getline(line_str, s, '\t');
  result.end = atol(s.c_str());
  getline(line_str, result.name, '\t');
  getline(line_str, s, '\t');
  // NOTE: assumes score is integer
  result.score = atoi(s.c_str());
  getline(line_str, result.strand, '\t');
}

int main(int argc, char **argv) {
  if (argc < 4) {
    cerr << "USAGE: " << argv[0]
	 << " locus_bed lenvector_prefix min_read_len\n";
    return(1);
  }

  string bed_fn(argv[1]);
  // length vector filenames
  string lv_fn[2] = {string(argv[2]) + ".plus",
		     string(argv[2]) + ".minus"};
  string lvidx_fn[2] = {lv_fn[0] + ".idx",
			lv_fn[1] + ".idx"};

  int min_read_len(atoi(argv[3]));

  ifstream bed_file(bed_fn.c_str());
  if (!bed_file.is_open()) {
    cerr << "Could not open BED file " << bed_fn << "\n";
    return(1);
  }
  
  // open lenvec files
  ifstream lv_files[2];
  for(int i=0; i < 2; ++i) {
    lv_files[i].open(lv_fn[i].c_str());
    if (!lv_files[i].is_open()) {
      cerr << "Could not open lenvector file " << lv_fn[i] << "\n";
      return(1);
    }
  }
  // open lenvec chromosome indices (lets us seek to a chromsome in O(1) time)
  ifstream lvidx_files[2];
  for(int i=0; i < 2; ++i) {
    lvidx_files[i].open(lvidx_fn[i].c_str());
    if (!lvidx_files[i].is_open()) {
      cerr << "Could not find index file " << lvidx_fn[i] << "\n";
      return(1);
    }
  }

  // build chromosome indices with key "chr;strand"
  // note: this is a std::map but accesses will be infrequent
  map<string, size_t> chr_strand_idx;
  string line;
  for (int i=0; i < 2; ++i) {
    char strand = ( (i == 0) ? '+' : '-' );
    while(getline(lvidx_files[i], line)) {
      istringstream line_str(line);
      string key;
      getline(line_str, key, '\t');
      key += string(";") + strand;
      
      string field;
      getline(line_str, field, '\t');
      size_t pos = atol(field.c_str());
      
      chr_strand_idx[key] = pos;
    }
  }

  // the intrachromsomal index maps genomic bp to byte in the lenvector file
  // we build this at runtime
  map<size_t, size_t> intra_chr_idx;
  // index every N bp;
  // lower = faster lookups+more mem usage
  size_t intra_chr_idx_stride = 100;
  
  size_t n_lengths = 0; // size of length vectors

  // this assumes chromsome_strands are contiguous
  // in the BED file; if not, it will be slower
  string prev_chr_strand;
  while(getline(bed_file, line)) {
    BEDEntry bed;
    parse_bed_line(line, bed);

    ifstream &lv_file = (bed.strand == "+") ? lv_files[0] : lv_files[1];
    
    string chr_strand(bed.chr + ";" + bed.strand);
    if (chr_strand != prev_chr_strand) {
      // we've encountered a new chr_strand, so build the intra_chr index
      if (verbose)
	cerr << chr_strand << "... ";

      intra_chr_idx.clear();
      if (chr_strand_idx.find(chr_strand) == chr_strand_idx.end()) {
	cerr << "Could not find start position for " << chr_strand << " in index!\n";
	return(1);
      }
      size_t chr_strand_start = chr_strand_idx[chr_strand];

      lv_file.seekg(chr_strand_start);
      string lv_line;
      LengthVector lenvec;
      size_t prev_lv_pos = chr_strand_start;
      size_t prev_lv_bp=0;
      while(getline(lv_file, lv_line)) {
	parse_lenvec_line(lv_line, lenvec, false);
	// stop if we hit another chr_strand
	if ( (lenvec.chr + ";" + lenvec.strand) != chr_strand)
	  break;
	// record the curr position if the stride has elapsed
	if (prev_lv_bp == 0 ||
	    (lenvec.bp >= (prev_lv_bp+intra_chr_idx_stride)))
	  intra_chr_idx[lenvec.bp] = prev_lv_pos;
	prev_lv_bp = lenvec.bp;
	// track beginning of the following line
	prev_lv_pos = lv_file.tellg();
      }
      lv_file.clear(); // clear eof flag just in case
    }

    // intra_chr index is built, so use it
    // find the first position occurring before our query pos
    map<size_t, size_t>::iterator it = intra_chr_idx.upper_bound(bed.start);
    if (it != intra_chr_idx.begin())
      --it;
    lv_file.seekg(it->second);
    //if (verbose) {
    //  cerr << "Locus at " << bed.start << ": found "
    //	   << it->first << "/" << it->second << "\n";
    //}
    vector<double> lenvec_sum(n_lengths, 0);
    while(getline(lv_file, line)) {
      LengthVector lenvec;
      parse_lenvec_line(line, lenvec);
      if (n_lengths == 0) {
	n_lengths = lenvec.data.size();
	lenvec_sum.resize(n_lengths, 0);

	// print header
	cout << "name";
	for(size_t i=0; i < n_lengths; ++i)
	  cout << "\tE" << (min_read_len+i);
	cout << "\n";
      }

      // skip positions before the locus start
      if (lenvec.bp < bed.start)
	continue;
      // stop when we pass the locus' end
      if (lenvec.bp >= bed.end)
	break;

      // accumulate sum
      for(size_t i=0; i < n_lengths; ++i)
	lenvec_sum[i] += lenvec.data[i];
    } // for each lenvec line
    lv_file.clear(); // clear eof flag just in case

    // normalize sum by locus length
    double row_sum(0);
    for(size_t i=0; i < n_lengths; ++i) {
      lenvec_sum[i] /= double(bed.end-bed.start); // count per bp
      ++lenvec_sum[i];   // laplace smoothing
      row_sum += lenvec_sum[i];
    }

    double expected_prob = 1.0/double(lenvec_sum.size());

    // print output

    cout << bed.name;

    for(size_t i=0; i < n_lengths; ++i) {
      // convert to probability
      lenvec_sum[i] /= row_sum;
      // convert to log2 odds ratio versus expected probability
      lenvec_sum[i] = log2(lenvec_sum[i] / expected_prob);
      cout << "\t" << lenvec_sum[i];
    }
    cout << "\n";

    prev_chr_strand = chr_strand;
  } // for each bed line
  if (verbose)
    cerr << "\n";

  return(0);
}

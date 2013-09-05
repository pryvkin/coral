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

// index the length vector files by chromosome for faster access

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char **argv) {

  if (argc < 2) {
    cerr << "USAGE: " << argv[0] << " lenvector_file\n";
    return(1);
  }

  string in_fn(argv[1]);
  string out_fn(in_fn + ".idx");

  ifstream infile(in_fn.c_str());
  if (!infile.is_open()) {
    cerr << "Could not open input file " << in_fn << "\n";
    return(1);
  }
  ofstream outfile(out_fn.c_str());
  if (!outfile.is_open()) {
    cerr << "Could not open output file " << out_fn << "\n";
    return(1);
  }

  // assumes chromosomes are contiguous in the file
  string line;              // current line
  string prev_chr;          // previously read chromosome
  size_t prev_filepos = 0;  // position in file before getline()
  while(getline(infile, line)) {
    string chr;
    istringstream iss(line);
    getline(iss, chr, '\t');

    // output to the index on encountering a new chromosmoe
    if (chr != prev_chr)
      outfile << chr << "\t" << prev_filepos << "\n";

    prev_filepos = infile.tellg();
    prev_chr = chr;
  }

  return 0;
}

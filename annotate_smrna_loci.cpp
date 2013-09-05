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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#ifdef __GNUC__
#include <ext/hash_map>
#else
#include <hash_map>
#endif

// allow hash_map to use std::string keys
namespace __gnu_cxx {
  template<> struct hash< std::string > {
    size_t operator()( const std::string& x ) const {
      return hash< const char* >()( x.c_str() );
    }
  };
}
namespace std {  using namespace __gnu_cxx; }
using namespace std;

/////////////////////////////

class Overlap {
public:
  string locus_id;
  string type;
  int amount;
  string desc;

  Overlap(const string& line) {
    istringstream iss(line);
    string tmp;
    getline(iss, locus_id, '\t');
    getline(iss, type, '\t');
    getline(iss, tmp, '\t');
    istringstream iss2(tmp);
    iss2 >> amount;
    getline(iss, desc, '\t');
  }
};

class OverlapComparator {
  hash_map<string,int> &clspri;
public:
  OverlapComparator(hash_map<string,int> &cp) : clspri(cp) { }
  bool operator()  (const Overlap &a, const Overlap &b) {
    return( clspri[a.type] < clspri[b.type] );
  }
};

hash_map<string,int> class_priorities;
OverlapComparator compare_overlaps_by_priority(class_priorities);

void process_buffer(const string& b) {
  string line;
  istringstream iss(b);
  vector<Overlap> overlaps;

  while(getline(iss, line))
    overlaps.push_back(Overlap(line));

  hash_map<string, int> types;
  vector<string> type_list;
  string basic_type;
  string basic_desc;
  for(size_t i=0; i < overlaps.size(); ++i) {
    if (types.find(overlaps[i].type) == types.end())
      types[ overlaps[i].type ] = 0;

    types[ overlaps[i].type ] = types[overlaps[i].type] + 1;
    if (overlaps[i].type != "intergenic")
      basic_type = overlaps[i].type;
    basic_desc += overlaps[i].desc;
    basic_desc += ",";
  }
  //cout << types.size() << "\n";

  if (types.size() > 2) {   // includes intergenic type
    basic_type = "multi";
  } else if (types.size() == 1)  {
    basic_type = "intergenic";
  }

  // Output: basic_class prioritized_class all_classes  basic_desc  prioritized_desc

  sort(overlaps.begin(), overlaps.end(), compare_overlaps_by_priority);
  string pri_type(overlaps[0].type);
  string pri_desc(overlaps[0].desc);

  cout << overlaps[0].locus_id << "\t" << basic_type << "\t" << pri_type << "\t";
  for( hash_map<string,int>::iterator it=types.begin(); it != types.end(); ++it) {
    for(int i=0; i < (it->second); ++i) 
      cout << (it->first) << ",";
  }
  cout << "\t" << basic_desc << "\t" << pri_desc << "\n";
}

/////////////////////

int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "USAGE: " << argv[0] << " in.intersect3 class_pri.txt\n";
    return 1;
  }

  ifstream intersect_file(argv[1]);
  if (!intersect_file.is_open()) {
    cerr << "Could not open intersect file " << argv[1] << "\n";
    return 1;
  }

  ifstream clspri_file(argv[2]);
  if (!clspri_file.is_open()) {
    cerr << "Could not open class priority file " << argv[2] << "\n";
    return 1;
  }

  string line;

  // read class priorities
  while(getline(clspri_file, line)) {
    istringstream iss(line);
    string k;
    int v;
    iss >> k;
    iss >> v;
    class_priorities[k] = v;
  }

  // process all overlaps
  string prev_loc_id("");
  string buffer;

  while(getline(intersect_file, line)) {
    string curr_loc_id;
    string dummy;

    istringstream iss(line);
    getline(iss, curr_loc_id, '\t');
    
    if ((!prev_loc_id.empty()) && (curr_loc_id != prev_loc_id)) {
      process_buffer(buffer);
      buffer.clear();
    }
    
    buffer += line;
    buffer += "\n";

    prev_loc_id = curr_loc_id;
  }

  process_buffer(buffer);
  return 0;

}

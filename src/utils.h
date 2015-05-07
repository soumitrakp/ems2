#ifndef __UTILS_
#define __UTILS_

#include<string>
#include<vector>
#include<list>
#include<algorithm>
#include<sstream>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<csignal>
#include<cstring>
#include<cstdint>
#include<cmath>
#include<ctime>
#include<cassert>


typedef std::vector<std::string> Reads;
typedef std::vector<std::string> Motifs;

const std::string domain="ACGT";
const size_t domain_size = domain.size();

std::string to_str(uint64_t id, int l) {
  std::string s;
  do {
    int c = id & 0x3;
    if (c==0) { s.insert(0,1,'A'); }
    else if (c==1) { s.insert(0,1,'C'); }
    else if (c==2) { s.insert(0,1,'G'); }
    else if (c==3) { s.insert(0,1,'T'); }
    id >>= 2; l--;
  } while (id>0);
  while (l>0) {l--; s.insert(0,1,'A');}
  return s;
}

uint64_t to_int(std::string s) {
  uint64_t idx = 0;
  int len = s.size();
  for(int i=0;i<len;i++) {
    switch(s[i]) {
      case 'A':
      case 'a':
        idx=idx<<2;
        break;
      case 'C':
      case 'c':
      idx=idx<<2;
      idx+=1;
        break;
      case 'G':
      case 'g':
        idx=idx<<2;
      idx+=2;
        break;
      case 'T':
      case 't':
      idx=idx<<2;
      idx+=3;
        break;
      default:
        break;
    }
  }
  return idx;
}


void read_file(const char *fname, Reads &reads) {
  std::ifstream f1(fname);
  if (!f1.is_open()) {
    std::cerr << "ERROR: could not open file " << fname << "\n";
    exit(-1);
  }
  std::string p1;
  while(std::getline(f1,p1)) {
    if ((p1.size() > 0) && (p1[0] != '>')) reads.push_back(p1);
  }
}



bool diff_motifs(Motifs & m1, Motifs & m2) {
  if (m1.size() != m2.size()) return true;
  std::sort(m1.begin(), m1.end());
  std::sort(m2.begin(), m2.end());
  for (size_t i=0; i<m1.size(); i++) {
    if (m1[i] != m2[i]) return true;
  }
  return false;
}

void write_to_file(int l, int d, Motifs &m, std::string fname) {
  std::ofstream myfile;
  myfile.open (fname);
  myfile << l << " " << d << std::endl;
  for (size_t i=0; i<m.size(); i++) {
    myfile << m[i] << std::endl;
  }
  myfile.close();
}

void write_to_file(Motifs &m, std::string fname) {
  std::ofstream myfile;
  std::sort(m.begin(), m.end());
  myfile.open (fname);
  for (size_t i=0; i<m.size(); i++) {
    myfile << m[i] << std::endl;
  }
  myfile.close();
}


int edist(std::string& s1, std::string& s2) {
  //std::cout << "in edit dist " << s1 << ", " << s2 << std::endl;
  const size_t len1 = s1.size(), len2 = s2.size();
  std::vector<std::vector<unsigned int> > d(len1 + 1, std::vector<unsigned int>(len2 + 1));

  d[0][0] = 0;
  for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
  for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

  for(unsigned int i = 1; i <= len1; ++i)
    for(unsigned int j = 1; j <= len2; ++j)
      d[i][j] = std::min( std::min(d[i - 1][j] + 1, d[i][j - 1] + 1),
          d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) );
  //std::cout << "in edit dist done" << std::endl;
  return d[len1][len2];
}

bool has_overlap(std::string x, std::string y, int common) {
  int i = x.size() - common;
  int j = 0;
  while (common>0) {
    if (x[i] != y[j]) return false;
    i++; j++; common--;
  }
  return true;
}

bool found_in_seq(std::string& candi, const std::string& seq, int l, int d) {
  int m = seq.size();
  for (int q=-d; q<=+d; q++) {
    int k = l+q;
    for (int i=0; i<m-k+1; ++i) {
      std::string x = seq.substr(i, k);
      if (edist(x, candi) <= d) {
        return true;
      }
    }
  }
  return false;
}

inline double diffclock(clock_t clock1,clock_t clock2)
{
  double diffticks=clock1-clock2;
  double diffsec=(diffticks)/CLOCKS_PER_SEC;
  return diffsec;
}

inline void show_progress(size_t done, size_t todo, clock_t start_clk) {
  double sec = diffclock(clock(), start_clk); 
  std::cout << "\rProcessed " << done << " of " << todo << " candidates (" << std::fixed << std::setprecision(1) << (100.0*done)/todo << "%) in " << sec << " seconds.     ";
  std::cout.flush();
}

std::string removeExtension(const std::string & filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
}

std::string get_out_file(const std::string& input, int l, int d, const std::string & prefix="") {
  std::ostringstream oss;
  oss << removeExtension(input) << "_" << prefix << "_l" << l << "_d" << d << ".txt";
  return oss.str();
}

std::string get_out_file(const std::string& input, int l, int d, int ldash, const std::string & prefix="") {
  std::ostringstream oss;
  oss << removeExtension(input) << "_" << prefix << "_l" << l << "_d" << d << "_ld" << ldash << ".txt";
  return oss.str();
}

template<typename T>
void printList(const std::vector<std::vector<T>* >& l, const std::string& msg="", int m=-1) {
    std::cout << "================ " << msg << " ==================" << std::endl;
    //if (m==-1) m = l.size();
    for (auto &v : l) {
      for (auto const &i : *v) {
        std::cout << std::setw(4) << (int)i;
      }
      std::cout << std::endl;
      //if (--m == 0) return;
    }
}

#endif // __UTILS_

#ifndef __EMS_HPP__
#define __EMS_HPP__

#include <sys/time.h>
#include <sys/resource.h>
#include "utils.h"

template <class DerivedMotifFinder>
class MotifFinder {

protected:
  std::string name;
  std::string input;
  Reads reads;
  int l, d;
  Motifs motifs;

public:

  MotifFinder(const std::string& _name, const Reads &_reads, int _l, int _d) : name(_name), reads(_reads), l(_l), d(_d) { }
  MotifFinder(const std::string& _name, const std::string &_input, int _l, int _d) : name(_name), input(_input), l(_l), d(_d) {
    read_file(input.c_str(), reads);
  }
  ~MotifFinder() { }

  void search() { std::cout << "Please implement search function in the derived class." << std::endl; }

  Motifs& searchGetMotifs() {
    static_cast<DerivedMotifFinder*>(this)->search();
    return motifs;
  }

  void searchWriteMotifs() {
    std::string output = get_out_file(input, l, d, name);
    std::cout << "l      = " << l << ", d = " << d << std::endl;
    std::cout << "input  = " << input << std::endl;
    std::cout << "output = " << output << std::endl;
    clock_t begin=clock();
    static_cast<DerivedMotifFinder*>(this)->search();
    clock_t end=clock();
    double elapsed = diffclock(end,begin);
    struct rusage usage;
    getrusage( RUSAGE_SELF, &usage );
    std::ofstream out(output);
    std::cout << "\r" << name << ": (" << l << "," << d << ") Edited Motifs found (in "<< elapsed << " sec, using " << (size_t)usage.ru_maxrss << " KB): " << motifs.size() << "       " << std::endl;
    for (size_t i=0; i<motifs.size(); ++i) {
      out << motifs[i] << std::endl;
    }
    out.close();
  }
};

#endif // __EMS_HPP__


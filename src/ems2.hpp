#ifndef __EMS2_H_
#define __EMS2_H_

#include "ems.hpp"
#include "motif_tree.hpp"

class Ems2 : public MotifFinder<Ems2> {
  MotifTree main_tree;
  MotifTree *curr_tree;
  size_t count = 0;
  int leftmost;
  int rightmost;
  std::string x;
  std::string mark;

  void gen_nbrhood3(int start, int alpha) {
    int len = x.size();
    if (alpha > 0) {
      for (int j=start; j<len+1; j++) {
        if ((mark[j] == 'Y') || (mark[j] == 'Z')) continue;
        x.insert(j,1,'*');
        mark.insert(j,1,'N');
        gen_nbrhood3(j+1, alpha-1); 
        x.erase(j,1);
        mark.erase(j,1);
      }
    } else {
      curr_tree->insert(x);
    }
  }

  void gen_nbrhood2(int start, bool stars_before, int sigma, int alpha) {
    int len = x.size();
    if (sigma > 0) {
      for (int j=start; j<len; j++, stars_before=false) {
        if ((mark[j] == 'Y')) continue;
        char rr = mark[j+1];
        if (!rightmost && stars_before && (rr == 'Y')) continue;
        char t = x[j];
        char r = mark[j];
        x[j] = '*';
        mark[j] = 'Y';
        if (!leftmost && stars_before && (rr == 'N')) mark[j+1] = 'Z';
        gen_nbrhood2(stars_before, j+1, sigma-1, alpha); 
        x[j] = t;
        mark[j] = r;
        mark[j+1] = rr;
      }
    } else {
      gen_nbrhood3(0, alpha);
    }
  }

  void gen_nbrhood(int start, int end, int delta, int sigma, int alpha) {
    if (delta > 0) {
      for (int j=start; j<end; j++) {
        char t = x[j];
        char rr = mark[j];
        x.erase(j,1);
        mark.erase(j,1);
        char r = mark[j];
        mark[j] = 'Y';
        gen_nbrhood(j, end-1, delta-1, sigma, alpha);
        x.insert(j,1,t);
        mark[j] = r;
        mark.insert(j,1, rr);
      }
    } else {
      gen_nbrhood2(0, true, sigma, alpha);
    }
  }

  void gen_all(std::string& seq, int l, int d) {
    int m = seq.size();
    for (int q=-d; q<=+d; q++) {
      int k = l+q;
      for (int delta = std::max(0,q); delta <= (d+q)/2; delta++) {
        int alpha = delta - q;
        int sigma = d - alpha - delta;
        for (int i=0; i<m-k+1; i++) {
          x = seq.substr(i, k);
          mark = std::string(k+1, 'N');
          int start;
          if (i>0) { mark[0] = 'Z'; leftmost = 0;} else {leftmost = 1;}
          if (i+k<m) { mark[k] = 'Z'; start = 1; rightmost = 0;} else { start = 0; rightmost = 1;}
          gen_nbrhood(start, k, delta, sigma, alpha);
        }
      }
    }
  }
  void reset() {
    for (size_t i=0; i< motifs.size(); i++) {
      std::cout << motifs[i] << std::endl;
    }
  }

  public:

  Ems2(const std::string &input, int l, int d):
    MotifFinder("Ems2", input, l, d), main_tree(l, motifs, "main"), curr_tree(&main_tree) {
    }

  void search() {
    std::cout << "Processing sequence " << 0 << "..." << std::endl;
    std::cout.flush();
    gen_all(reads[0], l, d);
    for (size_t i=1; i<reads.size(); i++) {
      std::cout << "Processing sequence " << i << "..." << std::endl;
      std::cout.flush();
      Motifs tmp;
      MotifTree tmp_tree(l, tmp, "tmp");
      curr_tree = &tmp_tree;
      gen_all(reads[i], l, d);
      main_tree.intersect(curr_tree);
    }
    main_tree.traverse();
  }

};

#endif // __EMS2_H_


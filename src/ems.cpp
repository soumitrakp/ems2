#include "ems.hpp"
#include "ems1.hpp"
#include "ems2.hpp"

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cout << "Usage: " << argv[0] << " <version> <input-sequence-file> <l> <d>" << std::endl;
    std::cout << "       versions: 1, 2" << std::endl;
    exit(-1);
  }
  std::string version(argv[1]);
  std::string input(argv[2]);
  int l = atoi(argv[3]);
  int d = atoi(argv[4]);

  if (version == "1") {
    Ems1 ems(input, l, d);
    ems.searchWriteMotifs();
  } else if (version == "2") {
    Ems2 ems(input, l, d);
    ems.searchWriteMotifs();
  } else {
    std::cout << "Unknown version: " << version << std::endl;
    exit(-1);
  }
}


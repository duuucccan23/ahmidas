#include <iostream>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <ildg_conf>\n";
    exit(EXIT_FAILURE);
  }
  std::string const filename = argv[1];
  int status = Tool::IO::checkILDG(filename);
  std::cout << filename << (status < 0 ? " BROKEN\n" : " CORRECT\n");
  return status;
}

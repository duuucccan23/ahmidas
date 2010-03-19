#include "FileReader.ih"

namespace Input
{

  FileReader::FileReader(std::string file)
  {
    std::ifstream fin(file.c_str());
    if (!fin.is_open())
    {
      std::cerr << "ERROR: File " << file << " not found" << std::endl;
      exit(1);
    }
    else if (!fin.good())
    {
      std::cerr << "ERROR: Something went wrong when trying to read " << file << std::endl;
      std::cerr << "Check file permissions!" << file << std::endl;
      exit(1);
    }

    Parser parser;
    std::string line;
    std::string ID;
    std::string contents;

    line.reserve(buf_size);

    char * buffer = new char[buf_size];

    while(fin.getline(buffer, buf_size))
    {
      line = buffer;
      parser.parseLine(line, ID, contents);
    }
    delete [] buffer;
    fin.close();
  }
}

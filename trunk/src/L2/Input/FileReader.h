#ifndef __GUARD_FILE_READER__
#define __GUARD_FILE_READER__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstring>
#include <stack>
#include <utility>

// SD: Actually, this seems to be incompatible, and since boost uses libxml as well,
// I think it is better not to use it at all and write an own simple parser.
// #include <libxml2/libxml/xmlmemory.h>
// #include <libxml2/libxml/parser.h>

#include <L2/Input/Parser.h>


namespace Input
{

  class FileReader
  {

    enum ValType
    {
      type_FLOAT = 1,
      type_INTEGER = 2,
      type_SIZE_T = 3,
      type_STRING = 4
    };

    enum InputTag
    {
      type_FILES = 100,
      type_PARAM = 200
    };

    static const size_t buf_size = 1000;


//     std::multimap< int, size_t >      mm_size_ts;
//     std::multimap< int, int >         mm_integers;
//     std::multimap< int, double >      mm_floats;
//     std::multimap< int, std::string > mm_strings;

    std::multimap< InputTag, std::pair< std::string, std::string > > input;

    std::vector< std::vector< std::string > > files;


    public:

    FileReader(std::string file);

    size_t getL();
    size_t getT();

    std::vector< std::vector< std::string > > const &getFiles();

    size_t * getPosition(int tag);

  };
}
#endif

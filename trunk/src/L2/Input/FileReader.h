#ifndef __FileReader_h__
#define __FileReader_h__

#include <string>
#include <vector>
#include <map>
#include <cstring>
#include <utility>

// SD: Actually, this seems to be incompatible, and since boost uses libxml as well,
// I think it is better not to use it at all and write an own simple parser.
// #include <libxml2/libxml/xmlmemory.h>
// #include <libxml2/libxml/parser.h>

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

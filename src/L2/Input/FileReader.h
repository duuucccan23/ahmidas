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

  enum FileTag
  {
    tag_NOFILE,
    tag_FILE
  };

  class FileReader;

  class File
  {
    friend class FileReader;

    std::string d_type;
    std::string d_directory;
    std::string d_filenameBase;
    std::string d_filenameEnding;
    int d_firstIndex;
    int d_lastIndex;
    size_t d_indexWidth;
    
    public:
    
      File(std::string const type, std::string const directory,
           std::string const filenameBase, std::string const filenameEnding,
           std::string const firstIndex, std::string const lastIndex, std::string const indexWidth);
      
  };

  class FileReader
  {
  


//     enum ValType
//     {
//       type_FLOAT = 1,
//       type_INTEGER = 2,
//       type_SIZE_T = 3,
//       type_STRING = 4
//     };
// 
//     enum InputTag
//     {
//       type_FILES = 100,
//       type_PARAM = 200
//     };

    static const size_t buf_size = 1000;


    std::stack< std::string > open_tags;

    std::multimap< std::string, std::string > input;
    
    // this one does not store the filenames but the raw data
    std::vector< File > files;
    
    void line_error(size_t const line, std::string const message = "");
    

    public:

    FileReader(std::string file);
    
    void initializeParameters(size_t &L, size_t &T, 
                         std::vector< std::vector< std::string > > &filenames, 
                         std::vector< double > &floats, 
                         std::vector< size_t * > &positions, 
                         std::vector< int > &operators) const;

  };
}

#include "FileReader/FileReader.inlines"

#endif

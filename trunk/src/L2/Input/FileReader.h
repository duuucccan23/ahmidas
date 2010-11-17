#pragma once

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
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
    size_t d_firstIndex;
    size_t d_lastIndex;
    size_t d_indexWidth;

    public:

      File(std::string const type, std::string const directory,
           std::string const filenameBase, std::string const filenameEnding,
           std::string const firstIndex, std::string const lastIndex, std::string const indexWidth);

  };

  class FileReader
  {

    static const size_t buf_size = 1000;


    std::stack< std::string > open_tags;

    std::multimap< std::string, std::string > input;

    // this one does not store the filenames but the raw data
    std::vector< File > files;

    void line_error(size_t const line, std::string const message = "");

    public:

    FileReader(std::string const file);
    ~FileReader();

    void initializeParameters(size_t &L, size_t &T,
                              std::vector< std::vector< std::string > > &filenames,
                              std::map< std::string, double > &floats,
                              std::vector< size_t * > &positions,
                              std::vector< int > &operators) const;

    // shorter version
    void initializeParameters(size_t &L, size_t &T,
                              std::vector< std::vector< std::string > > &filenames,
                              std::map< std::string, double > &floats) const;

    // longer version 
    void initializeParameters(size_t &L, size_t &T,
                              std::vector< std::vector< std::string > > &filenames,
                              std::map< std::string, double > &floats,
                              std::vector< size_t * > &positions,
                              std::vector< int > &operators,
                              std::vector< int > &rcombinations) const;

    // longer version with output filename
    void initializeParameters(size_t &L, size_t &T,
                              std::vector< std::vector< std::string > > &filenames,
                              std::map< std::string, double > &floats,
                              std::vector< size_t * > &positions,
                              std::vector< int > &operators,
                              std::vector< int > &rcombinations,
			      std::string &outputname) const;


  };
}

#include "FileReader/FileReader.inlines"

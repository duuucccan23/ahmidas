#ifndef __GUARD_PARSER__
#define __GUARD_PARSER__

#include <iostream>
#include <string>
#include <utility>

namespace Input
{
  enum TagState
  {
    tag_OPENING   =  -1,
    tag_COMPLETE  =   0,
    tag_CLOSING   =   1,
    tag_EMPTY     =   2,
    tag_UNMATCHED =   3,
    tag_NOCONTENT =   4,
    tag_ERROR_A   =  10,
    tag_ERROR_B   =  11,
    tag_COMMENT   = 100,
    tag_HEADLINE  = 101
  };

  class Parser
  {
    public:

      /* This function reads a single line containing one of the following:
          - opening and closing tag with some content (characters) in between
            (return value: tag_COMPLETE, tag contains identifyer and data contains the content)
          - opening tag (return value: tag_OPENING, tag contains identifyer, data is empty)
          - closing tag (return value: tag_CLOSING, tag contains identifyer, data is empty)
         If anything else occurs, it returns another TagState indicating the error
         and the strings tag and data are empty on return.
         An XML comment is simply ignored, the function returns tag_COMMENT.
      */
      int parseLine(std::string const &line, std::string& tag, std::string& data) const;
  };
}

#include "Parser/Parser.inlines"

#endif

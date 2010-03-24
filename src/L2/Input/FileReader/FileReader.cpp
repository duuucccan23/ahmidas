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
    int state(0), global_state(0);
    size_t line_counter(0);

    while (fin.getline(buffer, buf_size))
    {
      line_counter++;
      line = buffer;
      state = parser.parseLine(line, ID, contents);
      global_state += state;
      switch (state)
      { 
        case tag_OPENING:
          open_tags.push(ID);
          break;

        case tag_CLOSING:
          if(open_tags.top() == ID)
          {
            open_tags.pop();
          }
          else
          {
            line_error(line_counter, "Unmatched opening tag!");
          }
          break;

        case tag_COMPLETE:
          input.insert(std::make_pair(ID, contents));
          break;

        default:
          line_error(line_counter, "Unmatched opening tag!");
      }
      
    } // while (fin.getline(buffer, buf_size))
    
    if (global_state != 0)
    {
      std::cerr << "Unknown error occurred after file " << file << " had been read."  << std::endl;
      exit(100);
    }
    
    delete [] buffer;
    fin.close();
  }
  
  
  // there should be an error handling here
  void FileReader::initializeParameters(size_t &L, size_t &T, 
                                        std::vector< std::vector< std::string > > &filenames, 
                                        std::vector< double > &floats, 
                                        std::vector< size_t * > &positions, 
                                        std::vector< int > &operators) const
  {
    std::multimap< std::string, std::string >::iterator It; 
    for (It=input.begin(); It!=input.end(); It++)
    {
         
    }

    
  }                                   

}

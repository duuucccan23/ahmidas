#include "FileReader.ih"

namespace Input
{

  FileReader::FileReader(std::string const file)
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

    bool fileMode = false;
    std::string * fileData = NULL;

    while (fin.getline(buffer, buf_size))
    {
      line_counter++;
      line = buffer;
      state = parser.parseLine(line, ID, contents);
      switch (state)
      {
        case tag_OPENING:
          global_state += state;
          if (ID == "file")
          {
            fileMode = true;
            fileData = new std::string[7];
          }
          open_tags.push(ID);
          break;

        case tag_CLOSING:
          global_state += state;
          if(open_tags.top() == ID)
          {
            if (ID == "file")
            {
              files.push_back(File(fileData[0], fileData[1], fileData[2], fileData[3],
                                   fileData[4], fileData[5], fileData[6]));
              fileMode = false;
              delete [] fileData;
              fileData = NULL;
            }
            open_tags.pop();
          }
          else
          {
            line_error(line_counter, "Unmatched opening tag!");
          }
          break;

        case tag_COMPLETE:
          if (!fileMode)
          {
            input.insert(std::make_pair(ID, contents));
          }
          else
          {
            if (ID == "type")
            {
              fileData[0] = contents;
            }
            else if (ID == "directory")
            {
              fileData[1] = contents;
            }
            else if (ID == "filenameBase")
            {
              fileData[2] = contents;
            }
            else if (ID == "filenameEnding")
            {
              fileData[3] = contents;
            }
            else if (ID == "firstIndex")
            {
              fileData[4] = contents;
            }
            else if (ID == "lastIndex")
            {
              fileData[5] = contents;
            }
            else if (ID == "indexWidth")
            {
              fileData[6] = contents;
            }
          }
          break;

        case tag_COMMENT:
          break;
        case tag_HEADLINE:
          break;

        default:
          std::cerr << state << std::endl;
          line_error(line_counter, "Unmatched opening tag!");
      }// switch (state)

    } // while (fin.getline(buffer, buf_size))

    if (global_state != 0 || !open_tags.empty())
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
                                        std::map< std::string, double > &floats,
                                        std::vector< size_t * > &positions,
                                        std::map< std::string, int > &operators) const
  {
    for (size_t iF=0; iF<filenames.size(); iF++)
      filenames[iF].clear();
    filenames.clear();
    floats.clear();
    positions.clear();
    operators.clear();
    std::multimap< std::string, std::string >::const_iterator It;
    for (It=input.begin(); It!=input.end(); It++)
    {
      if((*It).first == "L")
      {
        std::istringstream iss((*It).second);
        iss >> L;
        // std::cout << "L = " << L << std::endl;
        continue;
      }
      if((*It).first == "T")
      {
        std::istringstream iss((*It).second);
        iss >> T;
        // std::cout << "T = " << T << std::endl;
        continue;
      }
      if((*It).first == "position")
      {
        size_t *pos = new size_t[4];
        std::istringstream iss((*It).second);
        iss >> pos[0];
        iss >> pos[1];
        iss >> pos[2];
        iss >> pos[3];
        positions.push_back(pos);
        continue;
      }
      double tmp;
      std::istringstream iss((*It).second);
      iss >> tmp;
      floats.insert(make_pair((*It).first, tmp));

//       if((*It).first == "kappa")
//       {
//         double tmp;
//         std::istringstream iss((*It).second);
//         iss >> tmp;
//         floats.insert(make_pair((*It).first, tmp));
//         continue;
//       }
//       if((*It).first == "mu")
//       {
//         double tmp;
//         std::istringstream iss((*It).second);
//         iss >> tmp;
//         floats.insert(make_pair((*It).first, tmp));
//         continue;
//       }
    }
    for (size_t iF=0; iF<files.size(); iF++)
    {
      std::vector< std::string > tmp;
      //NOTE: some checks are due here!!!
      for (size_t idx=(files[iF]).d_firstIndex; idx<=(files[iF]).d_lastIndex; idx++)
      {
        std::ostringstream oss;
        if((files[iF]).d_directory[((files[iF]).d_directory).size()-1] == '/')
          oss  <<  (files[iF]).d_directory <<  (files[iF]).d_filenameBase;
        else
          oss  <<  (files[iF]).d_directory << "/" <<  (files[iF]).d_filenameBase;
        oss.fill('0');
        oss.width((files[iF]).d_indexWidth);
        oss << idx;
        oss <<  (files[iF]).d_filenameEnding;
        oss.flush();
        tmp.push_back(oss.str());
      }
      filenames.push_back(tmp);
    }
  }
}

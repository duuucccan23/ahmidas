#include "../IO.h"

void Tool::IO::load(Core::Field< QCD::Gauge > *field, std::string const &filename, Tool::IO::filetype type)
{
   switch(type) {
   case Tool::IO::fileILDG :
      Tool::IO::loadILDG(field, filename);
      break;
   default :
      break;
   }

}

void Tool::IO::load(Core::Field< QCD::Spinor > *field, std::string const &filename, Tool::IO::filetype type)
{
   switch(type) {
   case Tool::IO::fileSCIDAC :
      Tool::IO::loadScidac(field, filename);
      break;
   default :
      break; 
   }
}


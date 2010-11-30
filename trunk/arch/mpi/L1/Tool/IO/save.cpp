#include "../IO.h"

namespace Tool
{
  namespace IO
  {
    void saveScidac(Core::Propagator const *propagator, std::vector< std::string > const &filenames);
    void saveScidac(Core::StochasticPropagator< 4 > const *spropagator, std::vector< std::string > const &filenames);
    void saveScidac(Core::StochasticPropagator< 1 > const *spropagator, std::vector< std::string > const &filenames);
  }
}


void Tool::IO::save(Core::Field< QCD::Gauge > const *field, std::string const &filename, filetype type)
{
  switch(type) {
  case fileILDG :
      // at current status, this will give an error message and
      // exit, since not implemented yet
      saveILDG(*field, filename);
      break;
  default :
      break;
  }
}

void Tool::IO::save(Core::Field< QCD::Spinor > const *field, std::string const &filename, filetype type)
{
  switch(type) {
  case fileSCIDAC :
      saveScidac(*field, filename);
      break;
  default :
      break;
  }
}

void Tool::IO::save(Core::Propagator const *propagator, std::vector< std::string> const &filenames, filetype type)
{
  switch(type) {
  case fileSCIDAC :
      saveScidac(propagator, filenames);
      break;
  default :
      break;
  }
}

void Tool::IO::save(Core::StochasticPropagator< 4 > const *sPropagator, std::vector< std::string> const &filenames, filetype type)
{
  switch(type) {
    case fileSCIDAC :
      saveScidac(sPropagator, filenames);
      break;
    default :
      break;
  }
}

void Tool::IO::save(Core::StochasticPropagator< 1 > const *sPropagator, std::vector< std::string> const &filenames, filetype type)
{
  switch(type) {
    case fileSCIDAC :
      saveScidac(sPropagator, filenames);
      break;
    default :
      break;
  }
}


// a lot of data is copied back and forth - this has to be reviewed if more efficiency is desired
void Tool::IO::saveScidac(Core::Propagator const *propagator, std::vector< std::string> const &filenames)
{
  if (filenames.size() < 12)
  {
    std::cerr << "Error in void Tool::IO::loadScidac(Core::Propagator *, std::vector< std::string> const &):"
              << std::endl;
    std::cerr << "filenames.size() should be 12" << std::endl;
    exit(1);
  }

  savePropagator(*propagator, filenames); // Using Lemon
}

void Tool::IO::saveScidac(Core::StochasticPropagator< 4 > const *spropagator, std::vector< std::string> const &filenames)
{
  if (filenames.size() < 4)
  {
    std::cerr << "Error in void Tool::IO::loadScidac(Core::StochasticPropagator<4> *spropagator, std::vector< std::string> const &filenames):"
              << std::endl;
    std::cerr << "filenames.size() should be 4" << std::endl;
    exit(1);
  }

  savePropagator(*spropagator, filenames); // Using Lemon
}

void Tool::IO::saveScidac(Core::StochasticPropagator< 1 > const *spropagator, std::vector< std::string> const &filenames)
{
  if (filenames.size() < 1)
  {
    std::cerr << "Error in void Tool::IO::loadScidac(Core::StochasticPropagator<1> *spropagator, std::vector< std::string> const &filenames):"
              << std::endl;
    std::cerr << "filenames.size() should be 1" << std::endl;
    exit(1);
  }

  savePropagator(*spropagator, filenames); // Using Lemon
}

#include "../IO.h"

namespace Tool
{
  namespace IO
  {
    void saveScidac(Core::Propagator *propagator, std::vector< std::string > const &filenames);
  }
}

void Tool::IO::save(Core::Field< QCD::Gauge > *field, std::string const &filename, filetype type)
{
   switch(type) {
   case Tool::IO::fileSCIDAC :
      Tool::IO::saveScidac(*field, filename);
      break;
   case Tool::IO::fileILDG :
      Tool::IO::saveILDG(*field, filename);
      break;
   default :
      break;
   }
}

void Tool::IO::save(Core::Field< QCD::Spinor > *field, std::string const &filename, filetype type)
{
   switch(type) {
   case Tool::IO::fileSCIDAC :
      Tool::IO::saveScidac(*field, filename);
      break;
   default :
      break;
   }
}

void Tool::IO::save(Core::Propagator *propagator, std::vector< std::string> const &filenames, filetype type)
{
  switch(type) {
  case fileSCIDAC :
      saveScidac(propagator, filenames);
      break;
  default :
      break;
  }
}


// a lot of data is copied back and forth - this has to be reviewed if more efficiency is desired
void Tool::IO::saveScidac(Core::Propagator *propagator, std::vector< std::string> const &filenames)
{
  // would like to do this, but isolate is private
  // propagator->isolate();

  if (filenames.size() == 12)
  {
    Core::Field< QCD::Spinor > tmp [12] =
    {
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T()),
      Core::Field< QCD::Spinor > (propagator->L(), propagator->T())
    };

    Core::Propagator::iterator itTensor = propagator->begin();
    Core::Field< QCD::Spinor >::iterator itsSpinor [12] =
    {
      tmp[ 0].begin(),tmp[ 1].begin(),tmp[ 2].begin(),
      tmp[ 3].begin(),tmp[ 4].begin(),tmp[ 5].begin(),
      tmp[ 6].begin(),tmp[ 7].begin(),tmp[ 8].begin(),
      tmp[ 9].begin(),tmp[10].begin(),tmp[11].begin()
    };

    QCD::Spinor **spinors = new QCD::Spinor *[12];
    for (size_t i=0; i<12; i++)
      spinors[i] = NULL;

    // would like to see for loop here
    // but postfix Field::iterator operator++(int) is not implemented yet
    while (itTensor != propagator->end())
    {
      for (size_t i=0; i<12; i++)
      {
        (*(itsSpinor[i])) = (*itTensor)[i];
        ++(itsSpinor[i]);
      }
      ++itTensor;
    }

    for (size_t i=0; i<12; i++)
    {
      Tool::IO::save(tmp+i, filenames[i], Tool::IO::fileSCIDAC);
    }

  }
  else
  {
    std::cerr << "Error in void Tool::IO::saveScidac(Core::Propagator *, std::vector< std::string> const &):"
              << std::endl;
    std::cerr << "filenames.size() should be 12" << std::endl;
    exit(1);
  }
}

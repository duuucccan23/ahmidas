#include "Tool.ih"

#include <fstream>

// print correlation functions in a quasi-standard ETMC format
void Tool::printLightMesonCorrelator(std::vector< Core::Correlator > const &correlator, std::string const filename)
{

  if (correlator.size() == 0)
  {
    std::cerr << "No data to print. Aborting..." << std::endl;
    exit(1);
  }

  // for compatibility with parallel architecture: serial output
  if (!correlator[0].isRoot())
    return;


  std::ofstream fout(filename.c_str());
  if (!fout.is_open())
  {
    std::cerr << "Could not open file " << filename << ". Aborting..." << std::endl;
    exit(1);
  }

  size_t T = correlator[0].getT();
  size_t L = correlator[0].getL();
  double convention_factor = 1.0/double(L*L*L);

  for  (size_t idx=0; idx<correlator.size(); idx++)
  {

    //fout << "\nGamma combination index " << idx << " \n" << std::endl;
    for (size_t t=0; t<T; t++)
    {
      fout.width(3);
      fout << idx+1 << " ";
      // this is for compatibility with the contraction package light routine
      // and means local-local (no smeared or fuzzed fields involved)
      fout << "  1 ";
      fout.width(3);
      fout << t << " ";
      fout << std::scientific << std::setprecision(10) << std::showpos
                << convention_factor*(correlator[idx].getTrSum(t).real()) << "  "
                << convention_factor*(correlator[idx].getTrSum(t).imag()) << std::endl;
    }
  }
  fout<< std::endl;

}

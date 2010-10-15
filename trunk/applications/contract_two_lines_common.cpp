#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L1/Tool/IO.h>
#include <L2/Contract/Meson.h>
#include <L2/Input/FileReader.h>

int main(int argc, char **argv)
{
  // Ahmidas start(&argc, &argv);

  size_t L_tmp(0);
  size_t T_tmp(0);

  //read inputs from parameter file
  std::map< std::string, double > floats;
  std::vector< std::vector< std::string > > files;
  std::vector< size_t * > positions; // not actually used
  std::vector< int > operators;

  {
    Input::FileReader reader("./contract_two_lines_input.xml");
    reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);
  }

  size_t const L = L_tmp;
  size_t const T = T_tmp;

  Base::Weave weave(L, T);

  double const kappa   = floats["kappa"];
  double const mu_1    = floats["mu_1"];
  double const mu_2    = floats["mu_2"];
  double const thetax  = floats["thetax"];
  double const thetay  = floats["thetay"];
  double const thetaz  = floats["thetaz"];
  double const thetat  = floats["thetat"];

  size_t const timeslice_source  = (size_t) floats["timeslice"];

  bool const phys_base = bool(floats["phys_base"] != 0.0); // 0=no,!=0 =yes

  {
    // write some output concerning what is going to be read
    std::ostringstream out;
    out << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;
    out << "kappa = "  << kappa  << ", mu_1 = "<< mu_1 << ", mu_2 = "<< mu_2 << std::endl;
    out << "thetaX = " << thetax << ", ";
    out << "thetaY = " << thetay << ", ";
    out << "thetaZ = " << thetaz << ", ";
    out << "thetaT = " << thetat << std::endl;
    out << "timeslice of source(s) is " << timeslice_source << std::endl;
    out << "will ";
    if (!phys_base) out << "not ";
    out << "rotate propagators to physical base" << std::endl;
    Print(out.str());
  }


  std::string const &gaugeFieldFile = files[0][0];
  std::vector< std::string > const &propFiles1 = files[1];
  std::vector< std::string > const &propFiles2 = files[2];

  PropagatorType prop1(L, T);
  PropagatorType prop2(L, T);

  Tool::IO::load(&prop1, propFiles1, Tool::IO::fileSCIDAC);
  Print("First propagator successfully loaded\n");

  if (propFiles2 != propFiles1)
  {
    Tool::IO::load(&prop2, propFiles2, Tool::IO::fileSCIDAC);
    Print("Second propagator successfully loaded\n");
  }
  else if (!phys_base)
  {
    Print("Second propagator is the same as first one\n");
    assert(mu_1 == mu_2);
    prop2 = prop1;
  }


  if(phys_base)
  {
    //read the gauge configuration
    Core::Field<QCD::Gauge> gauge_field(L,T);
    Tool::IO::load(&gauge_field, gaugeFieldFile, Tool::IO::fileILDG);
    Print("\ngauge field successfully loaded\n");
    prop1.rotateToPhysicalBasis(mu_1 > 0.0);
    if (propFiles2 == propFiles1)
    {
      Print("Second propagator is the same as first one\n");
      assert(mu_1 == mu_2);
      prop2 = prop1;
    }
    else
    {
      prop2.rotateToPhysicalBasis(mu_2 > 0.0);
    }
  }

  assert(operators.size() % 2 == 0);

  // now we want to tell the contraction routine which operator (gamma) combinations to use
  std::vector< std::pair< Base::Operator, Base::Operator > > operator_combinations;

  {
    for (size_t i = 0; i < operators.size(); i += 2)
    {
      std::ostringstream out;
      out << "operator combination no ";
      out.width(3);
      out << i/2 << ": ";
      out.width(3);
      out << std::showpos << operators[i] << " and ";
      out.width(3);
      out << std::showpos << operators[i+1] << std::endl;
      Print(out.str());
      operator_combinations.push_back(std::make_pair(Tool::convertIntToOperator(operators[i]),
                                                     Tool::convertIntToOperator(operators[i+1])));
    }
  }

#ifdef StocCase
  std::vector< Core::Correlator< Dirac::Matrix > > C2(Contract::light_meson_twopoint_stochastic(prop1, prop2, operator_combinations));
#else
  std::vector< Core::Correlator< Dirac::Matrix > > C2(Contract::light_meson_twopoint(prop1, prop2, operator_combinations));
#endif

  {
    std::ofstream fout("correlators.dat");
    std::ostringstream out;
    out.flush();
    for (size_t iCorr = 0; iCorr < C2.size(); iCorr++)
    {
      C2[iCorr].setOffset(timeslice_source);
      out << C2[iCorr];
    }
    out.flush();
    Print(out.str(), fout);
    fout.close();
  }

  return EXIT_SUCCESS;
}

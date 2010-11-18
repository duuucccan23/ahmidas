
#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Ahmidas.h>
#include <L0/Base/Weave.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Matrix.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>


int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  double const tolerance(1.e-12);

  const size_t L = 4;
  const size_t T = 4;

  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;

  const std::string filename_base1("../../test/source4x4");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    propfilesU.push_back(std::string(filename_base1).append("_u").append(oss.str()));
    propfilesD.push_back(std::string(filename_base1).append("_d").append(oss.str()));
  }

  Core::Propagator uProp(L, T);
  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC);
  std::cout << "u quark propagator successfully loaded\n";

  Core::Propagator dProp(L, T);
  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC);
  std::cout << "d quark propagator successfully loaded\n";


  size_t timeslice_source(0);
  size_t timeslice_sink(T/2);
  size_t timeslice_boundary(T - 1);
  uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);

  // generate sequential sources

  Core::Propagator sequentialSource_fixedProjector_d(L, T);
  Core::Propagator sequentialSource_fixedProjector_u(L, T);
  sequentialSource_fixedProjector_d *= std::complex< double >(0, 0);
  sequentialSource_fixedProjector_u *= std::complex< double >(0, 0);

  Contract::create_sequential_source_proton_d(sequentialSource_fixedProjector_d,
                                              uProp, uProp, timeslice_sink, Base::proj_PARITY_PLUS_TM);

  Contract::create_sequential_source_proton_u(sequentialSource_fixedProjector_u,
                                              dProp, uProp, timeslice_sink, Base::proj_PARITY_PLUS_TM);


  // now we would have to invert the sequential sources and then use the backward propagator for the contraction

  // backward propagators
  Core::Propagator uBwProp(L, T);
  Core::Propagator dBwProp(L, T);

  // have to specify operators
  std::vector< Base::Operator > operators;
  operators.push_back(Base::op_GAMMA_4);

  // the actual contraction
  std::vector< Core::BaryonCorrelator > threepoint =
    Contract::proton_threepoint_sequential(uBwProp, uProp, dBwProp, dProp, NULL, operators);


  if(false)
  {
    std::cout << "\nSUCCESS: threepoint matches reference result!\n" << std::endl;
  }
  else
  {
    std::cerr << "\nFAILURE: threepoint matches reference result!\n" << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}

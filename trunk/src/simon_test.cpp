
#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Meson.h>
#include <L1/Tool/IO.h>



// #define __PRINT__PROPS__

#define __MPI_ARCH__


int main(int argc, char **argv)
{

#ifdef __MPI_ARCH__
  MPI::Init(argc, argv);
  int numprocs(MPI::COMM_WORLD.Get_size());
  int myid(MPI::COMM_WORLD.Get_rank());
#endif

  const size_t L = 4;
  const size_t T = 4;

  std::vector<std::string> propfilesU;

  const std::string filename_base("../test/source.9999.01");
  for (int f=0; f<4; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f%4;
    oss << ".inverted";
    oss.flush();
#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << oss.str() << std::endl;
    propfilesU.push_back(oss.str());
  }

  Core::StochasticPropagator< 4 > *uProp = new Core::StochasticPropagator< 4 >(L, T);

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
    std::cout << "u quark propagator successfully loaded\n" << std::endl;


#ifdef __PRINT__PROPS__

  Core::Propagator::iterator my_iterator = uProp->begin();

  int count(0);
  while(++my_iterator != uProp->end())
  {
 #ifdef __MPI_ARCH__
    if (myid == 0)
    {
 #endif
      std::cout << "Element no. " << count++ << std::endl;
      std::cout << std::endl;
      std::cout << *(my_iterator);
 #ifdef __MPI_ARCH__
    }
 #endif
    ++my_iterator
  }
#endif


  //double normFactor = 1.0/double(L*L*L);
  std::vector< Core::Correlator > C2 = Contract::light_meson_twopoint_stochastic(*uProp, *uProp);

  delete uProp;


  propfilesU.clear();
  const std::string filename_base1("../test/source4x4");
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
#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << propfilesU[f] << std::endl;
  }

  Core::Propagator *uProp12 = new Core::Propagator(L, T);

  Tool::IO::load(uProp12, propfilesU, Tool::IO::fileSCIDAC);

#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << "u quark propagator successfully loaded\n" << std::endl;

  //charged_pion
  Dirac::Gamma5 gamma5;
  Core::Correlator C2_naive = Contract::light_meson_twopoint(uProp12, 0, gamma5, gamma5);


  delete uProp12;

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

#ifdef __MPI_ARCH__
  if (!MPI::Is_finalized())
    MPI::Finalize();
#endif

  return EXIT_SUCCESS;
}

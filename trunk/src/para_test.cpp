

#include <mpi.h>

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


#define __LOAD_PROPAGATORS_ONLY__

int main(int argc, char **argv)
{

  MPI::Init(argc, argv);
  int numprocs(MPI::COMM_WORLD.Get_size());
  int myid(MPI::COMM_WORLD.Get_rank());

  if (myid==0)
    std::cout << "\nprogramm is running on " << numprocs << "cpu(s)\n" << std::endl;

  const size_t L = 16;
  const size_t T = 32;

  std::vector<std::string> propfilesU;

  //const std::string filename_base("../test/source.9999.01");
  const std::string filename_base("/usr1/scratch/dinter/ahmidas_test/source");
  for (int f=0; f<4; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    if (myid==0)
      std::cout << oss.str() << std::endl;
    propfilesU.push_back(oss.str());
  }


  Core::StochasticPropagator< 4 > *uProp = new Core::StochasticPropagator< 4 >(L, T);

  if (myid==0)
    std::cout << "\nmemory for Propagator structure allocated\n" << std::endl;

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  if (myid==0)
    std::cout << "propagator successfully loaded\n" << std::endl;

#ifdef __LOAD_PROPAGATORS_ONLY__
  if (myid==0)
    std::cout << "Done! Finalizing now.\n" << std::endl;
  MPI::Finalize();
  return EXIT_SUCCESS;
#endif


  std::vector< Core::Correlator > C2 = Contract::light_meson_twopoint_stochastic(*uProp, *uProp);

  delete uProp;

  if (myid==0)
    std::cout << "contractions done\n" << std::endl;

  Tool::printLightMesonCorrelator(C2, "./lightMesonCorrelator.txt");

  if (myid==0)
    std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  if (!MPI::Is_finalized())
    MPI::Finalize();

  return EXIT_SUCCESS;
}

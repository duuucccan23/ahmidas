

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


int main(int argc, char **argv)
{

  MPI::Init(argc, argv);
  int numprocs(MPI::COMM_WORLD.Get_size());
  int myid(MPI::COMM_WORLD.Get_rank());

  if (myid==0) 
    std::cout << "\nprogramm is running on " << numprocs << "cpu(s)\n" << std::endl;

  const size_t L = 4;
  const size_t T = 4;

  std::vector<std::string> propfilesU;

  const std::string filename_base("../test/source.9999.01");
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss << filename_base << ".";
    oss.fill('0');
    oss.width(2);
    oss << f%4;
    oss << ".inverted";
    oss.flush();
    if (myid==0)
      std::cout << oss.str() << std::endl;
    propfilesU.push_back(oss.str());
  }


  Core::Propagator *uProp = new Core::Propagator(L, T);

  if (myid==0)
    std::cout << "\nmemory for Propagator structure allocated\n" << std::endl;

  if (uProp->load(propfilesU, "Scidac"))
  {
    if (myid==0)
      std::cout << "u quark propagator successfully loaded\n" << std::endl;
  }
  else
  {
    if (myid==0)
    {
      std::cout << "error reading u quark  propagator\n" << std::endl;
      exit(EXIT_FAILURE);
    }
  }


  std::vector< Core::Correlator > C2 = Contract::light_meson_twopoint_stochastic(*uProp, *uProp);

  delete uProp;

  Tool::printLightMesonCorrelator(C2, "./lightMesonCorrelator.txt");

  if (myid==0)
    std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  if (!MPI::Is_finalized())
    MPI::Finalize();

  return EXIT_SUCCESS;
}
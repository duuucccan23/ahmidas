

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
  int numprocs;
  int myid;
  MPI_Status stat; 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

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
    std::cout << oss.str() << std::endl;
    propfilesU.push_back(oss.str());
  }


  Core::Propagator *uProp = new Core::Propagator(L, T);


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


  double normFactor(1.0/double(L*L*L));
  std::vector< Core::Correlator > C2 = Contract::light_meson_twopoint_stochastic(*uProp, *uProp, normFactor);

  delete uProp;

  if (myid==0) 
    std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  MPI_Finalize();

  return EXIT_SUCCESS;
}

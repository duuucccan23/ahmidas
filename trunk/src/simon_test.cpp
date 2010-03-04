
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
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>


// #define __MPI_ARCH__


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

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
  std::cout << "The following files are going to be read:" << std::endl;

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

  Core::Propagator *uProp = new Core::Propagator(L, T);

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  Core::Propagator *dProp = new Core::Propagator(*uProp);

#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Correlator C2_P = Contract::proton_twopoint(*uProp, *dProp, Base::proj_PARITY_PLUS_TM);

  for (size_t t=0; t<C2_P.getT(); t++)
  {
    std::cout << t << " " << tr(C2_P[t]) << std::endl;
    std::cout << C2_P[t] << std::endl;
  }

  std::cout << "that is supposed to be the result:\n"
    << "0   1.60903097e-03  -7.64473273e-05\n"
    << "1  -9.80937450e-04  -9.66418283e-04\n"
    << "2  -4.03095025e-08   4.29027982e-05\n"
    << "3   9.64650917e-04  -9.44681511e-04"
    << std::endl;

  delete uProp;
  delete dProp;

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


#include <string>
#include <vector>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Core/Field.h>
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
  std::vector<std::string> propfilesD;

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
  std::cout << "The following files are going to be read:" << std::endl;

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
#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << propfilesU[f] << std::endl;
      std::cout << propfilesD[f] << std::endl;
  }

  Core::Propagator *uProp = new Core::Propagator(L, T);

  Tool::IO::load(uProp, propfilesU, Tool::IO::fileSCIDAC);

  #ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << "u quark propagator successfully loaded\n" << std::endl;

  Core::Propagator *dProp = new Core::Propagator(L, T);

  Tool::IO::load(dProp, propfilesD, Tool::IO::fileSCIDAC);

#ifdef __MPI_ARCH__
    if (myid == 0)
#endif
      std::cout << "d quark propagator successfully loaded\n" << std::endl;


  size_t timeslice_source(0);
  size_t timeslice_boundary(T-1);
  uProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  dProp->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);

  Core::Correlator C2_P = Contract::proton_twopoint(*uProp, *dProp, Base::proj_PARITY_PLUS_TM);


  std::cout << "\nreliable code gives this result:\n"
    << "0   1.16145514e-03  -7.64689531e-05\n"
    << "1   1.37746719e-03  -1.02786839e-05\n"
    << "2   4.29020380e-05   3.89420319e-08\n"
    << "3   1.34970449e-03   1.41088319e-05"
    << std::endl;

  delete uProp;
  delete dProp;

  double tolerance = 1.e-9;

  if  (C2_P[0].trace().imag() < (-7.64689531e-05 + tolerance) && C2_P[0].trace().imag() > (-7.64689531e-05 - tolerance)
    && C2_P[1].trace().imag() < (-1.02786839e-05 + tolerance) && C2_P[1].trace().imag() > (-1.02786839e-05 - tolerance)
    && C2_P[2].trace().imag() < ( 3.89420319e-08 + tolerance) && C2_P[2].trace().imag() > ( 3.89420319e-08 - tolerance)
    && C2_P[3].trace().imag() < ( 1.41088319e-05 + tolerance) && C2_P[3].trace().imag() > ( 1.41088319e-05 - tolerance)
    && C2_P[0].trace().real() < ( 1.16145514e-03 + tolerance) && C2_P[0].trace().real() > ( 1.16145514e-03 - tolerance)
    && C2_P[1].trace().real() < ( 1.37746719e-03 + tolerance) && C2_P[1].trace().real() > ( 1.37746719e-03 - tolerance)
    && C2_P[2].trace().real() < ( 4.29020380e-05 + tolerance) && C2_P[2].trace().real() > ( 4.29020380e-05 - tolerance)
    && C2_P[3].trace().real() < ( 1.34970449e-03 + tolerance) && C2_P[3].trace().real() > ( 1.34970449e-03 - tolerance))
  {
  #ifdef __MPI_ARCH__
    if (!MPI::Is_finalized())
      MPI::Finalize();
  #endif
    return EXIT_SUCCESS;
  }



#ifdef __MPI_ARCH__
  if (!MPI::Is_finalized())
    MPI::Finalize();
#endif

  return EXIT_FAILURE;
}

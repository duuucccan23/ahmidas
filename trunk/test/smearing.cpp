#include <cstring>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L1/Tool/IO.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/Weave.h>

int main(int argc, char **argv)
{

  // this test directly compares the result of a smearing procedure to some reference result
  // one could compare the average plaquette as well

  size_t const L = 4;
  size_t const T = 8;

  Base::Weave weave(L, T);

  double const precision_APE    = 1.e-14;
  double const precision_Jacobi = 1.e-14;

  double const APE_alpha = 0.4;
  size_t const APE_iterations = 7;

  double const Jac_alpha = 0.5;
  size_t const Jac_iterations = 5;

  Core::Field< QCD::Gauge > my_gauge_field(L, T);
  Tool::IO::load(&my_gauge_field, "../../test/conf.48", Tool::IO::fileILDG);
  Core::Field< QCD::Spinor > my_spinor_field(L, T);
  Tool::IO::load(&my_spinor_field, "../../test/point_src.48", Tool::IO::fileSCIDAC);

  weave.barrier();

  Smear::APE my_APE_tool(APE_alpha);

//   my_APE_tool.smear(my_gauge_field, APE_iterations, 5);
//   my_APE_tool.smear(my_gauge_field, APE_iterations, 0);

  my_APE_tool.smear(my_gauge_field, APE_iterations);

  weave.barrier();

  Smear::Jacobi my_Jacobi_tool(Jac_alpha);
//   my_Jacobi_tool.smear(&my_spinor_field, my_gauge_field, Jac_iterations, 0);

  my_Jacobi_tool.smear(&my_spinor_field, my_gauge_field, Jac_iterations);

  // test APE smearing
  /* ----------------------------------------------------------------------------- */

  bool testPassed(false);

  size_t const idx_T = 5;
  size_t const idx_X = 2;
  size_t const idx_Y = 1;
  size_t const idx_Z = 3;

  size_t site[4];
  site[Base::idx_X] = idx_X;
  site[Base::idx_Y] = idx_Y;
  site[Base::idx_Z] = idx_Z;
  site[Base::idx_T] = idx_T;

  size_t localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);

  size_t data_rank = weave.rank(site);
//   std::cout << "my rank = "<< weave.rank() << ", data rank = " << data_rank  << ", localIndex = " << localIndex << std::endl;

//   size_t const *coords = weave.d_grid.coords();

//   std::cout << weave.rank() << ": " <<coords[0] <<  coords[1] << coords[2] <<  coords[3] << std::endl;
//   std::cout << weave.rank() <<": dims: " << weave.d_grid.dim(0) <<  weave.d_grid.dim(1) << weave.d_grid.dim(2) <<  weave.d_grid.dim(3) << std::endl;

  /* if data locally available */
  if (data_rank == weave.rank())
  {

    QCD::Gauge const testLinks = my_gauge_field[localIndex];

    // QCD::Gauge consists of 4 SU3::Matrix-es that in turn consist of of 9 std::complex< double >
    std::complex< double > const referenceValues[36] = {
     // mu = Base::idx_X
      std::complex< double >(+3.4668310958013260e-01, -1.9537274568075569e-01),
      std::complex< double >(-2.3880598524204788e-01, -5.0755019874785234e-01),
      std::complex< double >(-7.2216761643778793e-01, +7.4018529498815058e-02),
      std::complex< double >(-5.1255498427208335e-01, -1.1318276494815946e-01),
      std::complex< double >(+6.7380058064865644e-01, -2.7250810532068637e-01),
      std::complex< double >(-2.0656395062741395e-01, +3.9184243538671593e-01),
      std::complex< double >(+7.0724169299244544e-01, +2.5676165852152877e-01),
      std::complex< double >(+3.6943967912988218e-01, +1.4356507479210701e-01),
      std::complex< double >(+9.9929974390575743e-02, +5.1652689320992629e-01),
      // mu = Base::idx_Y
      std::complex< double >(+4.8152211313262633e-01, -6.2530493580424040e-01),
      std::complex< double >(+2.4439045567515388e-02, +4.9612385598841346e-01),
      std::complex< double >(+2.6403695281549389e-01, +2.4633012796243686e-01),
      std::complex< double >(+6.3474816743433499e-02, +5.1513433485282190e-01),
      std::complex< double >(+2.5421825129519737e-01, +8.5134629741020651e-02),
      std::complex< double >(+8.0150234289120958e-01, +1.2777611107832770e-01),
      std::complex< double >(-8.7469747883267987e-02, -3.1636497866153068e-01),
      std::complex< double >(-5.7130253648657214e-01, -5.9582081190758762e-01),
      std::complex< double >(+4.2837745332166677e-01, +1.6542663893448939e-01),
      // mu = Base::idx_Z
      std::complex< double >(+1.1273224981999440e-02, -3.3874845575700502e-02),
      std::complex< double >(-8.6541876004954776e-01, +1.8024993255237626e-01),
      std::complex< double >(+1.1865168963739542e-01, -4.5078544491929390e-01),
      std::complex< double >(-1.1911462424999303e-01, +3.5715355606328147e-01),
      std::complex< double >(+2.7128524522814523e-01, -3.2412892869008259e-01),
      std::complex< double >(-3.1303561480323450e-01, -7.6263130073032070e-01),
      std::complex< double >(+5.2789636636365289e-01, -7.6046293749065441e-01),
      std::complex< double >(+1.7691571818651852e-01, -9.2768063122174199e-02),
      std::complex< double >(-4.4293299640043901e-02, -3.1804805582371165e-01),
      // mu = Base::idx_T
      std::complex< double >(-2.0691771313139423e-01, -4.8421056686791669e-01),
      std::complex< double >(+2.9684056168679718e-01, -5.4882296015843846e-01),
      std::complex< double >(-3.7216664852272946e-01, +4.4147051089928990e-01),
      std::complex< double >(-2.4218049858647311e-01, +1.6102805095450634e-01),
      std::complex< double >(+6.2413621653960649e-01, +4.3941933130055100e-01),
      std::complex< double >(+3.3339666447715788e-01, +4.7077581873829305e-01),
      std::complex< double >(+7.8361159217834719e-01, +1.5523090040273169e-01),
      std::complex< double >(-1.3992658210723971e-01, -9.2001225226543623e-02),
      std::complex< double >(+6.7139748529505897e-02, +5.7385086957141673e-01)};

    QCD::Gauge const referenceLinks(referenceValues);

    if (testLinks.equals(referenceLinks, precision_APE))
    {
      testPassed = true;
      std::cout << "SUCCESS: Gauge links equal reference values at particular lattice site!" << std::endl;
      std::cout << "SUCCESS: APE smearing works!" << std::endl;
    }
    else
    {
      std::cerr << "FAILURE: Gauge links do not equal reference values at particular lattice site!" << std::endl;
      std::cerr << "reference links:\n"        << referenceLinks << std::endl;
      std::cerr << "gauge links calculated:\n" << testLinks      << std::endl;
    }
  }

  weave.barrier();
  weave.broadcast(&testPassed, 1, data_rank);

  if(!testPassed)
    return EXIT_FAILURE;


  // test Jacobi smearing
  /* ----------------------------------------------------------------------------- */

  testPassed = false;

  site[Base::idx_X] = 1;
  site[Base::idx_Y] = 0;
  site[Base::idx_Z] = 2;
  site[Base::idx_T] = 0;

  localIndex = weave.globalCoordToLocalIndex(site[Base::idx_X], site[Base::idx_Y], site[Base::idx_Z], site[Base::idx_T]);

  data_rank = weave.rank(site);

  /* if data locally available */
  if (data_rank == weave.rank())
  {

    QCD::Spinor const testSpinor = my_spinor_field[localIndex];

    // QCD::Gauge consists of 4 SU3::Matrix-es that in turn consist of of 9 std::complex< double >
    std::complex< double > const referenceValues[12] = {
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+0.0, +0.0),
      std::complex< double >(+6.1161734806701711e-05, +1.6548301203740594e-03),
      std::complex< double >(-1.0343940850567506e-03, +1.5980327730610931e-04),
      std::complex< double >(+5.0470395898905528e-04, +1.4899974883655540e-04)};

    QCD::Spinor const referenceSpinor(referenceValues);

    if (testSpinor.equals(referenceSpinor, precision_Jacobi))
    {
      testPassed = true;
      std::cout << "SUCCESS: Spinor equals reference Spinor at particular lattice site!" << std::endl;
      std::cout << "SUCCESS: Jacobi smearing works!" << std::endl;
    }
    else
    {
      std::cerr << "FAILURE: Spinor does not equal reference Spinor at particular lattice site!" << std::endl;
      std::cerr << "reference Spinor:\n"  << referenceSpinor << std::endl;
      std::cerr << "Spinor calculated:\n" << testSpinor      << std::endl;
    }
  }

  weave.barrier();
  weave.broadcast(&testPassed, 1, data_rank);

  if(testPassed)
   return EXIT_SUCCESS;
  else
  return EXIT_FAILURE;
}

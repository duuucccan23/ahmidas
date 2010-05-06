
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

  double const precision_APE    = 1.e-5;
  double const precision_Jacobi = 1.e-5;

  double const APE_alpha = 0.4;
  size_t const APE_iterations = 7;

  double const Jac_alpha = 0.5;
  size_t const Jac_iterations = 5;


  Core::Field< QCD::Gauge > my_gauge_field(L, T);
  Tool::IO::load(&my_gauge_field, "../../test/conf.48", Tool::IO::fileILDG);
  Core::Field< QCD::Spinor > my_spinor_field(L, T);
  Tool::IO::load(&my_spinor_field, "../../test/point_src.48", Tool::IO::fileSCIDAC, 64);

  Smear::APE my_APE_tool(APE_alpha);
  my_APE_tool.smear(my_gauge_field, APE_iterations);


  Smear::Jacobi my_Jacobi_tool(Jac_alpha);
  my_Jacobi_tool.smear(&my_spinor_field, my_gauge_field, Jac_iterations);

  // test APE smearing
  /* ----------------------------------------------------------------------------- */

  bool testPassed(false);

  Base::Weave weave(L, T);

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

  /* if data locally available */
  if (data_rank == weave.rank())
  {

    QCD::Gauge const testLinks = my_gauge_field[localIndex];

    // QCD::Gauge consists of 4 SU3::Matrix-es that in turn consist of of 9 std::complex< double >
    std::complex< double > const referenceValues[36] = {
     // mu = Base::idx_X
      std::complex< double >(+3.4477029755409028e-01, -1.9888641589140504e-01),
      std::complex< double >(-2.4382977498541045e-01, -5.0512987727987058e-01),
      std::complex< double >(-7.2135340394014014e-01, +8.1349554590636572e-02),
      std::complex< double >(-5.1366970220748243e-01, -1.0805505511375929e-01),
      std::complex< double >(+6.7108235741315003e-01, -2.7912274789207958e-01),
      std::complex< double >(-2.0267938681495259e-01, +3.9386237372556521e-01),
      std::complex< double >(+7.0975195341343900e-01, +2.4959435534999930e-01),
      std::complex< double >(+3.7087157542856919e-01, +1.3994314426545071e-01),
      std::complex< double >(+1.0518940994962575e-01, +5.1551935109160107e-01),
      // mu = Base::idx_Y
      std::complex< double >(+4.7977975611333301e-01, -6.2659800402423693e-01),
      std::complex< double >(+2.5765088495681366e-02, +4.9612932608136556e-01),
      std::complex< double >(+2.6468188888744859e-01, +2.4560471634010850e-01),
      std::complex< double >(+6.4838961420946251e-02, +5.1497588721321275e-01),
      std::complex< double >(+2.5456280916019464e-01, +8.4404057504865015e-02),
      std::complex< double >(+8.0175406042085084e-01, +1.2593610498455721e-01),
      std::complex< double >(-8.8471461144364258e-02, -3.1615640466464051e-01),
      std::complex< double >(-5.7289548750810093e-01, -5.9418544552934038e-01),
      std::complex< double >(+4.2884555475076419e-01, +1.6445009953470077e-01),
      // mu = Base::idx_Z
      std::complex< double >(+1.1194757744614696e-02, -3.3994697992655774e-02),
      std::complex< double >(-8.6493202735282881e-01, +1.8242113031241258e-01),
      std::complex< double >(+1.1753030903372362e-01, -4.5113277828232695e-01),
      std::complex< double >(-1.1833840603534554e-01, +3.5745791040506641e-01),
      std::complex< double >(+2.7049028374277129e-01, -3.2482260851296552e-01),
      std::complex< double >(-3.1486378927118375e-01, -7.6184377502952427e-01),
      std::complex< double >(+5.2605949734513935e-01, -7.6170880729204093e-01),
      std::complex< double >(+1.7677610370129626e-01, -9.3223429644285150e-02),
      std::complex< double >(-4.5160974669460931e-02, -3.1793267501023753e-01),
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
      std::complex< double >(+4.79248576360632e-05, +1.65762039468414e-03),
      std::complex< double >(-1.03132081855303e-03, +1.65690419920917e-04),
      std::complex< double >(+5.12548418071813e-04, +1.47231641929049e-04)};

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

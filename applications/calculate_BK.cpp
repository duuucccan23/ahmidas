/*
  This code calculates all three-point and two-point functions to measure B_K
  It uses the one-end-trick, assuming that the propagators are stochastic timeslice propagators
*/

// NOTE:
// the multi-mass solver of the tmLQCD package writes D^{dagger}D to the files.
// the isospin components have to be reconstructed 
// by application of D or D^{dagger}, respectively.


#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>
#include <L0/Core/Correlator.h>
#include <L2/Contract/Meson.h>

int main(int argc, char **argv)
{

  size_t L_tmp = 0;
  size_t T_tmp = 0;

  Input::FileReader reader("./calculate_BK_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);

  size_t const L(L_tmp);
  size_t const T(T_tmp);

  Base::Weave weave(L, T);

  if (weave.isRoot())
    std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu_d    = floats["mu_d"];
  double mu_s    = floats["mu_s"];
  if (weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu_d = " << mu_d << ", mu_s = " << mu_s << std::endl;

  size_t const t_L = size_t(floats["t_L"]);
  assert(t_L >= 0 && t_L < T);
  size_t const t_R = size_t(floats["t_R"]);
  assert(t_R >= 0 && t_R < T);

  if (weave.isRoot())
    std::cout << "timeslice (left)  = " << t_L << std::endl;
  if (weave.isRoot())
    std::cout << "timeslice (right) = " << t_R << std::endl;

  std::vector< std::string > const &propfilesD_L(files[0]);
  std::vector< std::string > const &propfilesS_L(files[1]);
  std::vector< std::string > const &propfilesD_R(files[2]);
  std::vector< std::string > const &propfilesS_R(files[3]);
  std::vector< std::string > const &gaugeFieldFiles(files[4]);



/* ******************************************************************************* */
/* ******** load gauge field ***************************************************** */
/* ******************************************************************************* */


  Core::Field< QCD::Gauge > gauge_field(L, T);
  if (weave.isRoot())
    std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "done.\n" << std::endl;



//TODO: Use de #define directive and one common template file to produce the versions with stochastic, ultrastochastics and normal propagators 

/* ******************************************************************************* */
/* ******** load and prepare propagators from timeslice t_L ********************** */
/* ******************************************************************************* */

  // temporary Propagator just just for loading the MMSolver output
  Core::StochasticPropagator< 4 >  *tmpProp = new Core::StochasticPropagator< 4 > (L, T);

  Tool::IO::load(tmpProp, propfilesD_L, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
    std::cout << "1st propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  dProp_plus_L(L, T);
  Core::StochasticPropagator< 4 >  dProp_minus_L(L, T);
  // reconstruct the full doublet from the combined MMS output, being ~ (DD^dagger)^-1
  tmpProp->reconstruct_doublet(dProp_plus_L, dProp_minus_L, gauge_field, kappa, mu_d);
  if (weave.isRoot())
    std::cout << "1st propagator successfully reconstructed\n" << std::endl;

  Tool::IO::load(tmpProp, propfilesS_L, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
    std::cout << "2nd propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  sProp_plus_L(L, T);
  Core::StochasticPropagator< 4 >  sProp_minus_L(L, T);
  tmpProp->reconstruct_doublet(sProp_plus_L, sProp_minus_L, gauge_field, kappa, mu_s);
  if (weave.isRoot())
    std::cout << "2nd propagator successfully reconsructed\n" << std::endl;

/* ******************************************************************************* */
/* ******** load and prepare propagators from timeslice t_R ********************** */
/* ******************************************************************************* */


  Tool::IO::load(tmpProp, propfilesD_R, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
      std::cout << "3rd propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  dProp_plus_R(L, T);
  Core::StochasticPropagator< 4 >  dProp_minus_R(L, T);
  tmpProp->reconstruct_doublet(dProp_plus_R, dProp_minus_R, gauge_field, kappa, mu_d);
  if (weave.isRoot())
    std::cout << "3rd propagator successfully reconstructed\n" << std::endl;

  Tool::IO::load(tmpProp, propfilesS_R, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
    std::cout << "4th propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  sProp_plus_R(L, T);
  Core::StochasticPropagator< 4 >  sProp_minus_R(L, T);
  tmpProp->reconstruct_doublet(sProp_plus_R, sProp_minus_R, gauge_field, kappa, mu_s);
  if (weave.isRoot())
    std::cout << "4th propagator successfully reconstructed\n" << std::endl;

  delete tmpProp;


// notes: construct objects O as in notes: multiply Dirac structure, dagger, and Tensor multiplication involved

  //gamma structure is: Gamma = gamma_mu (1 - gamma_5)

  Dirac::Gamma< 4 > gamma4;
  Dirac::Gamma< 54 > gamma5gamma4;
  Dirac::Gamma< 1 >  gamma1;
  Dirac::Gamma< 2 >  gamma2;
  Dirac::Gamma< 3 >  gamma3;
  Dirac::Gamma< 45 > gamma0gamma5;
  Dirac::Gamma< 15 > gamma1gamma5;
  Dirac::Gamma< 25 > gamma2gamma5;
  Dirac::Gamma< 35 > gamma3gamma5;
  Dirac::Gamma< 54 > gamma5gamma0;

  //Contractions;

//TODO include the complete gamma structure

//we always need the dagger d propagator
  Core::StochasticPropagator< 4 > dProp_plus_R_dagger(dProp_plus_R);
  dProp_plus_R_dagger.dagger();
  Core::StochasticPropagator< 4 > dProp_minus_L_dagger(dProp_minus_L);
  dProp_minus_L_dagger.dagger();

  // The structure is
 //  phi_s*phi_d ^dagger * gamma5*gammamu(1-gamma5)= phi_s*phi_d ^dagger *(gamma5gammamu+gammamu)


  //RIGHT SIDE
  Core::StochasticPropagator< 4 > *psiR;
  Core::StochasticPropagator< 4 > *O_R;
  Core::Correlator< Dirac::Matrix > *twopointR;
 
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
  (*psiR).rightMultiply(gamma4);
  twopointR = new Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  O_R = new Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma0gamma5);
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma1) ;
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma1gamma5);
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma2);
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma2gamma5);
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma3);
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;
  psiR = new Core::StochasticPropagator< 4 >(dProp_plus_R_dagger);
 (*psiR).rightMultiply(gamma3gamma5);
  *twopointR += Core::Correlator< Dirac::Matrix >(sProp_plus_R*(*psiR));
  *O_R += Core::StochasticPropagator< 4 >(*psiR);
  delete psiR;

  //LEFT SIDE
  Core::StochasticPropagator< 4 > *psiL;
  Core::StochasticPropagator< 4 > *O_L;
  Core::Correlator< Dirac::Matrix > *twopointL;

  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
  (*psiL).rightMultiply(gamma4);
  twopointL = new Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  O_L = new Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply(gamma0gamma5);
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply(gamma1) ;
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply(gamma1gamma5);
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply(gamma2);
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply( gamma2gamma5);
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply(gamma3);
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;
  psiL = new Core::StochasticPropagator< 4 >(dProp_minus_L_dagger);
 (*psiL).rightMultiply(gamma3gamma5);
  *twopointL += Core::Correlator< Dirac::Matrix >(sProp_plus_L*(*psiL));
  *O_L += Core::StochasticPropagator< 4 >(*psiL);
  delete psiL;


  Core::Correlator< Dirac::Matrix >threepoint_Disconnected(*twopointL);
  threepoint_Disconnected *=(*twopointR);


  (*O_R).leftMultiply(sProp_plus_R);
  (*O_L).leftMultiply(sProp_plus_L);
  Core::Correlator< Dirac::Matrix >threepoint_Otto((*O_R) * (*O_L));

  delete O_R;
  delete O_L;

  Core::Correlator< Dirac::Matrix >threepoint(threepoint_Otto);
  threepoint +=threepoint_Disconnected;

  (*twopointL).sumOverSpatialVolume(); // this is the projection to zero momentum
  (*twopointR).sumOverSpatialVolume();
  threepoint_Otto.sumOverSpatialVolume(); 
  threepoint_Disconnected.sumOverSpatialVolume(); 
  threepoint.sumOverSpatialVolume();




  if (weave.isRoot())
  {
    std::cout << "Ahmidas result for twopointL:\n" << std::endl;
           std::cout << (*twopointL) << std::endl;
   std::cout << "Ahmidas result for twopointR:\n" << std::endl;
            std::cout << (*twopointR) << std::endl;
    std::cout << "Ahmidas result for threepoint_Otto:\n" << std::endl;
            std::cout << threepoint_Otto << std::endl;
    std::cout << "Ahmidas result for threepoint_Disconnected:\n" << std::endl;
          std::cout << threepoint_Disconnected << std::endl;
    std::cout << "Ahmidas result for threepoint:\n" << std::endl;
        std::cout << threepoint << std::endl;

  }
  delete twopointL;
  delete twopointR;
  if (weave.isRoot())
    std::cout << "contractions performed successfully\n" << std::endl;
  
  return EXIT_SUCCESS;

  /////////////////////////////
  ////// OLD VERSION //////////
 //////////////////////////////

/*Core::StochasticPropagator< 4 >  O_R = gamma4 * dProp_plus_R;
  O_R += gamma5gamma4 * dProp_plus_R;
  O_R.dagger();
  O_R.rightMultiply(sProp_plus_R); // actually, is it right or left here? check again!


  Core::StochasticPropagator< 4 >  O_L = gamma4 * dProp_minus_L;
  O_L += gamma5gamma4 * dProp_minus_L;
  O_L.dagger();
  O_L.rightMultiply(sProp_plus_L); // actually, is it right or left here? check again!
  
  
  Core::Field < std::complex < double > > correlator = O_R.trace();
  correlator *= O_L.trace();
  correlator *= -1.0; // note that the sign is negative
  
  O_R.leftMultiply(O_L);
  correlator += O_R.trace();*/
  


}

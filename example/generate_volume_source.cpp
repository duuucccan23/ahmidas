/*
   This piece of code generates stochastic volume sources
*/

//ahmidas
#include <L0/Ahmidas.h>

// C++ complex library
#include <complex>

// C++ string library
#include <cstring>

// this we need for storing the input we read from an input file
#include <map>
#include <vector>

// this we need for C++ style output
#include <iomanip>
#include <iostream>


/* *** ahmidas interfaces *** */

// representation of Dirac gamma matrices
#include <L0/Dirac/Gamma.h>

// lattice structure
#include <L0/Core/Field.h>

// propagator structure
#include <L0/Core/Propagator.h>

// IO interface
#include <L1/Tool/IO.h>

// input file reader interface
#include <L2/Input/FileReader.h>

// needed for initialization of random number generator
#include <L0/Base/Random.h>
#include <L0/Base/Z2.h>

// interfaces needed for smearing
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>

#include <L0/Core/Correlator.h>
#include <L2/Contract/Disconnected.h>
#include <L0/Base/Base.h>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  // lattice size, to be read from input file (and therefore initialized to zero)
  size_t L_tmp(0);
  size_t T_tmp(0);


  /* those are the containers to be filled by the input file reader */

  // this one contains floating point variables like kappa, mu, ...
  std::map< std::string, double > floats;
  // this one contains the file names (themselves stored in containers) of propagators (or sources, gauge fields, ...)
  std::vector< std::vector< std::string > > files;


  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */
clock_t start_inp, finish_inp;
start_inp = clock();

  // create input file reader, the name of the input file has to be passed as a parameter
  Input::FileReader reader("./generate_stochastic_source_input.xml");

  // get input parameters
  // note: this is how to invoke a member function of an object in C++:
  // <name of object>.<name of function>(<parameter list>)
  reader.initializeParameters(L_tmp, T_tmp, files, floats);

  const size_t L(L_tmp);
  const size_t T(T_tmp);

  /* ****************************************** */
finish_inp = clock();

  // this is needed if we want to have the output (i.e. to the standard output) done by only
  // one process in the parallel version
  Base::Weave weave(L, T);

if (weave.isRoot())
	        std::cout << "read input in "<< double(finish_inp - start_inp)/CLOCKS_PER_SEC  << "seconds." << std::endl;



  // ##########################################################################################
  // ##########################################################################################


  // that's how writing to the standard output works in C++
  // note: the "<<" operator also works for most of the ahmidas objects like SU3::Spinor or QCD::Tensor

  // that's how one can access the values in the map "floats"
  double kappa = floats["kappa"];
  double mu    = floats["mu"];

  double thetax=floats["thetax"];
  double thetay=floats["thetay"];
  double thetaz=floats["thetaz"];
  double thetat=floats["thetat"];


  //  Prepare for vvSum method
  bool const flag_vvSum = bool(floats["vvSum"] != 0.0); // 0=no,!=0 =yes
  //  need to generate ~ g5 (1 + H)^\dag \xi
  //  maybe better to create a separate executable.

  if(weave.isRoot())
  {
	  std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;
	  std::cout<<"kappa="<<kappa<<", mu="<<mu<<std::endl;
	  std::cout<<"thetax="<<thetax<<", ";
	  std::cout<<"thetay="<<thetay<<", ";
	  std::cout<<"thetaz="<<thetaz<<", ";
	  std::cout<<"thetat="<<thetat<<std::endl;
  }



  uint64_t const rSeed = uint64_t(floats["seed"]);
  if(weave.isRoot())
	  std::cout << "\nrandom seed: " << rSeed << std::endl;


  // this is necessary to have reproducible results for different parallelizations
  Core::Field < uint64_t > seeds(L, T);


  //   if (rSeed > 0)
  //   {
  //     double *tmp = new double[512];
  //     std::generate_n(tmp, 512, Base::Random::Z2);
  //     delete [] tmp;
  //   }
clock_t start_ini, finish_ini;
start_ini = clock();


  // here we assign a seed to each lattice site
  {
	  uint64_t const increment(17);
	  size_t localIndex(0);
	  size_t globalIndex(0);
	  size_t localVolume = weave.localVolume();
	  size_t ctr(0);

	  for(size_t idx_T = 0; idx_T < T; idx_T++)
	  {
		  for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
		  {
			  for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
			  {
				  for(size_t idx_X = 0; idx_X < L; idx_X++)
				  {
					  globalIndex += increment;
					  localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
					  if (localIndex == localVolume)
						  continue;

					  seeds[localIndex] = /*uint64_t(1.e18) **/ (rSeed + uint64_t(globalIndex));
					  // std::cout << "I am lattice site no. " << ctr++ << " and I get a seed of " << seeds[localIndex] << std::endl;
				  }
			  }
		  }
	  }
  }
  // ###############################################################################################
  // let's do it in scalar, even though this is not very elegant!
finish_ini = clock();

if (weave.isRoot())      std::cout <<  "ini seed "<< double(finish_ini - start_ini)/CLOCKS_PER_SEC  << "seconds." << std::endl;



  
  if(weave.isRoot())
  {
	  std::cout << "\nThe following files are going to be created:" << std::endl;

	  // there should only be one container in files, which can be accessed by files[0]
	  // (similar to accessing an object in a C array)
	  for (size_t fileIndex=0; fileIndex<files[0].size(); fileIndex++)
	  {
		  std::cout << (files[0])[fileIndex] << std::endl;
	  }
  }

  Base::SourcePolarization pol_tmp;
  if(int(floats["SourcePolarization"]) == Base::sou_FULLY_POLARIZED)
  {
	  pol_tmp = Base::sou_FULLY_POLARIZED;
  }
  else if(int(floats["SourcePolarization"]) == Base::sou_PARTLY_POLARIZED)
  {
	  pol_tmp = Base::sou_PARTLY_POLARIZED;
  }
  else
  {
	  if(weave.isRoot())
		  std::cerr << "source polarization " << floats["SourcePolarization"] << " unknown or not implemented" << std::endl;
	  exit(1);
  }
  Base::SourcePolarization const polarization(pol_tmp);

  Base::SourceColorState col_tmp;
  if(int(floats["SourceColorState"]) == Base::sou_PURE)
  {
	  col_tmp = Base::sou_PURE;
  }
  else if(int(floats["SourceColorState"]) == Base::sou_GENERIC)
  {
	  col_tmp = Base::sou_GENERIC;
  }
  else
  {
	  if(weave.isRoot())
		  std::cerr << "source color state " << floats["SourceColorState"] << " unknown or not implemented" << std::endl;
	  exit(1);
  }
  Base::SourceColorState const colorState(col_tmp);

  Base::SourceStochasticTypeFlag type_tmp;
  if(int(floats["SourceStochasticTypeFlag"]) == Base::sou_Z4)
  {
	  type_tmp = Base::sou_Z4;
  }
  else if(int(floats["SourceStochasticTypeFlag"]) == Base::sou_Z2)
  {
	  type_tmp = Base::sou_Z2;
  }
  else if(int(floats["SourceStochasticTypeFlag"]) == Base::sou_P1)
  {
	  type_tmp = Base::sou_P1;
  }
  else if(int(floats["SourceStochasticTypeFlag"]) == Base::sou_M1)
  {
	  type_tmp = Base::sou_M1;
  }
  else
  {
	  if(weave.isRoot())
		  std::cerr << "SourceStochasticTypeFlag " << floats["SourceStochasticTypeFlag"]
			  << " unknown or not implemented" << std::endl;
	  exit(1);
  }
  Base::SourceStochasticTypeFlag type(type_tmp);


  double const APE_alpha      = floats["APE_param"];
  size_t const APE_iterations = size_t(floats["APE_steps"]);
  double const Jac_alpha      = floats["Jac_param"];
  size_t const Jac_iterations = size_t(floats["Jac_steps"]);

  if (weave.isRoot())
  {
	  std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
	  std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
  }

  // for simpler use let us reference the file names via convenient names
  std::vector< std::string > const & stochasticSourceFiles(files[0]);
  std::vector< std::string > const & gaugeFieldFiles(files[1]);

  std::vector< std::string > const & stochasticSourceFiles_vvSum(files[2]);

  Core::Field< QCD::Gauge > *gauge_field;

  if (Jac_iterations > 0)
  {
	  gauge_field = new Core::Field< QCD::Gauge > (L, T);
	  if (weave.isRoot())
		  std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
	  Tool::IO::load(gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
	  if (weave.isRoot())
		  std::cout << "done.\n" << std::endl;

	  // smear the gauge field
	  Smear::APE APE_tool(APE_alpha);
	  //APE_tool.smear(*gauge_field, APE_iterations, t_src);
	  APE_tool.smear(*gauge_field, APE_iterations);
	  if (weave.isRoot())
		  std::cout << "gauge field smeared successfully\n" << std::endl;
  }

  // version 1: spin (Dirac) and color dilution
  if(polarization == Base::sou_FULLY_POLARIZED && colorState == Base::sou_PURE)
  {
	  Core::StochasticSource< 12 > stochastic_source(L, T, polarization, colorState, seeds, type);

	  if (Jac_iterations > 0)
	  {
		  //stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field, t_src);
		  stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field);
		  delete gauge_field;
		  if(weave.isRoot())
			  std::cout << "stochastic source smeared successfully\n" << std::endl;
	  }
	  Tool::IO::save(&stochastic_source, stochasticSourceFiles, Tool::IO::fileSCIDAC);
	  if (weave.isRoot())
		  std::cout << "stochastic source saved successfully\n" << std::endl;
  }
  // version 2: spin (Dirac) dilution only
  else if(polarization == Base::sou_FULLY_POLARIZED && colorState == Base::sou_GENERIC)
  {
	  Core::StochasticSource< 4 > stochastic_source(L, T, polarization, colorState, seeds, type);

	  if (Jac_iterations > 0)
	  {
		  //stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field, t_src);
		  stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field);
		  delete gauge_field;
		  if(weave.isRoot())
			  std::cout << "stochastic source smeared successfully\n" << std::endl;
	  }
	  Tool::IO::save(reinterpret_cast< Core::StochasticPropagator< 4 > * >(&stochastic_source),
			  stochasticSourceFiles, Tool::IO::fileSCIDAC);
	  if (weave.isRoot())
		  std::cout << "stochastic source saved successfully\n" << std::endl;
  }
  // version 2: spin (Dirac) dilution only
  else if(polarization == Base::sou_PARTLY_POLARIZED && colorState == Base::sou_GENERIC)
  {

	  clock_t start_source, finish_source;
	  start_source = clock();

	  Core::StochasticSource< 1 > stochastic_source(L, T, polarization, colorState, seeds, type);
	  finish_source = clock();

	  if (weave.isRoot())
		  std::cout <<  "source generated "<< double(finish_source - start_source)/CLOCKS_PER_SEC  << "seconds." << std::endl;



	  if (Jac_iterations > 0)
	  {
		  //stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field, t_src);
		  stochastic_source.smearJacobi(Jac_alpha, Jac_iterations, *gauge_field);
		  delete gauge_field;
		  if(weave.isRoot())
			  std::cout << "stochastic source smeared successfully\n" << std::endl;

	  }

	  clock_t start_save,finish_save;
	  start_save = clock();

	  Tool::IO::save(reinterpret_cast< Core::StochasticPropagator< 1 > * >(&stochastic_source), stochasticSourceFiles, Tool::IO::fileSCIDAC);

	  finish_save = clock();

	  if (weave.isRoot())
		  std::cout << " sources saved in "<< double(finish_save - start_save)/CLOCKS_PER_SEC  << "seconds." << std::endl;



	  if (weave.isRoot())
		  std::cout << "stochastic source saved successfully\n" << std::endl;


	  //bool const flag_vvSum = bool(floats["vvSum"] != 0.0); // 0=no,!=0 =yes
	  //  need to generate ~ g5 (1 + H)^\dag \xi
	  //  maybe better to create a separate executable.
	  if (flag_vvSum == true)
	  {
		  std::vector< std::string > const &gaugeFieldFiles(files[1]); 

		  //read and declare gauge field 
		  Core::Field< QCD::Gauge > gauge_field(L, T);

		  if(weave.isRoot())
		  {
			  std::cout << "\n Now compute and save (1 + H) g5 xi for vvSum" << std::endl;
			  std::cout << "\nThe following files are going to be created:" << std::endl;

			  // there should only be one container in files, which can be accessed by files[0]
			  // (similar to accessing an object in a C array)
			  for (size_t fileIndex=0; fileIndex<files[2].size(); fileIndex++)
			  {
				  std::cout << (files[2])[fileIndex] << std::endl;
			  }
		  }


		  if (weave.isRoot())
			  std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
		  Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
		  if (weave.isRoot())
			  std::cout << "done.\n" << std::endl;

		  Dirac::Gamma<5> gamma5;

		  Core::StochasticPropagator< 1 > xi(stochastic_source);

		  xi.rightMultiply(gamma5);
		  xi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::C);

		  Tool::IO::save(reinterpret_cast< Core::StochasticPropagator< 1 > * > (&xi), stochasticSourceFiles_vvSum, Tool::IO::fileSCIDAC);

		  if (weave.isRoot())
			  std::cout << "Store  (1 + H) g5 xi to test vvSum. \n" << std::endl;

		  double norm_xi(xi.norm());
		  if (weave.isRoot())
			  std::cout << std::scientific << "norm of (1 + H) g5 xi: " << norm_xi << std::endl;

		  // debug  compute xi_tilde^dag xi -> should be real ?

		  std::vector< Base::HermitianBilinearOperator > my_operators;

		  // define a basis of Hermitian Dirac matrices 

		  my_operators.push_back(Base::op_G_0);
		  my_operators.push_back(Base::op_G_1);
		  my_operators.push_back(Base::op_G_2);
		  my_operators.push_back(Base::op_G_3);
		  my_operators.push_back(Base::op_G_4);
		  my_operators.push_back(Base::op_G_5);
		  my_operators.push_back(Base::op_G_6);
		  my_operators.push_back(Base::op_G_7);
		  my_operators.push_back(Base::op_G_8);
		  my_operators.push_back(Base::op_G_9);
		  my_operators.push_back(Base::op_G_10);
		  my_operators.push_back(Base::op_G_11);
		  my_operators.push_back(Base::op_G_12);
		  my_operators.push_back(Base::op_G_13);
		  my_operators.push_back(Base::op_G_14);
		  my_operators.push_back(Base::op_G_15);


		  Core::StochasticPropagator< 1 > phi(L,T);
		  Core::StochasticPropagator< 1 > phi_vv(L,T);

		  Tool::IO::load(&phi,stochasticSourceFiles, Tool::IO::fileSCIDAC);
		  Tool::IO::load(&phi_vv,stochasticSourceFiles_vvSum, Tool::IO::fileSCIDAC);
		  std::vector< Core::Correlator< Dirac::Matrix > > check_source = Contract::compute_loop(phi_vv,phi,my_operators);

		  if(weave.isRoot())
		  {
			  for(size_t i=0; i<my_operators.size(); i++)
			  {
				  for(size_t t = 0; t < T; t++)
				  {

					  std::cout << t << "  "  << i<< "  "<<check_source[i][t].trace().real() <<"  "<< check_source[i][t].trace().imag() << std::endl;

				  }
			  }

		  }
	  }


  }
  else
  {
	  if(weave.isRoot())
		  std::cerr << "source polarization and color state combination not implemented" << std::endl;
	  exit(1);
  }


  // leave main function
  return EXIT_SUCCESS;
}

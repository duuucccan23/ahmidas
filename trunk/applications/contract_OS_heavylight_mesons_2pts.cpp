/*
   This code performs contraction of disconnected for all the element of an hermitian basis of Dirac matrices.
   The code only perform local contraction. Fuzzing has to be include in the future
 */
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


//fstream provides an interface to read and write data from files as input/output streams.
#include <fstream>
#include <sstream>


#include <time.h>

/* *** ahmidas interfaces *** */

// representation of Dirac gamma matrices
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Field.h>
#include <L0/Base/Base.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L2/Contract/Meson.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>


#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Fuzz.h>
#include <L1/Smear/Jacobi.h>

#include <L0/Debug.h>                                                                                                                                    
#include <L0/Print.h>

#include <L0/Ahmidas.h>
#include <L2/Contract/Baryon.h>

#define _with_theta_
#define _with_smearing_
#define _with_momentum_projection


int main(int argc, char **argv)
{

	Ahmidas my_ahmidas(&argc, &argv);

	size_t L = 0;
	size_t T = 0;

	/* ****************************************** */
	/* ****** reading the input file ************ */
	/* ****************************************** */

	Input::FileReader reader("./contract_OS_heavylight_mesons_2pts_input.xml");

	std::map< std::string, double > floats;
	std::vector< size_t * > positions;
	std::vector< int > operators;
	std::vector< std::vector< std::string > > files;

	reader.initializeParameters(L, T, files, floats, positions, operators);



	Base::Weave weave(L, T);


	double kappa = floats["kappa"];
	double mu = floats["mu"];
	double mu_heavy = floats["mu_heavy"];

	double thetax=floats["thetax"];
	double thetay=floats["thetay"];
	double thetaz=floats["thetaz"];
	double thetat=floats["thetat"];


	bool const flag_mms = bool(floats["MMS"] != 0.0); // 0=no,!=0 =yes
	bool debug=(bool)floats["debug"]; // 0=no,!=0 =yes

	if(weave.isRoot() and debug==true) std::cout<<"Debug mode !" <<std::endl;
	if(weave.isRoot() and flag_mms==true) std::cout<<"Input propagators are assumed to be of the form (D^dag D)^-1" <<std::endl;
	if(weave.isRoot() and flag_mms==false) std::cout<<"Input propagators are assumed to be of the form D^-1" <<std::endl;

	if(weave.isRoot())
	{
		std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;
		std::cout<<"kappa="<<kappa<<", mu="<<mu<<", mu_heavy="<<mu_heavy<<std::endl;
		std::cout<<"thetax="<<thetax<<", ";
		std::cout<<"thetay="<<thetay<<", ";
		std::cout<<"thetaz="<<thetaz<<", ";
		std::cout<<"thetat="<<thetat<<std::endl;
	}

	Dirac::Gamma<5> gamma5;
	Dirac::Gamma<1> gamma1;
	Dirac::Gamma<2> gamma2;
	Dirac::Gamma<3> gamma3;
	Dirac::Gamma<4> gamma4;



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

	size_t const timeslice_boundary(T - 1);

	//read the gauge configuration, needed to create u and d propagators


	std::vector< std::string > const &propaLightFiles(files[0]);
	std::vector< std::string > const &gaugeFieldFiles(files[1]); 
	std::vector< std::string > const &propaOSFiles(files[2]);

	//read and declare gauge field 
	Core::Field< QCD::Gauge > gauge_field(L, T);


// because the code is not safe at all  ! exit !
	if (weave.isRoot())
		std::cout << "Not tested  !! exiting...\n" << std::endl;
	exit(1);
// 
	if (weave.isRoot())
		std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
	Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
	if (weave.isRoot())
		std::cout << "done.\n" << std::endl;



	double const APE_alpha      = floats["APE_param"];
	size_t const APE_iterations = size_t(floats["APE_steps"]);
	double const Jac_alpha      = floats["Jac_param"];
	size_t const Jac_iterations = size_t(floats["Jac_steps"]);

	if (weave.isRoot())
	{

#ifdef _with_smearing_
		std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
		std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
#endif
	}



	// read and declare Stochastic prop
	size_t const Nsample=propaLightFiles.size();

	std::vector<std::string> filename;
	for (size_t k=0; k<12;k++) filename.push_back(propaLightFiles[k]);

	if(weave.isRoot())
	{
		std::cout << "\nThe following light propagator is going to be read:" << std::endl;
		for (size_t I=0; I<filename.size(); I++)	std::cout << filename[I] << std::endl;

	}

	Core::Propagator phi(L, T);
	Core::Propagator phi2(L, T);

	if (flag_mms==false) Tool::IO::load(&phi,filename, Tool::IO::fileSCIDAC);


	if (flag_mms==true)
	{
		{
			//	Core::Propagator prop_out(L, T);
			Core::Propagator DD_prop(L, T);
			Tool::IO::load(&DD_prop,filename, Tool::IO::fileSCIDAC,32);

			phi=DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz);
			phi.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
			//				phi = prop_out;
			//  Core::Propagator prop_out(L, T);
			Tool::IO::load(&DD_prop,filename, Tool::IO::fileSCIDAC,32);

			phi2=DD_prop.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz);
			phi2.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
			//              phi = prop_out;

		}
	}
	std::vector<std::string> filename_OS;
	for (size_t k=0; k<12;k++) filename_OS.push_back(propaOSFiles[k]);

	if(weave.isRoot())
	{
		std::cout << "\nThe following light propagator is going to be read:" << std::endl;
		for (size_t I=0; I<filename_OS.size(); I++)	std::cout << filename_OS[I] << std::endl;

	}

	Core::Propagator phi_os(L, T);

	if (flag_mms==false) Tool::IO::load(&phi_os,filename_OS, Tool::IO::fileSCIDAC);


	if (flag_mms==true)
	{
		{
			//	Core::Propagator prop_out(L, T);
			Core::Propagator DD_prop(L, T);
			Tool::IO::load(&DD_prop,filename_OS, Tool::IO::fileSCIDAC,32);

			phi_os=DD_prop.applyDiracOperator(gauge_field,kappa,-mu_heavy,thetat,thetax,thetay,thetaz);
			phi_os.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
			//				phi = prop_out;
		}
	}


	if (weave.isRoot())
		std::cout << "\nReading done.\n" << std::endl;



	if(weave.isRoot()) 
		std::cout<<"Smear gauge field   ... "; 

#ifdef _with_smearing_
	Smear::APE APE_tool(APE_alpha);
	APE_tool.smear(gauge_field, APE_iterations);
#endif


	if(weave.isRoot())	
		std::cout << "done.\n" << std::endl;



#ifdef _with_smearing_
	Smear::Jacobi Jacobi_tool(Jac_alpha);
	phi.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);
	if (flag_mms==true) phi2.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

	if (weave.isRoot())	
		std::cout << "\n smear1 ok\n" << std::endl;

	phi_os.smearJacobi(Jac_alpha, Jac_iterations, gauge_field);

	if (weave.isRoot())	
		std::cout << "\n smear2 ok\n" << std::endl;
#endif

	if (weave.isRoot())
		std::cout << "Compute all bilinear 2pts  " << std::endl ;


	//		std::cout << "  bla bla \n" << std::endl;

	std::vector< std::complex <double>  > C_bilinear_meson1 = Contract::compute_contract_meson(phi,phi_os,my_operators);
	std::vector< std::complex <double>  > C_bilinear_meson2 = Contract::compute_contract_meson(phi_os,phi,my_operators);



	if (flag_mms==true)
	{
		Core::Propagator phi_tmp(phi_os);

		phi_tmp.rotateToPhysicalBasis(true);
		phi2.rotateToPhysicalBasis(false);

		Core::BaryonCorrelator C2_P=Contract::proton_twopoint(phi_tmp, phi_tmp, phi2, Base::proj_NO_PROJECTOR);

		C2_P.deleteField();

		//C2_P.setOffset(0);







		if (weave.isRoot())
		{
			std::ofstream fout("output_2point_unprojected.dat");
			for(size_t t = 0; t < T; t++)
			{
				fout.precision(8);
				fout.width(3);
				fout << t << "  " << std::scientific << std::showpos
					<< C2_P[t][ 0].real() << "  " << C2_P[t][ 0].imag() << "  "
					<< C2_P[t][ 1].real() << "  " << C2_P[t][ 1].imag() << "  "
					<< C2_P[t][ 2].real() << "  " << C2_P[t][ 2].imag() << "  "
					<< C2_P[t][ 3].real() << "  " << C2_P[t][ 3].imag() << "  "
					<< std::endl;
				fout.width(3);
				fout << t << "  " << std::scientific << std::showpos
					<< C2_P[t][ 4].real() << "  " << C2_P[t][ 4].imag() << "  "
					<< C2_P[t][ 5].real() << "  " << C2_P[t][ 5].imag() << "  "
					<< C2_P[t][ 6].real() << "  " << C2_P[t][ 6].imag() << "  "
					<< C2_P[t][ 7].real() << "  " << C2_P[t][ 7].imag() << "  "
					<< std::endl;
				fout.width(3);
				fout << t << "  " << std::scientific << std::showpos
					<< C2_P[t][ 8].real() << "  " << C2_P[t][ 8].imag() << "  "
					<< C2_P[t][ 9].real() << "  " << C2_P[t][ 9].imag() << "  "
					<< C2_P[t][10].real() << "  " << C2_P[t][10].imag() << "  "
					<< C2_P[t][11].real() << "  " << C2_P[t][11].imag() << "  "
					<< std::endl;
				fout.width(3);
				fout << t << "  " << std::scientific << std::showpos
					<< C2_P[t][12].real() << "  " << C2_P[t][12].imag() << "  "
					<< C2_P[t][13].real() << "  " << C2_P[t][13].imag() << "  "
					<< C2_P[t][14].real() << "  " << C2_P[t][14].imag() << "  "
					<< C2_P[t][15].real() << "  " << C2_P[t][15].imag() << "  "
					<< std::endl;

			}
			fout.close();
		}
	}
	if (weave.isRoot())
		std::cout << "contractions performed and saved successfully\n" << std::endl;

	// for "vv" correlators
	// the full Dirac structure is kept for a reason
	Core::Correlator< Dirac::Matrix > *pion_twopoint;
	phi_os.revert();

	phi *= gamma5;
	phi_os *= gamma5;




	pion_twopoint = new Core::Correlator< Dirac::Matrix >(phi*phi_os);


	pion_twopoint->sumOverSpatialVolume();


	if (weave.isRoot())
	{
		std::cout << "\nPion two point function:" << std::endl;
		for (size_t t=0; t<T; t++)
		{
			// this is the way formatted output works in C++
			std::cout.width(3);
			std::cout << t << "  ";
			std::cout.width(20);
			// since the Correlator stores complex numbers, we can access the real and imaginary parts using
			// the complex class member functions real() and imag()
			std::cout << std::fixed << std::setprecision(20) << std::showpos << ((*pion_twopoint)[t]).trace().real() << "  ";
			std::cout.width(20);
			std::cout << std::fixed << std::setprecision(20) << std::showpos << ((*pion_twopoint)[t]).trace().imag() << std::endl;
		}

	}


	if (weave.isRoot())	
	{
		std::ofstream fout_bilinear_meson1;
		std::ofstream fout_bilinear_meson2;


		fout_bilinear_meson1.open("output_os_meson1.dat");
		fout_bilinear_meson2.open("output_os_meson2.dat");

		clock_t start, finish;
		start = clock();


		for(size_t i=0; i<my_operators.size(); i++)
		{

			for(size_t t = 0; t < T; t++)
			{
				fout_bilinear_meson1 << t << std::scientific <<"  "
					<< i <<"  "<<C_bilinear_meson1[t +  T*i].real() <<"  "<< C_bilinear_meson1[t + T*i].imag() << std::endl;
				fout_bilinear_meson2 << t << std::scientific <<"  "
					<< i <<"  "<<C_bilinear_meson2[t +  T*i].real() <<"  "<< C_bilinear_meson2[t + T*i].imag() << std::endl;


				fout_bilinear_meson1.flush();
				fout_bilinear_meson2.flush();

			}
		}			
		fout_bilinear_meson1.close();
		fout_bilinear_meson2.close();

		finish = clock();
		std::cout << "done in "<< double(finish - start)/CLOCKS_PER_SEC  << std::endl;

	}/* end if weave.isRoot() */

#ifdef __MPI_ARCH__
	if (myid == 0)
#endif

		if(weave.isRoot())
			std::cout << "\nprogram is going to exit normally now\n" << std::endl;

	return EXIT_SUCCESS;

}



//check my functions :

/*	Core::StochasticPropagator< 1 >  tmp(L, T);
#ifdef _with_theta_
tmp = phi.applyHOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz);
#else
tmp = phi.applyHOperator(gauge_field,kappa,mu, timeslice_boundary);
#endif
Core::StochasticPropagator< 1 >  tmp2(L, T);
#ifdef _with_theta_
tmp2 = phi.applyBdaggerOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz);
#else
tmp2 = phi.applyBdaggerOperator(gauge_field,kappa,mu,timeslice_boundary);
#endif

tmp2 *= std::complex<double>(1.0 + 4.*kappa*kappa*mu*mu,0.0);

tmp += tmp2;


tmp *= 1./(2.*kappa);

Core::StochasticPropagator< 1 >  diff(L, T);
diff = tmp;
diff -=  source;

double const diffNorm = diff.norm();
if (weave.isRoot())
std::cout << "CHECK : norm of difference between source and source obtained with H and Bdagger operator: " << diffNorm << std::endl;
 */



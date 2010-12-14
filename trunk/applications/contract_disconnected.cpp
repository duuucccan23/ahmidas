
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
#include <L2/Contract/Disconnected.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>


#include <L1/Smear.h>
#include <L1/Smear/APE.h>
//#include <L1/Smear/Fuzz.h>
#include <L1/Smear/Jacobi.h>

#include <L0/Debug.h>                                                                                                                                    
#include <L0/Print.h>

#include <L0/Ahmidas.h>


#define _with_theta_
//#define _with_MMS_

int main(int argc, char **argv)
{

	Ahmidas my_ahmidas(&argc, &argv);

	size_t L = 0;
	size_t T = 0;

	/* ****************************************** */
	/* ****** reading the input file ************ */
	/* ****************************************** */

	Input::FileReader reader("./contract_disconnected_input.xml");

	std::map< std::string, double > floats;
	std::vector< size_t * > positions;
	std::vector< int > operators;
	std::vector< std::vector< std::string > > files;

	reader.initializeParameters(L, T, files, floats, positions, operators);



	Base::Weave weave(L, T);


	double kappa = floats["kappa"];
	double mu = floats["mu"];
	double thetax=floats["thetax"];
	double thetay=floats["thetay"];
	double thetaz=floats["thetaz"];
	double thetat=floats["thetat"];



	if(weave.isRoot())
	{
		std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;
		std::cout<<"kappa="<<kappa<<", mu="<<mu<<std::endl;
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



	std::vector< Base::Operator > my_operators;

// not exactly the convention used by carsten 
//  missing are : g5 g0 g1, g5 g0 g2, g5 g0 g2
	my_operators.push_back(Base::op_GAMMA_5);
	my_operators.push_back(Base::op_GAMMA_1);
	my_operators.push_back(Base::op_GAMMA_2);
	my_operators.push_back(Base::op_GAMMA_3);
	my_operators.push_back(Base::op_GAMMA_45);
	my_operators.push_back(Base::op_GAMMA_14);
	my_operators.push_back(Base::op_GAMMA_24);
	my_operators.push_back(Base::op_GAMMA_34);
	my_operators.push_back(Base::op_UNITY);
	my_operators.push_back(Base::op_GAMMA_15);
	my_operators.push_back(Base::op_GAMMA_25);
	my_operators.push_back(Base::op_GAMMA_35);
	my_operators.push_back(Base::op_GAMMA_4);

	size_t const timeslice_boundary(T - 1);

	//read the gauge configuration, needed to create u and d propagators


	std::vector< std::string > const &propaStochaFiles(files[0]);
	std::vector< std::string > const &gaugeFieldFiles(files[1]); 

	//read source if possible ?
	//std::vector< std::string > const &gaugeFieldFiles(files[1]); 

	//read and declare gauge field 
	Core::Field< QCD::Gauge > gauge_field(L, T);


	if (weave.isRoot())
		std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
	Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
	if (weave.isRoot())
		std::cout << "done.\n" << std::endl;



	Core::Field< QCD::Gauge > gauge_field_f(gauge_field);

	double const APE_alpha      = floats["APE_param"];
	size_t const APE_iterations = size_t(floats["APE_steps"]);
	double const Nlong    = floats["Fuzz_param"];
	double const Jac_alpha      = floats["Jac_param"];
	size_t const Jac_iterations = size_t(floats["Jac_steps"]);

	if (weave.isRoot())
	{
		std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
		std::cout << "Fuzzing parameter = " << Nlong <<  std::endl;
		std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;

	}


	if(weave.isRoot()) 
		std::cout<<"Smear gauge field   ... "; 

	Smear::APE APE_tool(APE_alpha);
	//	Smear::Fuzz Fuzz_tool(Nlong);
	APE_tool.smear(gauge_field_f, APE_iterations);
	//Fuzz_tool.smear(gauge_field_f<L,T>);
	if(weave.isRoot())	
		std::cout << "done.\n" << std::endl;



	// read and declare Stochastic prop
	size_t const Nsample=propaStochaFiles.size();
	std::cout << "Nsample = "<< Nsample << "std::endl";

	Core::StochasticPropagator<1> source(L,T);
	Core::StochasticPropagator< 1 >  phi(L, T);
	Core::StochasticPropagator< 1 >  tmp(L, T);

	for(size_t n=0;n<Nsample;n++)
	{  

		if(weave.isRoot()) std::cout<< "Stochastic propagator " <<n<< " to be read from " << propaStochaFiles[n] <<" ... "; 

		std::vector<std::string> filename;
		filename.push_back(propaStochaFiles[n]);

		Tool::IO::load(&phi,filename, Tool::IO::fileSCIDAC);
		if (weave.isRoot())
			std::cout << "done.\n" << std::endl;


		// if source have been produced using MMS ... need to apply Dirac operator to get 1/D
#ifdef _with_MMS_
		if(weave.isRoot())	
			std::cout << "MMS :" << propaStochaFiles[n] << "is 1/(Ddagger D) so apply Ddagger to get the propagator" <<" ...";
		exit(1);
		if(weave.isRoot())	
			std::cout << "done.\n" << std::endl;
#endif	
		Core::StochasticPropagator< 1 >  phi_f(phi);
		Smear::Jacobi Jacobi_tool(Jac_alpha);
		phi_f.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_f);



		//Now compute the source 


		if(weave.isRoot()) 
			std::cout << "Apply Dirac operator to get the source " << " ... ";

		source = phi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Full); 
		
		if (weave.isRoot())
			std::cout << "done.\n" << std::endl;

		// Now compute g5 [B^dagger H]^4 g5
		if(weave.isRoot()) 
			std::cout<<"Compute g5 [B^dagger H]^4 g5 times the source field "<< " ... ";

		Core::StochasticPropagator< 1 >  xi(source);
		xi.rightMultiply(gamma5);

		for(size_t i=0;i<4;i++)
		{	
			/* apply Hopping part ...*/
			tmp = xi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::H);
			/* apply Bdagger ... */
			xi = tmp.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Bdagger);
		}

		xi.rightMultiply(gamma5);

		//g5 [B^dagger H]^4 g5 computed

		//to have the same normalization than carsten . Origin ?
		xi *= 1./(4.*kappa);

		if (weave.isRoot())
			std::cout << "done.\n" << std::endl;

		// for "vv" correlators

		Core::StochasticPropagator< 1 >  psi(phi);
		Core::StochasticPropagator< 1 >  g5_phi(phi);

		g5_phi.rightMultiply(gamma5);


		if (weave.isRoot())
			std::cout << "Compute loops for vv and v4 method"<<" ... ";
		// v4 loop
		std::vector< Core::Correlator< Dirac::Matrix > > C_v4 = Contract::compute_loop(xi,phi,my_operators);
		//vv loop
		std::vector< Core::Correlator< Dirac::Matrix > > C_vv = Contract::compute_loop(phi,g5_phi,my_operators);


		for(size_t i=0; i<my_operators.size(); i++)
		{
			C_v4[i].sumOverSpatialVolume();
			C_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	
			C_vv[i].sumOverSpatialVolume();
		}


		//output
		if (weave.isRoot())	
		{
			std::ofstream fout_v4;
			std::ofstream fout_vv;

			if (n==0)
			{
				fout_v4.open("output_disc_v4.dat");
				fout_vv.open("output_disc_vv.dat");
			}
			else
			{
				fout_v4.open("output_disc_v4.dat",std::ios::app);
				fout_vv.open("output_disc_vv.dat",std::ios::app);

			}

			for(size_t i=0; i<my_operators.size(); i++)
			{
				for(size_t t = 0; t < T; t++)
				{
					fout_v4 << t << std::scientific <<"  "
						<< i <<"  "<< n <<"  "<<C_v4[i][t].trace().real() <<"  "<< C_v4[i][t].trace().real() << std::endl;
					fout_vv << t << std::scientific <<"  "
						<< i <<"  "<< n <<"  "<<C_vv[i][t].trace().real() <<"  "<< C_vv[i][t].trace().real() << std::endl;

				}			
			}
			fout_v4.close();
			fout_vv.close();


		}

		if (weave.isRoot())
			std::cout << "done.\n" << std::endl;

	}


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



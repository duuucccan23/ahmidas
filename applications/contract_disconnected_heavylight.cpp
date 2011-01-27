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
#include <L2/Contract/Disconnected.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>


#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Fuzz.h>
#include <L1/Smear/Jacobi.h>

#include <L0/Debug.h>                                                                                                                                    
#include <L0/Print.h>

#include <L0/Ahmidas.h>


#define _with_theta_
#define _with_MMS_
//#define _without_MMS_

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
	double mu_s = floats["mu_s"];

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


	std::vector< std::string > const &propaStochaFiles(files[0]);
	std::vector< std::string > const &gaugeFieldFiles(files[1]); 
#ifdef _with_heavy_light_
	std::vector< std::string > const &propaStochaFiles_light(files[2]); 
#endif
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
	Smear::Fuzz Fuzz_tool(Nlong);
	APE_tool.smear(gauge_field_f, APE_iterations);
	Fuzz_tool.smear(gauge_field_f);

	if(weave.isRoot())	
		std::cout << "done.\n" << std::endl;



	// read and declare Stochastic prop
	size_t const Nsample=propaStochaFiles.size();
	if (weave.isRoot())
		std::cout << "Number of stochastic propagator to read : "<< Nsample << std::endl;



	if (Nsample%12 != 0){ 

		if (weave.isRoot())
			std::cout << "reads only multiple of 12 sources for efficiency  ! " << std::endl;

		exit(1);	
	}
	else { //  Nsample==12

		if(weave.isRoot()) std::cout<< "The code will read sources by packet of 12 " <<" ... ";

		size_t N = size_t(Nsample/12);
		for (size_t n=0;n<N;n++)
		{

			Core::Propagator source(L,T);
			Core::Propagator tmp(L, T);

			std::vector<std::string> filename;
			for (size_t k=0; k<12;k++) filename.push_back(propaStochaFiles[k+12*n]);

			if(weave.isRoot())
			{
				std::cout << "\nThe following files are going to be read:" << std::endl;
				for (size_t I=0; I<filename.size(); I++)	std::cout << filename[I] << std::endl;

			}

#ifdef _with_heavy_light_

			std::vector<std::string> filename_light;
			for (size_t k=0; k<12;k++) filename_light.push_back(propaStochaFiles_light[k+12*n]);
			if(weave.isRoot())
			{
				std::cout << "\nThe following light propagator  are going to be read:" << std::endl;
				for (size_t I=0; I<filename_light.size(); I++)	std::cout << filename_light[I] << std::endl;

			}
			Core::Propagator light_DD_prop(L, T);
			Core::Propagator prop_light_out(L, T);
			Core::Propagator prop_light_out_opp(L, T);

			Tool::IO::load(&light_DD_prop,filename_light, Tool::IO::fileSCIDAC,32);

			prop_light_out=light_DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz);
			prop_light_out.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
			Core::Propagator phi_light(prop_light_out);

			prop_light_out_opp=light_DD_prop.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz);
			prop_light_out_opp.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator

			Core::Propagator phi_light_opp(prop_light_out_opp);

			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;
#endif


#ifdef _without_MMS_
			Core::Propagator phi(L, T);
			Tool::IO::load(&phi,filename, Tool::IO::fileSCIDAC);
#endif
			// if source have been produced using MMS ... read (D+D-)^-1
#ifdef _with_MMS_

			if(weave.isRoot()) 
				std::cout << "MMS compilation  1/(Ddagger D) is read so  Ddagger is applied to get the propagator" <<" ...";

			Core::Propagator DD_prop(L, T);
			Core::Propagator prop_out(L, T);

			Tool::IO::load(&DD_prop,filename, Tool::IO::fileSCIDAC,32);

			prop_out=DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz);
			prop_out.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator

			Core::Propagator phi(prop_out);
#ifdef _with_heavy_light_
			//build the propagator with opposite mass
			Core::Propagator prop_opp(L, T);
			prop_opp=DD_prop.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz);
			prop_opp.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
			Core::Propagator phi_opp(prop_opp);
#endif

#endif

			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;


			//Now compute the source 


			if(weave.isRoot()) 
				std::cout << "Apply Dirac operator to get the source " << " ... ";

			source = phi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Full); 

			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;

				// for "vv" correlators

			Core::Propagator  psi(phi);
			Core::Propagator  g5_phi(phi);

			g5_phi.rightMultiply(gamma5);


			if (weave.isRoot())
				std::cout << "Compute loops of the form 1/M_u + 1/M_d  - 2/M_s "<<" ... ";


#ifdef _with_heavy_light_


			Core::Propagator phi_light_conj(phi_light);
			Core::Propagator phi_light_opp_conj(phi_light);
			phi_light_conj.conjugate();
			phi_light_opp_conj.conjugate();
			Core::Propagator  g5_phi_opp(phi_opp);
			g5_phi_opp.rightMultiply(gamma5);

// compute for the bilinear operator first
			//1/M_u  - 1/M_s
			std::vector< std::complex <double>  > C_hl1 = Contract::compute_loop_new(phi_light_conj,g5_phi_opp,my_operators);
			//1/M_d  - 1/tilde{M}_s
			std::vector< std::complex <double>  > C_hl2 = Contract::compute_loop_new(phi_light_opp_conj,g5_phi,my_operators);

#endif

			//loop involving twist 2 operators
#ifdef	_with_Omunu_
			std::vector< std::complex<double> > C_twist2 = Contract::compute_loop_twist2_operator(gauge_field,xi,phi);
			std::vector< std::complex<double> > C_twist2_pol = Contract::compute_loop_twist2_operator(gauge_field,xi,g5_phi);
#endif
		//	for(size_t i=0; i < C_vv.size(); i++)
		//	{
				//C_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
		//	}


			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;


			if (weave.isRoot())	
			{

#ifdef _with_Omunu_
				std::ofstream fout_twist2;
				std::ofstream fout_twist2_pol;
#endif

				if (n==0)
				{

#ifdef _with_Omunu_
					fout_twist2.open("output_disc_twist2.dat");
					fout_twist2_pol.open("output_disc_twist2_pol.dat");
#endif
				}
				else	
				{

#ifdef _with_Omunu_
					fout_twist2.open("output_disc_twist2.dat",std::ios::app);
					fout_twist2_pol.open("output_disc_twist2_pol.dat",std::ios::app);
#endif
				}

				std::cout << "write output for sources from " << 12*n << " to " << 11+12*n<<" ... "  << std::endl;


#ifdef _with_Omunu_

				for(size_t j=0; j<12; j++)
				{
					for(size_t i=0; i<4; i++)
					{
						for(size_t k=0; k<4; k++)
						{
							for(size_t t = 0; t < T; t++)
							{



								fout_twist2 << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<< k <<"  " <<C_twist2[t  + T*k + 4*T*i + 4*T*4*j ].real() <<"  "<< C_twist2[t + T*k + 4*T*i +4*T*4*j ].imag() << std::endl;
								fout_twist2_pol << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<< k <<"  " <<C_twist2_pol[t  + T*k + 4*T*i + 4*T*4*j ].real() <<"  "<< C_twist2_pol[t + T*k + 4*T*i +4*T*4*j ].imag() << std::endl;

								fout_twist2.flush();
								fout_twist2_pol.flush();

							}/*t*/
						}/*end loop on the 4 terms contributing to a symmetrized covariant derivative */
					} /*i loop on the 4 operators : O_mumu (unpo1arized) and Omumu g5 (polarized) */
				} /*loop on the 12 sources*/

				fout_twist2.close();
				fout_twist2_pol.close();
#endif

				std::cout << "done.\n" << std::endl;

			}/* end if weave.isRoot() */
		} /* loop on Nsample */
	}/* end case N%12 = 0 */

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



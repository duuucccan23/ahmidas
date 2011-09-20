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


#define _with_Fuzzing_
#define _with_theta_
#define _with_Omunu_
#define _with_momentum_projection


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


    bool const flag_mms = bool(floats["MMS"] != 0.0); // 0=no,!=0 =yes
    bool debug=(bool)floats["debug"]; // 0=no,!=0 =yes

	 if(weave.isRoot() and debug==true) std::cout<<"Debug mode !" <<std::endl;
	 if(weave.isRoot() and flag_mms==true) std::cout<<"Input propagators are assumed to be of the form (D^dag D)^-1" <<std::endl;
	 if(weave.isRoot() and flag_mms==false) std::cout<<"Input propagators are assumed to be of the form D^-1" <<std::endl;
	
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
	std::vector< std::string > const &propaStochaFiles_vvSum(files[2]); 



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
#ifdef _with_Fuzzing__
		std::cout << "Fuzzing parameter = " << Nlong <<  std::endl;

#endif

#ifdef _with_APE_
		std::cout << "APE    smearing: parameter = " << APE_alpha << ", iterations = " << APE_iterations << std::endl;
		std::cout << "Jacobi smearing: parameter = " << Jac_alpha << ", iterations = " << Jac_iterations << std::endl;
#endif
	}


	if(weave.isRoot()) 
		std::cout<<"Smear gauge field   ... "; 

#ifdef _with_APE_
	Smear::APE APE_tool(APE_alpha);
	APE_tool.smear(gauge_field_f, APE_iterations);
#endif

#ifdef _with_Fuzzing__
	Smear::Fuzz Fuzz_tool(Nlong);
	Fuzz_tool.smear(gauge_field_f);
#endif

	if(weave.isRoot())	
		std::cout << "done.\n" << std::endl;



	// read and declare Stochastic prop
	size_t const Nsample=propaStochaFiles.size();
	if (weave.isRoot())
		std::cout << "Number of stochastic propagator to read : "<< Nsample << std::endl;



	if (Nsample%12 != 0){ 


	}
	else { //  Nsample==12

		if(weave.isRoot()) std::cout<< "The code will read sources by packet of 12 " <<" ... ";

		size_t N = size_t(Nsample/12);
		for (size_t n=0;n<N;n++)
		{


			// reading standard prop

			std::vector<std::string> filename;
			for (size_t k=0; k<12;k++) filename.push_back(propaStochaFiles[k+12*n]);

			if(weave.isRoot())
			{
				std::cout << "\nThe following files are going to be read:" << std::endl;
				for (size_t I=0; I<filename.size(); I++)	std::cout << filename[I] << std::endl;

			}

			Core::Propagator phi(L, T);

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
				}
			}

			std::vector<std::string> filename_vv;
			for (size_t k=0; k<12;k++) filename_vv.push_back(propaStochaFiles_vvSum[k+12*n]);

			if(weave.isRoot())
			{
				std::cout << "\nThe following files are going to be read:" << std::endl;
				for (size_t I=0; I<filename_vv.size(); I++)	std::cout << filename_vv[I] << std::endl;

			}

			Core::Propagator phi_vv(L, T);

			if (flag_mms==false) Tool::IO::load(&phi_vv,filename_vv, Tool::IO::fileSCIDAC);


			if (flag_mms==true)
			{
				{
					//	Core::Propagator prop_out(L, T);
					Core::Propagator DD_prop(L, T);
					Tool::IO::load(&DD_prop,filename_vv, Tool::IO::fileSCIDAC,32);

					phi_vv=DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz);
					phi_vv.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
					//				phi = prop_out;
				}
			}


			if (weave.isRoot())
				std::cout << "\nReading done.\n" << std::endl;

			// for "vv" correlators

			Core::Propagator  g5_phi_vv(phi_vv);
			g5_phi_vv.rightMultiply(gamma5);



		
		if (weave.isRoot())
			std::cout << "Compute loops for vv Sum method"<<" ... ";
			std::vector< std::complex <double>  > C_vv = Contract::compute_loop_new(g5_phi_vv,phi,my_operators);
			std::vector< std::complex <double> > C_conserved_vv = Contract::compute_loop_conserved_vector_current(gauge_field,g5_phi_vv,phi);


		for(size_t i=0; i < C_vv.size(); i++)
		{
			C_vv[i] *=  std::complex<double>(2.0/(16.*kappa*kappa),0);	 // factor = + 2 /( 4 kappa) ^2
		}

		for(size_t i=0; i < C_conserved_vv.size(); i++)
		{
			C_conserved_vv[i] *=  std::complex<double>(2.0/(16.*kappa*kappa),0);	 // factor + 2 /( 4 kappa) ^2
		}

		if (weave.isRoot()) std::cout << "done.\n" << std::endl;


		//sum_disc ? 
		//			if (weave.isRoot())
		//				std::cout << "Compute loops for vv_sum  method"<<" ..."; 

		//			Core::Propagator  xi(L,T);
		//			xi = phi.conjugate();
		//			xi.rightMultiply(gamma5);
		//			xi.leftMultiply(gamma5);
		//			xi = xi.applyDiracOperator(gauge_field,kappa,0,thetat,thetax,thetay,thetaz,Base::Full); //  D(mu=0) = 1+H
		//			xi *= 2;
		//			xi.conjugate();

		//			std::vector< std::complex <double>  > C_vv_sum = Contract::compute_loop_new(xi,phi,my_operators);

		//			if (weave.isRoot()) std::cout << "done.\n" << std::endl;


#ifdef _with_Omunu_
		if (weave.isRoot())
			std::cout << "Twist 2 operators  "<<" ... ";


		std::vector< std::complex<double> > C_twist2_vv = Contract::compute_loop_twist2_operator(gauge_field,g5_phi_vv,phi,false);
		std::vector< std::complex<double> > C_twist2_pol_vv = Contract::compute_loop_twist2_operator(gauge_field,g5_phi_vv,phi,true);
		for(size_t i=0; i < C_twist2_vv.size(); i++)
		{
			C_twist2_vv[i] *=  std::complex<double>(2.0/(16.*kappa*kappa),0);	 // factor = + 2 /( 4 kappa) ^2
			C_twist2_pol_vv[i] *=  std::complex<double>(2.0/(16.*kappa*kappa),0);	 // factor = + 2 /( 4 kappa) ^2
		}


		if (weave.isRoot())
			std::cout << "done.\n" << std::endl;


#endif



#ifdef _with_momentum_projection

		std::vector< std::complex <double>  > C_vv_mom;
		if (weave.isRoot())
			std::cout << "Non zero momentum "<<" ... ";


		size_t const * const source_position = positions[0];
		size_t const timeslice_source = source_position[Base::idx_T] % T;

		if(weave.isRoot())
			std::cout << " The position of the source is (x,y,z,t) = " << source_position[0] << " " << source_position[1] << " "
				<< source_position[2] << " " << source_position[3] << std::endl;


		// this is important!!!
		int const sourcePos[3] = {source_position[0], source_position[1], source_position[2]}; 


		std::vector< int* > momenta;
		for(size_t I=0; I<19; I++)
			momenta.push_back( new int[3]);
		{
			int momenta_raw[57] = {
				+0, +0, +0,
				1, +0, +0,
				-1, +0, +0,
				+0,  1, +0,
				+0, -1, +0,
				+0, +0,  1,
				+0, +0, -1,
				1, +0,  1,
				-1, +0, -1,
				1, +0, -1,
				-1, +0,  1,
				1,  1, +0,
				-1, -1, +0,
				1, -1, +0,
				-1,  1, +0,
				+0,  1,  1,
				+0, -1, -1,
				+0,  1, -1,
				+0, -1,  1};




			{

				for(size_t I=0; I<momenta.size(); I++)
					std::copy(&(momenta_raw[3*I]), &(momenta_raw[3*I]) + 3, momenta[I]);
			}

			C_vv_mom=Contract::compute_loop_new(g5_phi_vv,phi,my_operators,sourcePos,momenta,timeslice_source);


			for(size_t i=0; i < C_vv_mom.size(); i++)
			{
				C_vv_mom[i] *=  std::complex<double>(2.0/(16.*kappa*kappa)); 
			}



			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;
		}
#endif



		if (weave.isRoot())	
		{


			std::ofstream fout_vv;
#ifdef _with_momentum_projection
			std::ofstream fout_vv_mom;
#endif
			std::ofstream fout_conserved_vv;

#ifdef _with_Omunu_

			std::ofstream fout_twist2_vv;
			std::ofstream fout_twist2_pol_vv;

#endif

			if (n==0)
			{

				fout_vv.open("output_disc_vvSum.dat");
				//	fout_vv_sum.open("output_disc_vv_sum.dat");

#ifdef _with_momentum_projection
				fout_vv_mom.open("output_disc_vvSum_mom.dat");
#endif
				fout_conserved_vv.open("output_disc_conserved_vvSum.dat");

#ifdef _with_Omunu_
				fout_twist2_vv.open("output_disc_twist2_vvSum.dat");
				fout_twist2_pol_vv.open("output_disc_twist2_pol_vvSum.dat");

#endif
			}
			else	
			{

				fout_vv.open("output_disc_vvSum.dat",std::ios::app);

#ifdef _with_momentum_projection
				fout_vv_mom.open("output_disc_vvSum_mom.dat",std::ios::app);
#endif
				fout_conserved_vv.open("output_disc_conserved_vvSum.dat",std::ios::app);

#ifdef _with_Omunu_
				fout_twist2_vv.open("output_disc_twist2_vvSum.dat",std::ios::app);
				fout_twist2_pol_vv.open("output_disc_twist2_pol_vvSum.dat",std::ios::app);

#endif
			}

			std::cout << "write output for sources from " << 12*n << " to " << 11+12*n<<" ... "  << std::endl;

			for(size_t j=0; j<12; j++)
			{

				for(size_t i=0; i<my_operators.size(); i++)
				{

					for(size_t t = 0; t < T; t++)
					{
										fout_vv << t << std::scientific <<"  "
							<< i <<"  "<< j + 12*n <<"  "<<C_vv[t +  T*i + 16*T*j].real() <<"  "<< C_vv[t + T*i + 16*T*j].imag() << std::endl;



					}
				}			
			}
			fout_vv.close();


			for(size_t j=0; j<12; j++)
			{
				for(size_t i=0; i<4; i++)
				{
					for(size_t k=0; k<2; k++)
					{
						for(size_t t = 0; t < T; t++)
						{

							fout_conserved_vv << t << std::scientific <<"  "
								<< i <<"  "<< j+12*n <<"  "<<k <<"  "<< C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].real() <<"  "<< C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].imag() << std::endl;
						}
					}			
				}
			}

			fout_conserved_vv.close();

#ifdef _with_Omunu_

			for(size_t j=0; j<12; j++)
			{
				for(size_t i=0; i<4; i++)
				{
					for(size_t k=0; k<4; k++)
					{
						for(size_t t = 0; t < T; t++)
						{

							fout_twist2_vv << t << std::scientific <<"  "
								<< i <<"  "<< j + 12*n <<"  "<< k <<"  " <<C_twist2_vv[t  + T*k + 4*T*i + 4*T*4*j ].real() <<"  "<< C_twist2_vv[t + T*k + 4*T*i +4*T*4*j ].imag() << std::endl;
							fout_twist2_pol_vv << t << std::scientific <<"  "
								<< i <<"  "<< j + 12*n <<"  "<< k <<"  " <<C_twist2_pol_vv[t  + T*k + 4*T*i + 4*T*4*j ].real() <<"  "<< C_twist2_pol_vv[t + T*k + 4*T*i +4*T*4*j ].imag() << std::endl;


							fout_twist2_vv.flush();
							fout_twist2_pol_vv.flush();


						}/*t*/
					}/*end loop on the 4 terms contributing to a symmetrized covariant derivative */
				} /*i loop on the 4 operators : O_mumu (unpo1arized) and Omumu g5 (polarized) */
			} /*loop on the 12 sources*/

			fout_twist2_vv.close();
			fout_twist2_pol_vv.close();

#endif

#ifdef _with_momentum_projection
			int N=momenta.size();
			int index;

			for(size_t I = 0; I < momenta.size(); I++)
			{
				for(size_t j=0; j<12; j++)
				{
					for(size_t i=0; i<16; i++)
					{
						for(size_t t = 0; t < T; t++)
						{
							index = t + T*i  +T*16*j + T*16*12*I;

							fout_vv_mom << t << std::scientific <<"  " << momenta[I][0] << "  "<< momenta[I][1] << "  "	<< momenta[I][2] << "  "
								<< i <<"  "<< j + 12*n <<"  "<<C_vv_mom[index].real() <<"  "<< C_vv_mom[index].imag() << std::endl;

						}
					}
				}
			}
			fout_vv_mom.close();
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


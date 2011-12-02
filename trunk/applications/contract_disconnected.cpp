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


//#define _with_Fuzzing_

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

	//read and declare gauge field 
	Core::Field< QCD::Gauge > gauge_field(L, T);


	if (weave.isRoot())
		std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
	Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
	if (weave.isRoot())
		std::cout << "done.\n" << std::endl;



#ifdef _with_APE_
	Core::Field< QCD::Gauge > gauge_field_f(gauge_field);

	double const APE_alpha      = floats["APE_param"];
	size_t const APE_iterations = size_t(floats["APE_steps"]);
	double const Jac_alpha      = floats["Jac_param"];
	size_t const Jac_iterations = size_t(floats["Jac_steps"]);
#endif

#ifdef _with_Fuzzing__
	double const Nlong    = floats["Fuzz_param"];
#endif
	
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

#ifdef _with_APE_
	if(weave.isRoot()) 
		std::cout<<"Smear gauge field   ... "; 


	clock_t start_ape, finish_ape;
	start_ape = clock();


	Smear::APE APE_tool(APE_alpha);
	APE_tool.smear(gauge_field_f, APE_iterations);


	finish_ape = clock();

	if (weave.isRoot())
		std::cout << "APE smearing done in "<< double(finish_ape - start_ape)/CLOCKS_PER_SEC  << "seconds." << std::endl;


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

		Core::StochasticPropagator< 1 > source(L,T);
		Core::StochasticPropagator< 1 > tmp(L, T);

		for(size_t n=0;n<Nsample;n++)
		{  

			std::vector<std::string> filename;
			filename.push_back(propaStochaFiles[n]);

			Core::StochasticPropagator< 1 > phi(L, T);

			if (flag_mms==false)
			{
				if(weave.isRoot()) std::cout<< "Stochastic propagator " <<n<< " to be read from " << propaStochaFiles[n] <<" ... "; 


				Tool::IO::load(&phi,filename, Tool::IO::fileSCIDAC);
				if (weave.isRoot())
					std::cout << "done.\n" << std::endl;
			}
			// if source have been produced using MMS ... read (D+D-)^-1
			if (flag_mms==true)
			{
				if(weave.isRoot())	
					std::cout << "MMS :" << propaStochaFiles[n] << " is 1/(Ddagger D) so apply Ddagger to get the propagator" <<" ...";

				Core::StochasticPropagator< 1 > DD_prop(L, T);
				Core::StochasticPropagator< 1 > prop_out(L, T);

				Tool::IO::load(&DD_prop,filename, Tool::IO::fileSCIDAC,32);

				prop_out=DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz);
				prop_out.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator

				Core::StochasticPropagator< 1 > phi(prop_out);

				if(weave.isRoot())	
					std::cout << "done.\n" << std::endl;
			}
#ifdef _with_APE_

			Core::StochasticPropagator< 1 >  phi_f(phi);
			Smear::Jacobi Jacobi_tool(Jac_alpha);
			phi_f.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_f);
#endif

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


			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;


			//g5 [B^dagger H]^4 g5 computed

			//to have the same normalization than carsten . Origin ?
			xi *= 1./(4.*kappa);

			// for "vv" correlators

			Core::StochasticPropagator< 1 >  psi(phi);
			Core::StochasticPropagator< 1 >  g5_phi(phi);
			g5_phi.rightMultiply(gamma5);
			
#ifdef _with_APE_
			Core::StochasticPropagator< 1 >  g5_phi_f(phi_f);
			g5_phi_f.rightMultiply(gamma5);
#endif
			
			


			if (weave.isRoot())
				std::cout << "Compute loops for vv and v4 method"<<" ... ";
			// v4 & vv loop local quark fields
			std::vector< Core::Correlator< Dirac::Matrix > > C_v4 = Contract::compute_loop(xi,phi,my_operators);
			std::vector< Core::Correlator< Dirac::Matrix > > C_vv = Contract::compute_loop(g5_phi,phi,my_operators);

			// v4 & vv loop smeared quark fields
#ifdef _with_APE_
			std::vector< Core::Correlator< Dirac::Matrix > > C_v4_f = Contract::compute_loop(xi,phi_f,my_operators);
			std::vector< Core::Correlator< Dirac::Matrix > > C_vv_f = Contract::compute_loop(g5_phi,phi_f,my_operators);
#endif


#ifdef	_with_Omunu_
			std::vector< Core::Correlator< Dirac::Matrix > > C_twist2 = Contract::compute_loop_twist2_operator(gauge_field,xi,phi);
#endif





			for(size_t i=0; i<my_operators.size(); i++)
			{
				C_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	

				C_v4[i].deleteField();
				C_vv[i].deleteField();

			}


			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;

			// constant contributin for :
			//	i=0	addimag = 2.*kappa*mu/sqrt(1 + 4.*kappa*kappa*mu*mu)*L*L*L*3.*4.*2.*kappa*2 *4kappa.;
			//	i=8	addreal = 1./sqrt(1 + 4.*kappa*kappa*mu*mu)*L*L*L*3.*4.*2.*kappa*2. *4kappa;
			//output
			if(weave.isRoot()) 
				std::cout<<"Write results"<< " ... ";

			if (weave.isRoot())	
			{

				double addimag =  2.*kappa*mu*L*L*L*3.*4. / (sqrt(1 + 4.*kappa*kappa*mu*mu));
				double addreal = L*L*L*3.*4. / ( sqrt(1 + 4.*kappa*kappa*mu*mu));


				std::ofstream fout_v4;
				std::ofstream fout_vv;
#ifdef _with_Omunu_
				std::ofstream fout_twist2;
#endif
				if (n==0)
				{
					fout_v4.open("output_disc_v4.dat");
					fout_vv.open("output_disc_vv.dat");

#ifdef _with_Omunu_
					fout_twist2.open("output_disc_twist2.dat");
#endif				
				}
				else
				{
					fout_v4.open("output_disc_v4.dat",std::ios::app);
					fout_vv.open("output_disc_vv.dat",std::ios::app);

#ifdef _with_Omunu_
					fout_twist2.open("output_disc_twist2.dat",std::ios::app);
#endif
				}

				for(size_t i=0; i<my_operators.size(); i++)
				{
					for(size_t t = 0; t < T; t++)
					{

						if (i != 0 && i !=8)
						{
							fout_v4 << t << std::scientific <<"  "
								<< i <<"  "<< n <<"  "<<C_v4[i][t].trace().real() <<"  "<< C_v4[i][t].trace().imag() << std::endl;
						}
						if (i==0) 
						{
							fout_v4 << t << std::scientific <<"  "
								<< i <<"  "<< n <<"  "<<C_v4[i][t].trace().real() <<"  "<< C_v4[i][t].trace().imag() +  addimag << std::endl;
						}
						if (i==8) 
						{
							fout_v4 << t << std::scientific <<"  "
								<< i <<"  "<< n <<"  "<<C_v4[i][t].trace().real()  + addreal <<"  "<< C_v4[i][t].trace().imag()  << std::endl;
						}

						fout_vv << t << std::scientific <<"  "
							<< i <<"  "<< n <<"  "<<C_vv[i][t].trace().real() <<"  "<< C_vv[i][t].trace().imag() << std::endl;

					}			
				}

#ifdef _with_Omunu_
				for(size_t i=0; i<8; i++)
				{
					for(size_t j=0; j<4; j++)
					{
						for(size_t t = 0; t < T; t++)
						{

							fout_twist2 << t << std::scientific <<"  "
								<< i <<"  "<< j <<"  "<< n <<"  "<<C_twist2[4*i+j][t].trace().real() <<"  "<< C_twist2[4*i+j][t].trace().imag() << std::endl;

						}/*end loop on the 4 terms contributing to a symmetrized covariant derivative */
					} /*t*/
				} /*loop on the 8 operator : O_mumu (unpolarized) and Omumu g5 (polarized)*/

				fout_twist2.close();
#endif
				fout_v4.close();
				fout_vv.close();



			}

			C_vv.clear();
			C_v4.clear();

			if (weave.isRoot())

				std::cout << "done.\n" << std::endl;

		}

	}
	else { //  Nsample==12

		if(weave.isRoot()) std::cout<< "The code will read sources by packet of 12 " <<" ... ";

		size_t N = size_t(Nsample/12);
		for (size_t n=0;n<N;n++)
		{

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

			if (weave.isRoot())
				std::cout << "\nReading done.\n" << std::endl;

			// for "vv" correlators

			clock_t start, finish;
			start = clock();


			Core::Propagator  g5_phi(phi);
			g5_phi.rightMultiply(gamma5);



			//Now compute the source x gamma5 
			std::vector< std::complex <double>  > C_naive;
			std::vector< std::complex <double>  > C_v4;
			std::vector< std::complex <double> > C_conserved_v4;

#ifdef _with_APE_
			std::vector< std::complex <double> > C_v4_f;
			std::vector< std::complex <double> > C_vv_f;
#endif

#ifdef _with_Omunu_
			std::vector< std::complex<double> > C_twist2;
			std::vector< std::complex<double> > C_twist2_pol;
#endif

			{

#ifdef _with_APE_
				clock_t start_Jac, finish_Jac;
				start_Jac = clock();


				Core::Propagator  phi_f(phi);
				Smear::Jacobi Jacobi_tool(Jac_alpha);
				phi_f.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_f);


				finish_Jac = clock();

				if (weave.isRoot())
					std::cout << "Jacobi smearing done in "<< double(finish_Jac - start_Jac)/CLOCKS_PER_SEC  << "seconds." << std::endl;

#endif

				{
					if(weave.isRoot()) 
						std::cout << "Apply Dirac operator to get the source and multiply by gamma_5 " << " ... ";


					Core::Propagator  xi(L,T);
					xi = phi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Full); 



					C_naive = Contract::compute_loop_new(xi,phi,my_operators);
					
					// to check					
					if (weave.isRoot())
						std::cout << "done.\n" << std::endl;

					{
						double const xiNorm = xi.norm();
						if (weave.isRoot())
							std::cout << "norm of source calculated from  propagator: " << xiNorm << std::endl;

						double const phiNorm = phi.norm();
						if (weave.isRoot())
							std::cout << "norm of the propagator: " << phiNorm << std::endl;


					}


					xi.rightMultiply(gamma5);

					// Now compute g5 [B^dagger H]^4 g5
					if(weave.isRoot()) 
						std::cout<<"Compute g5 [B^dagger H]^4 g5 times the source field "<< " ... ";


					clock_t start_h4, finish_h4;
					start_h4 = clock();



					for(size_t i=0;i<4;i++)
					{	
						/* apply Hopping part ...*/
						xi = xi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::H);
						/* apply Bdagger ... */
						xi = xi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Bdagger);				
					}


					finish_h4 = clock();

					if (weave.isRoot())
						std::cout << "part H4 done in "<< double(finish_h4 - start_h4)/CLOCKS_PER_SEC  << "seconds." << std::endl;



					xi.rightMultiply(gamma5);


					if (weave.isRoot())
						std::cout << "done.\n" << std::endl;

					//to have the same normalization than carsten . Origin ?
					xi *= 1./(4.*kappa);

					finish = clock();

					if (weave.isRoot())
						std::cout << "part 0 before contract  done in "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;


					// v4  local quark fields
					C_v4 = Contract::compute_loop_new(xi,phi,my_operators);
					C_conserved_v4 = Contract::compute_loop_conserved_vector_current(gauge_field,xi,phi);

					//v4  smeared
					
#ifdef _with_APE_
					C_v4_f = Contract::compute_loop_new(xi,phi_f,my_operators);
#endif
					
#ifdef _with_Omunu_
					clock_t start_twist2, finish_twist2;
					start_twist2 = clock();


					C_twist2 = Contract::compute_loop_twist2_operator(gauge_field,xi,phi);
					C_twist2_pol = Contract::compute_loop_twist2_operator(gauge_field,xi,g5_phi);
					finish_twist2 = clock();

					if (weave.isRoot())
						std::cout << "part twist2 v4  done in "<< double(finish_twist2 - start_twist2)/CLOCKS_PER_SEC  << "seconds." << std::endl;
#endif

				}

				// vv loop smeared quark fields

#ifdef _with_APE_
				C_vv_f = Contract::compute_loop_new(g5_phi,phi_f,my_operators);
#endif
			}
			finish = clock();

			if (weave.isRoot())
				std::cout << "part 1 done in "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;



			if (weave.isRoot())
				std::cout << "Compute loops for vv and v4 method"<<" ... ";


			Core::Propagator g5_DW_phi(L,T);
			g5_DW_phi=phi.applyDiracOperator(gauge_field,kappa,0,thetat,thetax,thetay,thetaz); //D_Wilson = D_tm(mu=0)
			g5_DW_phi.rightMultiply(gamma5);

			std::vector< std::complex <double>  > C_vv = Contract::compute_loop_new(g5_phi,phi,my_operators);
			std::vector< std::complex <double>  > C_vvSum = Contract::compute_loop_new(g5_phi,g5_DW_phi,my_operators);

			std::vector< std::complex <double> > C_conserved_vv = Contract::compute_loop_conserved_vector_current(gauge_field,g5_phi,phi);


			for(size_t i=0; i < C_vv.size(); i++)
			{
				C_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
				C_vvSum[i] *=  std::complex<double>(2.0,0);	 
#ifdef _with_APE
				C_vv_f[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
#endif

			}

			for(size_t i=0; i < C_conserved_vv.size(); i++)
			{
				C_conserved_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
			}

			if (weave.isRoot()) std::cout << "done.\n" << std::endl;

			finish = clock();

			if (weave.isRoot())
				std::cout << "part 2 done in "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;


#ifdef _with_Omunu_
			if (weave.isRoot())
				std::cout << "Twist 2 operators  "<<" ... ";


			std::vector< std::complex<double> > C_twist2_vv = Contract::compute_loop_twist2_operator(gauge_field,g5_phi,phi,false);
			std::vector< std::complex<double> > C_twist2_pol_vv = Contract::compute_loop_twist2_operator(gauge_field,g5_phi,phi,true);
			std::vector< std::complex<double> > C_twist2_vvSum = Contract::compute_loop_twist2_operator(gauge_field,g5_phi,g5_DW_phi,false);
			std::vector< std::complex<double> > C_twist2_pol_vvSum = Contract::compute_loop_twist2_operator(gauge_field,g5_phi,g5_DW_phi,true);


			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;


#endif


			finish = clock();

			if (weave.isRoot())
				std::cout << "part 3 done in "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;




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

				C_vv_mom=Contract::compute_loop_new(g5_phi,phi,my_operators,sourcePos,momenta,timeslice_source);


				for(size_t i=0; i < C_vv_mom.size(); i++)
				{
					C_vv_mom[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa)); 
				}


				if (weave.isRoot())
					std::cout << "done  "<<  std::endl;
			}
#endif

			finish = clock();

			if (weave.isRoot())
				std::cout << "Computation of loops done in "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;


			if (weave.isRoot())	
			{
				double addimag =0; 
				double addreal =0;

				//double addimag =  2.*kappa*mu*L*L*L*3.*4. / (sqrt(1 + 4.*kappa*kappa*mu*mu));
				//double addreal =  L*L*L*3.*4./( sqrt(1. + 4.*kappa*kappa*mu*mu));

				FILE * fout_v4;
				FILE * fout_naive;
				//				am fout_vv;
				FILE * fout_vv;
				FILE * fout_vvSum;
#ifdef _with_APE
				FILE * fout_v4_f;
				FILE * fout_vv_f;
#endif

#ifdef _with_momentum_projection
				FILE * fout_vv_mom;
#endif
				FILE * fout_conserved_v4;
				FILE * fout_conserved_vv;

#ifdef _with_Omunu_
				FILE * fout_twist2;
				FILE * fout_twist2_pol;

				FILE * fout_twist2_vvSum;
				FILE * fout_twist2_pol_vvSum;
				FILE * fout_twist2_vv;
				FILE * fout_twist2_pol_vv;

#endif

				if (n==0)
				{

					fout_v4=fopen("output_disc_v4.dat","w");
					fout_naive=fopen("output_disc_naive.dat","w");
					fout_vvSum=fopen("output_disc_vvSum.dat","w");
					fout_vv =fopen("output_disc_vv.dat","w");
#ifdef _with_APE
					fout_v4_f=fopen("output_disc_v4_f.dat","w");
					fout_vv_f=fopen("output_disc_vv_f.dat","w");
#endif					

#ifdef _with_momentum_projection
					fout_vv_mom=fopen("output_disc_vv_mom.dat","w");
#endif
					fout_conserved_v4=fopen("output_disc_conserved_v4.dat","w");
					fout_conserved_vv=fopen("output_disc_conserved_vv.dat","w");

#ifdef _with_Omunu_
					fout_twist2=fopen("output_disc_twist2.dat","w");
					fout_twist2_pol=fopen("output_disc_twist2_pol.dat","w");
					fout_twist2_vvSum=fopen("output_disc_twist2_vvSum.dat","w");
					fout_twist2_pol_vvSum=fopen("output_disc_twist2_pol_vvSum.dat","w");
					fout_twist2_vv=fopen("output_disc_twist2_vv.dat","w");
					fout_twist2_pol_vv=fopen("output_disc_twist2_pol_vv.dat","w");


#endif
				}
				else	
				{

					fout_v4=fopen("output_disc_v4.dat","a");
					fout_naive=fopen("output_disc_naive.dat","a");
					fout_vvSum=fopen("output_disc_vvSum.dat","a");
					fout_vv =fopen("output_disc_vv.dat","a");
					//	fout_vv_sum.open("output_disc_vv_sum.dat",std::ios::app);
#ifdef _with_APE
					fout_v4_f=fopen("output_disc_v4_f.dat","a");
					fout_vv_f=fopen("output_disc_vv_f.dat","a");
#endif

#ifdef _with_momentum_projection
					fout_vv_mom=fopen("output_disc_vv_mom.dat","a");
#endif
					fout_conserved_v4=fopen("output_disc_conserved_v4.dat","a");
					fout_conserved_vv=fopen("output_disc_conserved_vv.dat","a");

#ifdef _with_Omunu_
					fout_twist2=fopen("output_disc_twist2.dat","a");
					fout_twist2_pol=fopen("output_disc_twist2_pol.dat","a");
					fout_twist2_vv=fopen("output_disc_twist2_vv.dat","a");
					fout_twist2_pol_vv=fopen("output_disc_twist2_pol_vv.dat","a");
					fout_twist2_vvSum=fopen("output_disc_twist2_vvSum.dat","a");
					fout_twist2_pol_vvSum=fopen("output_disc_twist2_pol_vvSum.dat","a");


#endif
				}

				std::cout << "write output for sources from " << 12*n << " to " << 11+12*n<<" ... "  << std::endl;

				//	clock_t start, finish;
				start = clock();

				for(size_t j=0; j<12; j++)
				{

					for(size_t i=0; i<my_operators.size(); i++)
					{

						for(size_t t = 0; t < T; t++)
						{
							if (i != 0 && i !=8)
							{
								fprintf(fout_v4,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_v4[t +  T*i + 16*T*j].real(),C_v4[t + T*i + 16*T*j].imag() );
								fprintf(fout_naive,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_naive[t +  T*i + 16*T*j].real(),C_naive[t + T*i + 16*T*j].imag() );
#ifdef _with_APE	
								fprintf(fout_v4_f,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_v4_f[t +  T*i + 16*T*j].real(),C_v4_f[t + T*i + 16*T*j].imag() );
#endif


							}
							if (i==0) 
							{
								fprintf(fout_v4,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_v4[t +  T*i + 16*T*j].real(),C_v4[t + T*i + 16*T*j].imag() +  addimag );
#ifdef _with_APE
								fprintf(fout_v4_f,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_v4_f[t +  T*i + 16*T*j].real(),C_v4_f[t + T*i + 16*T*j].imag() +  addimag );
#endif

							}
							if (i==8) 
							{
								fprintf(fout_v4,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_v4[t +  T*i + 16*T*j].real()+ addreal,C_v4[t + T*i + 16*T*j].imag() );
#ifdef _with_APE
								fprintf(fout_v4_f,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_v4_f[t +  T*i + 16*T*j].real()+ addreal,C_v4_f[t + T*i + 16*T*j].imag() );
#endif

							}


							fprintf(fout_vv,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_vv[t +  T*i + 16*T*j].real(),C_vv[t + T*i + 16*T*j].imag() );
							fprintf(fout_vvSum,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_vvSum[t +  T*i + 16*T*j].real(),C_vvSum[t + T*i + 16*T*j].imag() );
#ifdef _with_APE
							fprintf(fout_vv_f,"%3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n, C_vv_f[t +  T*i + 16*T*j].real(),C_vv_f[t + T*i + 16*T*j].imag() );
#endif

						}
					}			
				}
				fclose(fout_v4);				//				fout_vv.close();
				fclose(fout_vv);
				fclose(fout_vvSum);
				//	fout_vv_sum.close();
#ifdef _with_APE
				fclose(fout_v4_f);
				fclose(fout_vv_f);
#endif

				for(size_t j=0; j<12; j++)
				{
					for(size_t i=0; i<4; i++)
					{
						for(size_t k=0; k<2; k++)
						{
							for(size_t t = 0; t < T; t++)
							{
								fprintf(fout_conserved_v4,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k, C_conserved_v4[t+ T*k + 2*T*i + 2*T*4*j].real(),C_conserved_v4[t+ T*k + 2*T*i + 2*T*4*j].imag() );
								fprintf(fout_conserved_vv,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k, C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].real(),C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].imag() );


								//								fout_conserved_v4 << t << std::scientific <<"  "
								//									<< i <<"  "<< j+12*n <<"  "<<k <<"  "<< C_conserved_v4[t+ T*k + 2*T*i + 2*T*4*j].real() <<"  "<< C_conserved_v4[t+ T*k + 2*T*i + 2*T*4*j].imag() << std::endl;

								//								fout_conserved_vv << t << std::scientific <<"  "
								//									<< i <<"  "<< j+12*n <<"  "<<k <<"  "<< C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].real() <<"  "<< C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].imag() << std::endl;
							}
						}			
					}
				}

				fclose(fout_conserved_vv);
				fclose(fout_conserved_v4);

#ifdef _with_Omunu_

				for(size_t j=0; j<12; j++)
				{
					for(size_t i=0; i<4; i++)
					{
						for(size_t k=0; k<4; k++)
						{
							for(size_t t = 0; t < T; t++)
							{


								fprintf(fout_twist2,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k,C_twist2[t  + T*k + 4*T*i + 4*T*4*j ].real(),C_twist2[t + T*k + 4*T*i +4*T*4*j ].imag());
								fprintf(fout_twist2_pol,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k,C_twist2_pol[t  + T*k + 4*T*i + 4*T*4*j ].real(),C_twist2_pol[t + T*k + 4*T*i +4*T*4*j ].imag());
								fprintf(fout_twist2_vv,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k,C_twist2_vv[t  + T*k + 4*T*i + 4*T*4*j ].real(),C_twist2_vv[t + T*k + 4*T*i +4*T*4*j ].imag());
								fprintf(fout_twist2_pol_vv,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k,C_twist2_pol_vv[t  + T*k + 4*T*i + 4*T*4*j ].real(),C_twist2_pol_vv[t + T*k + 4*T*i +4*T*4*j ].imag());
								fprintf(fout_twist2_vvSum,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k,C_twist2_vvSum[t  + T*k + 4*T*i + 4*T*4*j ].real(),C_twist2_vvSum[t + T*k + 4*T*i +4*T*4*j ].imag());
								fprintf(fout_twist2_pol_vvSum,"%3d %3d %3d %3d %3.10e %3.10e\n",t,i,j+12*n,k,C_twist2_pol_vvSum[t  + T*k + 4*T*i + 4*T*4*j ].real(),C_twist2_pol_vvSum[t + T*k + 4*T*i +4*T*4*j ].imag());

							}/*t*/
						}/*end loop on the 4 terms contributing to a symmetrized covariant derivative */
					} /*i loop on the 4 operators : O_mumu (unpo1arized) and Omumu g5 (polarized) */
				} /*loop on the 12 sources*/

				fclose(fout_twist2);
				fclose(fout_twist2_pol);
				fclose(fout_twist2_vv);
				fclose(fout_twist2_pol_vv);
				fclose(fout_twist2_vvSum);
				fclose(fout_twist2_pol_vvSum);

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
								fprintf(fout_vv_mom,"%3d %3d %3d %3d %3d %3d %3.10e %3.10e \n",t,momenta[I][0],momenta[I][1],momenta[I][2],i,j+12*n,C_vv_mom[index].real(),C_vv_mom[index].imag());
								//								fout_vv_mom << t << std::scientific <<"  " << momenta[I][0] << "  "<< momenta[I][1] << "  "	<< momenta[I][2] << "  "
								//									<< i <<"  "<< j + 12*n <<"  "<<C_vv_mom[index].real() <<"  "<< C_vv_mom[index].imag() << std::endl;

							}
						}
					}
				}
				fclose(fout_vv_mom);
#endif				
				finish = clock();
				std::cout << "done in "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;

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



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
#define _with_MMS_
//#define _without_MMS_
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

		Core::StochasticPropagator< 1 > source(L,T);
		Core::StochasticPropagator< 1 > tmp(L, T);

		for(size_t n=0;n<Nsample;n++)
		{  

			std::vector<std::string> filename;
			filename.push_back(propaStochaFiles[n]);


#ifdef _without_MMS_
			if(weave.isRoot()) std::cout<< "Stochastic propagator " <<n<< " to be read from " << propaStochaFiles[n] <<" ... "; 

			Core::StochasticPropagator< 1 > phi(L, T);
			Tool::IO::load(&phi,filename, Tool::IO::fileSCIDAC);
			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;
#endif
			// if source have been produced using MMS ... read (D+D-)^-1
#ifdef _with_MMS_
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


			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;


			//g5 [B^dagger H]^4 g5 computed

			//to have the same normalization than carsten . Origin ?
			xi *= 1./(4.*kappa);

			// for "vv" correlators

			Core::StochasticPropagator< 1 >  psi(phi);
			Core::StochasticPropagator< 1 >  g5_phi(phi);
			Core::StochasticPropagator< 1 >  g5_phi_f(phi_f);

			g5_phi.rightMultiply(gamma5);
			g5_phi_f.rightMultiply(gamma5);


			if (weave.isRoot())
				std::cout << "Compute loops for vv and v4 method"<<" ... ";
			// v4 & vv loop local quark fields
			std::vector< Core::Correlator< Dirac::Matrix > > C_v4 = Contract::compute_loop(xi,phi,my_operators);
			std::vector< Core::Correlator< Dirac::Matrix > > C_vv = Contract::compute_loop(g5_phi,phi,my_operators);

			// v4 & vv loop smeared quark fields
			std::vector< Core::Correlator< Dirac::Matrix > > C_v4_f = Contract::compute_loop(xi,phi_f,my_operators);
			std::vector< Core::Correlator< Dirac::Matrix > > C_vv_f = Contract::compute_loop(g5_phi,phi_f,my_operators);



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


#ifdef _without_MMS_
			Core::Propagator phi(L, T);
			Tool::IO::load(&phi,filename, Tool::IO::fileSCIDAC);
#endif
			// if source have been produced using MMS ... read (D+D-)^-1
#ifdef _with_MMS_

			if(weave.isRoot()) 
				std::cout << "MMS compilation  1/(Ddagger D) is read so  Ddagger is applied to get the propagator" <<" ...";

			Core::Propagator phi(L ,T );

			{
			//	Core::Propagator prop_out(L, T);
				Core::Propagator DD_prop(L, T);
				Tool::IO::load(&DD_prop,filename, Tool::IO::fileSCIDAC,32);

				phi=DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz);
				phi.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator

//				phi = prop_out;
			}
#endif

			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;


				// for "vv" correlators

				Core::Propagator  g5_phi(phi);
				g5_phi.rightMultiply(gamma5);



			//Now compute the source x gamma5 
			std::vector< std::complex <double>  > C_v4;
			std::vector< std::complex <double> > C_conserved_v4;

			std::vector< std::complex <double> > C_v4_f;
			std::vector< std::complex <double> > C_vv_f;


#ifdef _with_Omunu_
			std::vector< std::complex<double> > C_twist2;
			std::vector< std::complex<double> > C_twist2_pol;
#endif

			{
				Core::Propagator  phi_f(phi);
				Smear::Jacobi Jacobi_tool(Jac_alpha);
				phi_f.smearJacobi(Jac_alpha, Jac_iterations, gauge_field_f);



				{
					if(weave.isRoot()) 
						std::cout << "Apply Dirac operator to get the source and multiply by gamma_5 " << " ... ";


					Core::Propagator  xi(L,T);
					xi = phi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Full); 
					xi.rightMultiply(gamma5);

					if (weave.isRoot())
						std::cout << "done.\n" << std::endl;

					// Now compute g5 [B^dagger H]^4 g5
					if(weave.isRoot()) 
						std::cout<<"Compute g5 [B^dagger H]^4 g5 times the source field "<< " ... ";


					for(size_t i=0;i<4;i++)
					{	
						/* apply Hopping part ...*/
						xi = xi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::H);
						/* apply Bdagger ... */
						xi = xi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Bdagger);				
					}
					xi.rightMultiply(gamma5);


					if (weave.isRoot())
						std::cout << "done.\n" << std::endl;

					//to have the same normalization than carsten . Origin ?
					xi *= 1./(4.*kappa);

					// v4  local quark fields
					C_v4 = Contract::compute_loop_new(xi,phi,my_operators);
					C_conserved_v4 = Contract::compute_loop_conserved_vector_current(gauge_field,xi,phi);

					//v4  smeared
					C_v4_f = Contract::compute_loop_new(xi,phi_f,my_operators);

#ifdef _with_Omunu_
					C_twist2 = Contract::compute_loop_twist2_operator(gauge_field,xi,phi);
					C_twist2_pol = Contract::compute_loop_twist2_operator(gauge_field,xi,g5_phi);
#endif


				}

				// vv loop smeared quark fields
				C_vv_f = Contract::compute_loop_new(g5_phi,phi_f,my_operators);
			}


			if (weave.isRoot())
				std::cout << "Compute loops for vv and v4 method"<<" ... ";
			std::vector< std::complex <double>  > C_vv = Contract::compute_loop_new(g5_phi,phi,my_operators);
			std::vector< std::complex <double> > C_conserved_vv = Contract::compute_loop_conserved_vector_current(gauge_field,g5_phi,phi);



			for(size_t i=0; i < C_vv.size(); i++)
			{
				C_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
				C_vv_f[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
			}

			for(size_t i=0; i < C_conserved_vv.size(); i++)
			{
				C_conserved_vv[i] *=  std::complex<double>(0,2.0*mu/(8.*kappa));	 // factor + 4 * i kappa mu /( 4 kappa) ^2
			}

			if (weave.isRoot())
				std::cout << "done.\n" << std::endl;



#ifdef _with_Omunu_
			if (weave.isRoot())
				std::cout << "Twist 2 operators  "<<" ... ";


			std::vector< std::complex<double> > C_twist2_vv = Contract::compute_loop_twist2_operator(gauge_field,g5_phi,phi,false);
			std::vector< std::complex<double> > C_twist2_pol_vv = Contract::compute_loop_twist2_operator(gauge_field,g5_phi,phi,true);

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
				std::cout << "\nsource position: " << source_position[0] << " " << source_position[1] << " "
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
					std::cout << "done.\n" << std::endl;
			}
#endif



			if (weave.isRoot())	
			{
				double addimag =0; 
				double addreal =0;

				//double addimag =  2.*kappa*mu*L*L*L*3.*4. / (sqrt(1 + 4.*kappa*kappa*mu*mu));
				//double addreal =  L*L*L*3.*4./( sqrt(1. + 4.*kappa*kappa*mu*mu));

				std::ofstream fout_v4;
				std::ofstream fout_vv;

				std::ofstream fout_v4_f;
				std::ofstream fout_vv_f;


#ifdef _with_momentum_projection
				std::ofstream fout_vv_mom;
#endif
				std::ofstream fout_conserved_v4;
				std::ofstream fout_conserved_vv;

#ifdef _with_Omunu_
				std::ofstream fout_twist2;
				std::ofstream fout_twist2_pol;

				std::ofstream fout_twist2_vv;
				std::ofstream fout_twist2_pol_vv;

#endif

				if (n==0)
				{

					fout_v4.open("output_disc_v4.dat");
					fout_vv.open("output_disc_vv.dat");
					fout_v4_f.open("output_disc_v4_f.dat");
					fout_vv_f.open("output_disc_vv_f.dat");

#ifdef _with_momentum_projection
					fout_vv_mom.open("output_disc_vv_mom.dat");
#endif
					fout_conserved_v4.open("output_disc_conserved_v4.dat");
					fout_conserved_vv.open("output_disc_conserved_vv.dat");

#ifdef _with_Omunu_
					fout_twist2.open("output_disc_twist2.dat");
					fout_twist2_pol.open("output_disc_twist2_pol.dat");
					fout_twist2_vv.open("output_disc_twist2_vv.dat");
					fout_twist2_pol_vv.open("output_disc_twist2_pol_vv.dat");

#endif
				}
				else	
				{

					fout_v4.open("output_disc_v4.dat",std::ios::app);
					fout_vv.open("output_disc_vv.dat",std::ios::app);
					fout_v4_f.open("output_disc_v4_f.dat",std::ios::app);
					fout_vv_f.open("output_disc_vv_f.dat",std::ios::app);

#ifdef _with_momentum_projection
					fout_vv_mom.open("output_disc_vv_mom.dat",std::ios::app);
#endif
					fout_conserved_v4.open("output_disc_conserved_v4.dat",std::ios::app);
					fout_conserved_vv.open("output_disc_conserved_vv.dat",std::ios::app);

#ifdef _with_Omunu_
					fout_twist2.open("output_disc_twist2.dat",std::ios::app);
					fout_twist2_pol.open("output_disc_twist2_pol.dat",std::ios::app);
					fout_twist2_vv.open("output_disc_twist2_vv.dat",std::ios::app);
					fout_twist2_pol_vv.open("output_disc_twist2_pol_vv.dat",std::ios::app);

#endif
				}

				std::cout << "write output for sources from " << 12*n << " to " << 11+12*n<<" ... "  << std::endl;

				for(size_t j=0; j<12; j++)
				{

					for(size_t i=0; i<my_operators.size(); i++)
					{

						for(size_t t = 0; t < T; t++)
						{
							if (i != 0 && i !=8)
							{
								fout_v4 << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<<C_v4[t +  t*i + 16*t*j].real() <<"  "<< C_v4[t + t*i + 16*t*j].imag() << std::endl;
								fout_v4_f << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<<C_v4_f[t +  t*i + 16*t*j].real() <<"  "<< C_v4_f[t + t*i + 16*t*j].imag() << std::endl;
							}
							if (i==0) 
							{
								fout_v4 << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<<C_v4[t +  T*i + 16*T*j].real() <<"  "<< C_v4[t + T*i + 16*T*j].imag() +  addimag << std::endl;
								fout_v4_f << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<<C_v4_f[t +  T*i + 16*T*j].real() <<"  "<< C_v4_f[t + T*i + 16*T*j].imag() +  addimag << std::endl;
							}
							if (i==8) 
							{
								fout_v4 << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<<C_v4[t +  T*i + 16*T*j].real() + addreal <<"  "<< C_v4[t + T*i + 16*T*j].imag() << std::endl;
								fout_v4_f << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<<C_v4_f[t +  T*i + 16*T*j].real() + addreal <<"  "<< C_v4_f[t + T*i + 16*T*j].imag() << std::endl;
							}


							fout_vv << t << std::scientific <<"  "
								<< i <<"  "<< j + 12*n <<"  "<<C_vv[t +  T*i + 16*T*j].real() <<"  "<< C_vv[t + T*i + 16*T*j].imag() << std::endl;
							fout_vv_f << t << std::scientific <<"  "
								<< i <<"  "<< j + 12*n <<"  "<<C_vv_f[t +  T*i + 16*T*j].real() <<"  "<< C_vv_f[t + T*i + 16*T*j].imag() << std::endl;


							fout_vv_f.flush();
							fout_v4_f.flush();

							fout_vv.flush();
							fout_v4.flush();

						}
					}			
				}
				fout_v4.close();
				fout_vv.close();

				fout_v4_f.close();
				fout_vv_f.close();

				for(size_t j=0; j<12; j++)
				{
					for(size_t i=0; i<4; i++)
					{
						for(size_t k=0; k<2; k++)
						{
							for(size_t t = 0; t < T; t++)
							{
								fout_conserved_v4 << t << std::scientific <<"  "
									<< i <<"  "<< j+12*n <<"  "<<k <<"  "<< C_conserved_v4[t+ T*k + 2*T*i + 2*T*4*j].real() <<"  "<< C_conserved_v4[t+ T*k + 2*T*i + 2*T*4*j].imag() << std::endl;

								fout_conserved_vv << t << std::scientific <<"  "
									<< i <<"  "<< j+12*n <<"  "<<k <<"  "<< C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].real() <<"  "<< C_conserved_vv[t+ T*k + 2*T*i + 2*T*4*j].imag() << std::endl;
							}
						}			
					}
				}

				fout_conserved_vv.close();
				fout_conserved_v4.close();

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

								fout_twist2_vv << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<< k <<"  " <<C_twist2_vv[t  + T*k + 4*T*i + 4*T*4*j ].real() <<"  "<< C_twist2_vv[t + T*k + 4*T*i +4*T*4*j ].imag() << std::endl;
								fout_twist2_pol_vv << t << std::scientific <<"  "
									<< i <<"  "<< j + 12*n <<"  "<< k <<"  " <<C_twist2_pol_vv[t  + T*k + 4*T*i + 4*T*4*j ].real() <<"  "<< C_twist2_pol_vv[t + T*k + 4*T*i +4*T*4*j ].imag() << std::endl;


								fout_twist2.flush();
								fout_twist2_pol.flush();
								fout_twist2_vv.flush();
								fout_twist2_pol_vv.flush();


							}/*t*/
						}/*end loop on the 4 terms contributing to a symmetrized covariant derivative */
					} /*i loop on the 4 operators : O_mumu (unpo1arized) and Omumu g5 (polarized) */
				} /*loop on the 12 sources*/

				fout_twist2.close();
				fout_twist2_pol.close();
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



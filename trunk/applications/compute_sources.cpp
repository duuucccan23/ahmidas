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

	Input::FileReader reader("./compute_sources_input.xml");

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
    bool gluon_loop=(bool)floats["gluon_loop"]; // 0=no,!=0 =yes

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
/*   Dirac::Gamma<1> gamma1;
   Dirac::Gamma<2> gamma2;
   Dirac::Gamma<3> gamma3;
   Dirac::Gamma<4> gamma4;

   size_t const timeslice_boundary(T - 1);*/

   //read the gauge configuration, needed to create u and d propagators


   std::vector< std::string > const &propaStochaFiles(files[0]);
   std::vector< std::string > const &gaugeFieldFiles(files[1]); 
	 std::vector< std::string > const & stochasticSourceFiles(files[2]); 
   
   //read and declare gauge field 
   Core::Field< QCD::Gauge > gauge_field(L, T);


   if (weave.isRoot())
     std::cout << "gauge field to be read from " << gaugeFieldFiles[0] << " ... ";
   Tool::IO::load(&gauge_field, gaugeFieldFiles[0], Tool::IO::fileILDG);
   if (weave.isRoot())
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

       //Now compute the source 

       if(weave.isRoot()) 
         std::cout << "Apply Dirac operator to get the source " << " ... ";

       source = phi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Full); 

       if (weave.isRoot())
         std::cout << "done.\n" << std::endl;


       if(weave.isRoot()) 
         std::cout << "Save the source " << " ... ";


       clock_t start_save,finish_save;
       start_save = clock();

       std::vector<std::string> tmp_filename;
       tmp_filename.push_back(stochasticSourceFiles[n]);


       Tool::IO::save(reinterpret_cast< Core::StochasticPropagator< 1 > * >(&source), tmp_filename, Tool::IO::fileSCIDAC);

       finish_save = clock();

       if (weave.isRoot())
         std::cout << " sources saved in "<< double(finish_save - start_save)/CLOCKS_PER_SEC  << "seconds." << std::endl;
     }

   }
   else { //  Nsample==12

     if(weave.isRoot()) std::cout<< "The code will read sources by packet of 12 " <<" ... ";

     Core::Propagator source(L, T);
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

       // now apply the tm Dirac operator to get the source and save it

       source = phi.applyDiracOperator(gauge_field,kappa,mu,thetat,thetax,thetay,thetaz,Base::Full);

       std::vector<std::string> tmp_filename;
       for (size_t k=0; k<12;k++) tmp_filename.push_back(stochasticSourceFiles[k+12*n]);

       if(weave.isRoot()) 
         std::cout << "Save the sources from " << 12*n <<" to " << 12+12*n << " ... ";

       clock_t start_save,finish_save;
       start_save = clock();

       Tool::IO::save(&source,tmp_filename, Tool::IO::fileSCIDAC);

       finish_save = clock();

       if (weave.isRoot())
         std::cout << " sources saved in "<< double(finish_save - start_save)/CLOCKS_PER_SEC  << "seconds." << std::endl;
 
     } /* loop on Nsample */
   }/* end case N%12 = 0 */

#ifdef __MPI_ARCH__
   if (myid == 0)
#endif

     if(weave.isRoot())
       std::cout << "\nprogram is going to exit normally now\n" << std::endl;

   return EXIT_SUCCESS;

}



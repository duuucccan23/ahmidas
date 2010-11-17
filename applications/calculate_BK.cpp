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
  std::vector< int > rcombinations;
  std::vector< std::vector< std::string > > files;
  std::string outputname;  

  reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators,rcombinations,outputname);

  size_t const L(L_tmp);
  size_t const T(T_tmp);

  Base::Weave weave(L, T);


  double kappa = floats["kappa"];
  double mu_d    = floats["mu_d"];
  double mu_s    = floats["mu_s"];
  bool const phys_base = bool(floats["phys_base"] != 0.0); // 0=no,!=0 =yes

  int connected = floats["connected"];  
  bool debug=(bool)floats["debug"];

  size_t const t_L = size_t(floats["t_L"]);
  assert(t_L >= 0 && t_L < T);
  size_t const t_R = size_t(floats["t_R"]);
  assert(t_R >= 0 && t_R < T);
  assert(rcombinations.size() % 4 == 0);


  {
    // write some output concerning what is going to be read
    std::ostringstream out;
    out << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;
    out << "kappa = "  << kappa  << ", mu_d = "<< mu_d << ", mu_s = "<< mu_s << std::endl;
    out << "timeslice (left)  = " << t_L << std::endl;
    out << "timeslice (right) = " << t_R << std::endl;
    out << "will ";
    if (!phys_base) out << "not ";
    out << "rotate propagators to physical base" << std::endl;
    Print(out.str());
  }


 if(weave.isRoot() and debug==true)
    {
      std::cout<<"The following files are going to be read:"<<std::endl;

      for(size_t i=0;i<files.size();i++)
        for(size_t j=0;j<files[i].size();j++)
          std::cout<<(files[i])[j]<<std::endl;

    }




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

  Core::StochasticPropagator< 4 >  *dProp_L[2];
  dProp_L[0] = new Core::StochasticPropagator< 4 >(L, T);
  dProp_L[1] = new Core::StochasticPropagator< 4 >(L, T);
  // reconstruct the full doublet from the combined MMS output, being ~ (DD^dagger)^-1
  tmpProp->reconstruct_doublet(*(dProp_L[0]),*(dProp_L[1]), gauge_field, kappa, mu_d);
  if (weave.isRoot())
    std::cout << "1st propagator successfully reconstructed\n" << std::endl;

  Tool::IO::load(tmpProp, propfilesS_L, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
    std::cout << "2nd propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  *sProp_L[2];
  sProp_L[0] = new Core::StochasticPropagator< 4 >(L, T);
  sProp_L[1] = new Core::StochasticPropagator< 4 >(L, T);

  tmpProp->reconstruct_doublet(*(sProp_L[0]), *(sProp_L[1]), gauge_field, kappa, mu_s);
  if (weave.isRoot())
    std::cout << "2nd propagator successfully reconsructed\n" << std::endl;

/* ******************************************************************************* */
/* ******** load and prepare propagators from timeslice t_R ********************** */
/* ******************************************************************************* */


  Tool::IO::load(tmpProp, propfilesD_R, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
      std::cout << "3rd propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  *dProp_R[2];
  dProp_R[0] = new Core::StochasticPropagator< 4 >(L, T);
  dProp_R[1] = new Core::StochasticPropagator< 4 >(L, T);
  
  tmpProp->reconstruct_doublet(*(dProp_R[0]),*(dProp_R[1]), gauge_field, kappa, mu_d);

  if (weave.isRoot())
    std::cout << "3rd propagator successfully reconstructed\n" << std::endl;

  Tool::IO::load(tmpProp, propfilesS_R, Tool::IO::fileSCIDAC,32);
  if (weave.isRoot())
    std::cout << "4th propagator successfully loaded\n" << std::endl;

  Core::StochasticPropagator< 4 >  *sProp_R[2];
  sProp_R[0] = new Core::StochasticPropagator< 4 >(L, T);
  sProp_R[1] = new Core::StochasticPropagator< 4 >(L, T);

  tmpProp->reconstruct_doublet(*(sProp_R[0]), *(sProp_R[1]), gauge_field, kappa, mu_s);

  if (weave.isRoot())
    std::cout << "4th propagator successfully reconstructed\n" << std::endl;

  delete tmpProp;

  if(phys_base)
  {
    dProp_R[0]->rotateToPhysicalBasis(0);
    dProp_R[1]->rotateToPhysicalBasis(1);
    dProp_L[0]->rotateToPhysicalBasis(0);
    dProp_L[1]->rotateToPhysicalBasis(1);
    sProp_R[0]->rotateToPhysicalBasis(0);
    sProp_R[1]->rotateToPhysicalBasis(1);
    sProp_L[0]->rotateToPhysicalBasis(0);
    sProp_L[1]->rotateToPhysicalBasis(1);


  }

  std::ostringstream out;
  out << "C3^{connected}    = sum_{vec{x}} Tr[s(x,tR) g5 revert(d'(x,tR)) Gamma1 s(x,t_L) g5 revert(d(x,t_L)) Gamma2] \n";
  out << "C3^{disconnected} = sum_{vec{x}} Tr[s(x,tR) g5 revert(d'(x,tR)) Gamma2] Tr[ s(x,t_L) g5 revert(d(x,t_L)) Gamma1] \n";
  out<<  "with Gamma1 = Gamma2 = Gamma \n";
  out<<  "where : \n";
  out<<  "        correlator    Gamma \n ";
  out<<  "           -1            I   \n ";
  out<<  "            1           g5   \n ";
  out<<  "            2           g1   \n ";
  out<<  "            3           g2   \n ";
  out<<  "            4           g3   \n ";
  out<<  "            5           g4   \n ";
  out<<  "            6          g1g5   \n ";
  out<<  "            7          g2g5   \n ";
  out<<  "            8          g3g5   \n ";
  out<<  "            9          g4g5   \n ";
  out<<  "           10          g1g2   \n ";
  out<<  "           11          g1g3   \n ";
  out<<  "           12          g1g4   \n ";
  out<<  "           13          g2g3   \n ";
  out<<  "           14          g2g4   \n ";
  out<<  "           15          g3g4   \n ";
  out<<  "Remember to reconstruct the Supersymmetric basis in the analysis program (see notes for the relation) \n";
  out<<  "and promediate over the different configurations with one r of opposite sign \n";

  Print(out.str());

  std::vector< std::pair< Base::Operator, Base::Operator > > *operator_combinations;
  std::vector< Core::Correlator< Dirac::Matrix > > *threepointsDis;
  std::vector< Core::Correlator< Dirac::Matrix > > *threepointsCon;

  //auxiliar index for r=0,1
  int idL, idR, isL, isR;
  char outputname_char[1024];
  strcpy(outputname_char, outputname.c_str());

  {
   // now we want to tell the contraction routine which operator (gamma) combinations to use
  // note that for Bk we always have gamma5 sources

    for (size_t i = 0; i < operators.size(); i += 1)
    {
      std::ostringstream out;
      out << "operator combination no ";
      out.width(3);
      out << i << ":  correlator ";
      out.width(3);
      out << std::showpos << operators[i] ;
      Print(out.str());
      operator_combinations= new std::vector< std::pair< Base::Operator, Base::Operator > >;
      (*operator_combinations).push_back(std::make_pair(Tool::convertIntToOperator(5),
                                                     Tool::convertIntToOperator(operators[i])));

      for (size_t j = 0; j < rcombinations.size(); j += 4){
		idL = rcombinations[j];
                isL = rcombinations[j+1]; 
 		idR = rcombinations[j+2];
  		isR = rcombinations[j+3];

		if (idL+isL+idR+isR==2 || idL+isL+idR+isR==0 || idL+isL+idR+isR==4) std::cout<<"WARNING: NOT A O(a) IMPROVED COMBINATION OF r"<<std::endl;
		else{
      			 if (connected==0 ){
      				threepointsDis = new std::vector< Core::Correlator< Dirac::Matrix > >(Contract::BK_threepoint_disconnected_stochastic(*(dProp_R[idR]), *(sProp_R[isR]), *(dProp_L[idL]), *(sProp_L[isL]), (*operator_combinations)));

				if(weave.isRoot())
   				 {
     					 std::ofstream fout(outputname_char);
	   				 std::ostringstream out;
					 out.flush();
					 out<< "RIGHT: (rd,rs)=("<<idR<<","<<isR<<")  LEFT: (rd,rs)=("<<idL<<","<<isL<<") \n";
					 out<<" Disconnected diagram : \n";
					 out << (*threepointsDis)[0];
      					 out.flush();
					 Print(out.str(), fout);
					fout.close();
	   			}
    				delete threepointsDis;
    			 }
    
 			 if (connected==1 || connected==2){
	    			  threepointsCon = new std::vector< Core::Correlator< Dirac::Matrix > >(Contract::BK_threepoint_connected_stochastic(*(dProp_R[idR]), *(sProp_R[isR]), *(dProp_L[idL]), *(sProp_L[isL]),(*operator_combinations)));
                                  if(weave.isRoot())
                                 {
                                         std::ofstream fout(outputname_char);
                                         std::ostringstream out;
                                         out.flush();
                                         out<< "RIGHT: (rd,rs)=("<<idR<<","<<isR<<")  LEFT: (rd,rs)=("<<idL<<","<<isL<<") \n";
                                         out<<" Connected diagram : \n";
                                         out << (*threepointsCon)[0];
                                         out.flush();
                                         Print(out.str(), fout);
                                         fout.close();
                                }


     	   			 delete threepointsCon;
			}
			if (connected==2){

				  threepointsDis = new std::vector< Core::Correlator< Dirac::Matrix > >(Contract::BK_threepoint_disconnected_stochastic(*(dProp_R[idR]), *(sProp_R[isR]), *(dProp_L[idL]), *(sProp_L[isL]), (*operator_combinations)));
				  threepointsCon = new std::vector< Core::Correlator< Dirac::Matrix > >(Contract::BK_threepoint_connected_stochastic(*(dProp_R[idR]), *(sProp_R[isR]), *(dProp_L[idL]), *(sProp_L[isL]),(*operator_combinations)));
				 if(weave.isRoot())
                                 {
                                         std::ofstream fout(outputname_char);
                                         std::ostringstream out;
                                         out.flush();
                                         out<< "RIGHT: (rd,rs)=("<<idR<<","<<isR<<")  LEFT: (rd,rs)=("<<idL<<","<<isL<<") \n";
					 out<<" Disconnected diagram : \n";
                                         out << (*threepointsDis)[0];
                                         out<<" Connected diagram : \n";
                                         out << (*threepointsCon)[0];
                                         out.flush();
                                         Print(out.str(), fout);
                                         fout.close();
                                }

  			 }
		}//end else
 	 } //end for r combination
   } // end for operator combination
 }



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

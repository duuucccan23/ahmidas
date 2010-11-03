//This program produce the sequential propagator 01 or 10 (charged
//ones) starting from the output of tmLQCD inverter over a source
//generated with "generate_meson_sequential_propagator"

//The propagator in input are in twisted basis so we have to rotate
//them to physical one, by multiplying each side for the appropriate
//(1+-ig5)/sqrt(2) factor. If the S0 original propagator is of the
//type "0", the second is the "1" beacuse we are studying charged case.

//In the parameter file you have to specify the following files:
//File0: gauge configuration, needed to produce the correct propagator
//File1: input propagator

//Output propagators file names are created automatically adding
//"final" to filenames

//Explanation of what is done. Let's take the S1 of of a '0' S0

// We want to compute D_up^-1*g5*D_dw^-1 , to obtain it we do:
// 0.25*(1+ig5)*(D_+)^-1*(1+ig5)*g5*(1-ig5)*(D_-)^-1*(1-ig5)
//
// this was the output of "generate_meson_sequential_source":
//                   0.5*(1+ig5)*g5*(1-ig5)*(D_-)^-1
//
// now we take as input the output of tmLQCD inverter over it:
//       0.5*(D_+D_-)^-1*(1+ig5)*g5*(1-ig5)*(D_-)^-1
//
// we get D_+^-1 by applying D_-:
//   0.5*D_-*(D_+D_-)^-1*(1+ig5)*g5*(1-ig5)*(D_-)^-1 =
//            0.5*D_+^-1*(1+ig5)*g5*(1-ig5)*(D_-)^-1
//
// and now we rotate left and right apprioprately to get what we want 


#include <complex>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>

int main(int narg,char **arg)
{
  size_t L,T;

  if(narg<2)
    {
      std::cerr<<"Use: "<<arg[0]<<" parameters_file"<<std::endl;
      exit(1);
    }

  //read inputs from parameter file
  std::map<std::string,double> floats;
  std::vector<std::vector<std::string> > files;
  std::vector<std::string> file_out;
  Input::FileReader reader(arg[1]);
  reader.initializeParameters(L,T,files,floats);
  double kappa=floats["kappa"];
  double mu=floats["mu"];
  double thetax=floats["thetax"];
  double thetay=floats["thetay"];
  double thetaz=floats["thetaz"];
  double thetat=floats["thetat"];
  bool s0_flav=(bool)floats["s0_flav"];
  bool debug=(bool)floats["debug"]; //0=no,!=0 =yes

  //the weave!
  Base::Weave weave(L,T);

  //Prepare output file names
  for(size_t j=0;j<files[1].size();j++)
    file_out.push_back(files[1][j]+".final");
  
  //write some output concerning what is going to be read
  if(weave.isRoot())
    {
      std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;
      std::cout<<"kappa="<<kappa<<", mu="<<mu<<std::endl;
      std::cout<<"thetax="<<thetax<<", ";
      std::cout<<"thetay="<<thetay<<", ";
      std::cout<<"thetaz="<<thetaz<<", ";
      std::cout<<"thetat="<<thetat<<std::endl;
      std::cout<<"flavor of the s0 line: "<<s0_flav<<std::endl;
    }

  if(weave.isRoot() and debug==true)
    {
      std::cout<<"The following files are going to be read:"<<std::endl;
      
      for(size_t i=0;i<files.size();i++)
	for(size_t j=0;j<files[i].size();j++)
	  std::cout<<(files[i])[j]<<std::endl;

      std::cout<<std::endl<<"Output files "<<std::endl;
      for(size_t j=0;j<files[1].size();j++)
	std::cout<<file_out[j]<<std::endl;
    }

  //read the gauge configuration, needed to cancel the D+ or D- from tmLQCD output
  Core::Field<QCD::Gauge> gauge_field(L,T);
  Tool::IO::load(&gauge_field,files[0][0],Tool::IO::fileILDG);
  if(weave.isRoot()) std::cout<<std::endl<<"gauge field successfully loaded"<<std::endl<<std::endl;

  //read the solution of (D+D-)^-1
  PropagatorType DD_prop(L,T);
  Tool::IO::load(&DD_prop,files[1],Tool::IO::fileSCIDAC);
  if(weave.isRoot()) std::cout<<"(D+D-)^-1 prop successfully loaded"<<std::endl<<std::endl;

  PropagatorType prop_out(L,T);
  Dirac::Gamma<5> gamma5;

  //Now apply the dirac operator. It has to be the opposite of what
  //you want in the s1 line, so the same of s0_flav (0->-,1->+)
  if(s0_flav==0) prop_out=DD_prop.applyDiracOperator(gauge_field,kappa,-mu,thetat,thetax,thetay,thetaz); 
  else           prop_out=DD_prop.applyDiracOperator(gauge_field,kappa,+mu,thetat,thetax,thetay,thetaz); 

  prop_out.rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
  
  //Now rotate on the source as s0 (0->-,1->+)
  PropagatorType temp(prop_out);
  temp.isolate();
  if(s0_flav==0) temp*=std::complex<double>(0,-1);
  else           temp*=std::complex<double>(0,+1);
  temp*=(gamma5);
  prop_out+=temp;

  //Now rotate on the sink as s1 (0->+,1->-)
  temp=prop_out;
  temp.isolate();
  if(s0_flav==0) temp*=std::complex<double>(0,+1);
  else           temp*=std::complex<double>(0,-1);
  temp.rightMultiply(gamma5);
  prop_out+=temp;

  //Multiply by 1/sqrt(2)**2
  prop_out*=0.5;
  
  Tool::IO::save(&prop_out,file_out,Tool::IO::fileSCIDAC);
  if(weave.isRoot()) std::cout<<"sequential propagator saved"<<std::endl;

  if(weave.isRoot()) std::cout<<std::endl<<"everything ok so exiting!"<<std::endl<<std::endl;

  return EXIT_SUCCESS;
}

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L1/Tool/IO.h>
#include <L2/Contract/Meson.h>
#include <L2/Input/FileReader.h>

using namespace std;

int main(int narg, char **arg)
{
  if(narg<2)
    {
      cerr<<"Use "<<arg[0]<<" inputfile"<<endl;
      cerr<<"with:"<<endl;
      cerr<<" L,T"<<endl;
      cerr<<" Kappa"<<endl;
      cerr<<" Theta x,y,z,t"<<endl;
      cerr<<" Tslice of source"<<endl;
      cerr<<" Norm of the source (for normalization)"<<endl;
      cerr<<" Gauge configuration file"<<endl;
      cerr<<" Number of contractions"<<endl;
      cerr<<" List of contraction pairs"<<endl;
      cerr<<" Number of masses"<<endl;
      cerr<<" List of masses, base name of the propagator, suffix of the propagator"<<endl;
      exit(1);
    }

  int L,T;
  double kappa;
  double thex,they,thez,thet;
  double norm;
  int tss;
  string gauge_file;
  int nmass,ncontr;
  std::vector<std::pair<Base::Operator,Base::Operator> > op_comb;
  ifstream ifile(arg[1]);

  ifile>>L>>T;
  Base::Weave weave(L,T);
  bool isr=weave.isRoot();

  ifile>>kappa;
  if(isr) cout<<"kappa: "<<kappa<<endl;
  ifile>>thex>>they>>thez>>thet;
  if(isr) cout<<"theta: "<<thex<<" "<<they<<" "<<thez<<" "<<thet<<endl;
  ifile>>tss;
  if(isr) cout<<"tss: "<<tss<<endl;
  ifile>>norm;
  if(isr) cout<<"norm: "<<norm<<endl;

  //Load gauge field
  Core::Field<QCD::Gauge> gauge_field(L,T);
  ifile>>gauge_file;
  Tool::IO::load(&gauge_field,gauge_file,Tool::IO::fileILDG);
  if(isr) cout<<"Gauge "<<gauge_file<<" loaded"<<endl;

  ifile>>ncontr;
  if(isr) cout<<"Ncontr: "<<ncontr<<endl;
  int op1,op2;
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      ifile>>op1>>op2;
      op_comb.push_back(std::make_pair(Tool::convertIntToOperator(op1),Tool::convertIntToOperator(op2)));
      if(isr) cout<<"Combination "<<op1<<" "<<op2<<" added"<<endl;
    }

  ifile>>nmass;
  if(isr) cout<<"Nmasses "<<nmass<<endl;
  PropagatorType ***prop=new PropagatorType**[nmass];
  string path_beg,path_end;
  char ind[3];
  vector<string> prop_path;
  double *mass=new double[nmass];
  if(isr) cout<<"We wait for "<<FilesPerProp<<" files per prop"<<endl;
  for(int imass=0;imass<nmass;imass++)
    {
      prop[imass]=new PropagatorType*[2];
      prop[imass][0]=new PropagatorType(L,T);
      prop[imass][1]=new PropagatorType(L,T);

      PropagatorType temp_prop(L,T);

      ifile>>mass[imass];
      ifile>>path_beg;
      //ifile>>path_end;
      path_end="";
      if(isr) cout<<"Mass "<<mass[imass]<<" "<<path_beg<<" added"<<endl;

      prop_path.clear();

      for(int ifi=0;ifi<FilesPerProp;ifi++)
	{
	  sprintf(ind,"%02d",ifi);
	  prop_path.push_back(path_beg+ind+path_end);
	}
      
      Tool::IO::load(&temp_prop,prop_path,Tool::IO::fileSCIDAC);

      temp_prop.reconstruct_doublet(*(prop[imass][0]),*(prop[imass][1]),gauge_field,kappa,mass[imass],thet,thex,they,thez);
      
      prop[imass][0]->rotateToPhysicalBasis(0);
      prop[imass][1]->rotateToPhysicalBasis(1);
    }

  ofstream fout("correlators.dat");
  
  std::vector<Core::Correlator<Dirac::Matrix> > C2;

  for(int im1=0;im1<nmass;im1++)
    for(int im2=im1;im2<nmass;im2++)
      for(int if1=0;if1<2;if1++)
	for(int if2=0;if2<2;if2++)
	  {
#ifdef UltraStocCase
	    C2=Contract::light_meson_twopoint_ultrastochastic(*(prop[im1][if1]),*(prop[im2][if2]),op_comb);
#elif defined StocCase
	    C2=Contract::light_meson_twopoint_stochastic(*(prop[im1][if1]),*(prop[im2][if2]),op_comb);
#else
	    C2=Contract::light_meson_twopoint(*(prop[im1][if1]),*(prop[im2][if2]),op_comb);
#endif
	    for(int icontr=0;icontr<ncontr;icontr++)
	      {
		C2[icontr].setOffset(tss);
		C2[icontr]*=1/norm;
		if(weave.isRoot())
		  {
		    fout<<mass[im1]<<" "<<mass[im2]<<" "<<if1<<" "<<if2<<" "<<op_comb[icontr].first<<endl;
		    for(int t=0;t<T;t++)
		      {
			complex<double> corr=C2[icontr][t].trace();
			fout<<corr.real()<<" "<<corr.imag()<<endl;
		      }
		  }
	      }
	  }
  
  return EXIT_SUCCESS;
}

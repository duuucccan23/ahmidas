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

int safe_mod(int x,int y)
{
  if(x>=0)
    return x%y;
  else
    return (y-(abs(x)%y))%y;
}

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
  Input::FileReader reader(arg[1]);
  reader.initializeParameters(L,T,files,floats);

  //the weave!
  Base::Weave weave(L,T);

  //write some output concerning what is going to be read
  if(weave.isRoot())
    {
      std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;
      std::cout<<files[0].size()<<std::endl;
    }

  //read the gauge configuration, needed to create u and d propagators
  Core::Field<QCD::Gauge> gauge_field(L,T);
  Tool::IO::load(&gauge_field,files[0][0].c_str(),Tool::IO::fileILDG);
  if(weave.isRoot()) std::cout<<std::endl<<"gauge field successfully loaded"<<std::endl<<std::endl;

  double plaqt=Tool::temporalPlaquette(gauge_field);
  double plaqs=Tool::spatialPlaquette(gauge_field);
  double plaq=(plaqs+plaqt)*3/2;
  std::cout<<plaqs<<" "<<plaqt<<" "<<plaq<<std::endl;

  {
    std::ofstream fdel("lattice2");
    fdel.precision(12);
    fdel<<1<<std::endl;
    for(int ic1=0;ic1<3;ic1++)
      for(int ic2=0;ic2<3;ic2++)
	for(int idir=0;idir<4;idir++)
	  {
	    for(int p=0;p<2;p++)
	      for(size_t t=+0;t<T;t++)
		for(size_t z=+0;z<L;z++)
		  for(size_t y=0;y<L;y++)
		    for(size_t x=0;x<L;x++)
		      {
			int iw=weave.globalCoordToLocalIndex(x%L,y%L,z%L,t%T);
			int sum=x+y+z+t;
			if(sum-2*(sum/2)==p)
			  fdel<<gauge_field[iw][idir](ic1,ic2).real()<<" "<<gauge_field[iw][idir](ic1,ic2).imag()<<std::endl;
		      }
	  }
  }
	    
  

  FILE *fout=NULL;

  fout=fopen((files[0][0]+"conv").c_str(),"wb");

  fwrite(&T,sizeof(int),1,fout);
  fwrite(&L,sizeof(int),1,fout);
  fwrite(&L,sizeof(int),1,fout);
  fwrite(&L,sizeof(int),1,fout);

  std::cout<<plaqs<<" "<<plaqt<<" "<<plaq<<std::endl;
  fwrite(&plaq,sizeof(double),1,fout);

  double temp[2];

  for(size_t t=0;t<T;t++)
    for(size_t z=0;z<L;z++)
      for(size_t y=0;y<L;y++)
	for(size_t x=0;x<L;x++)
	  {
	    int p=(x+y+z+t)%2;
	    if(p)
	      {
		int pos[4];
		int sizes[4]={T,L,L,L};
		
		pos[0]=t%T;
		pos[1]=z%L;
		pos[2]=y%L;
		pos[3]=x%L;

		int iw=weave.globalCoordToLocalIndex(x%L,y%L,z%L,t%T);
		
		for(int mu=0;mu<4;mu++)
		  {
		    int nei[4];
		    for(int nu=0;nu<4;++nu) nei[nu]=pos[nu];
		    nei[mu]=safe_mod(pos[mu]-1,sizes[mu]);
		    int iz=weave.globalCoordToLocalIndex(nei[3],nei[2],nei[1],nei[0]);

		    for(int ic1=0;ic1<3;ic1++)
		      for(int ic2=0;ic2<3;ic2++)
			{
			  temp[0]=gauge_field[iw][3-mu](ic1,ic2).real();
			  temp[1]=gauge_field[iw][3-mu](ic1,ic2).imag();
			  fwrite(&(temp),sizeof(double),2,fout);
			}
		    for(int ic1=0;ic1<3;ic1++)
		      for(int ic2=0;ic2<3;ic2++)
			{
			  temp[0]=gauge_field[iz][3-mu](ic1,ic2).real();
			  temp[1]=gauge_field[iz][3-mu](ic1,ic2).imag();
			  fwrite(&(temp),sizeof(double),2,fout);
			}
		    /*
		    //std::complex<double> *punt=(std::complex<double>*)&(gauge_field[iw][mu]);
		    //for(int ic1=0;ic1<9;ic1++) std::cout<<punt[ic1]<<std::endl;
		    
		    //for(int ic1=0;ic1<3;ic1++)
		    // for(int ic2=0;ic2<3;ic2++)
		    //std::cout<<gauge_field[iw][mu](ic1,ic2)<<std::endl;

		    fwrite(&(gauge_field[iw][3-mu]),sizeof(SU3::Matrix),1,fout);
		    fwrite(&(gauge_field[iz][3-mu]),sizeof(SU3::Matrix),1,fout);
		    */
		  }
	      }
	  }

  if(weave.isRoot()) std::cout<<std::endl<<"everything ok so exiting!"<<std::endl<<std::endl;

  return EXIT_SUCCESS;
}

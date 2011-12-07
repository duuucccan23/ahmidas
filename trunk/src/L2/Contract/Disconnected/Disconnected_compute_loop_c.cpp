#include "Disconnected.ih"
#include <L0/Base/Weave.h>
#include <L0/QCD/Tensor.h>

class complex12 
{
	private :
    std::complex<double> data[12];

	public :

	
    inline complex12() { std::fill_n(data,12,0.0); }
    inline complex12(std::complex<double> c) { std::fill_n(data,12,c); }
// inline ~complex12() { delete data ;}
	inline std::complex<double> &operator[](size_t index){ return data[index]; } 
	inline complex12 &operator+=(complex12 const &other){ std::transform (this->data, this->data + 12, other.data, this->data, std::plus< std::complex< double > >()); return *this; }
};

class complex192 
{
	private :
		std::complex<double> data[192];

	public :


		inline complex192() { std::fill_n(data,192,0.0); }
		inline complex192(std::complex<double> c) { std::fill_n(data,192,c); }
//		inline ~complex192() { delete data ;}
		inline std::complex<double> &operator[](size_t index){ return data[index]; } 
		inline complex192 &operator+=(complex192 const &other){ std::transform (this->data, this->data + 192, other.data, this->data, std::plus< std::complex< double > >()); return *this; }
		inline complex192 &operator=(complex192 const &other) {
			if (&other != this)
				std::copy(other.data, other.data + 192, data); 
			return *this;
		}
		inline complex192 operator*(std::complex<double> const &factor) const {
			complex192 tmp;

			std::transform (this->data, this->data + 192, tmp.data, std::bind2nd(std::multiplies< std::complex< double > >(), factor)); 
			return tmp;
		}

};


namespace Contract
{

	// compute xi^* x Gamma x psi
	std::vector< std::complex<double>  > compute_loop(
			Core::Propagator const &xi, Core::Propagator const &psi,
			std::vector< Base::HermitianBilinearOperator > ops)
	{

		assert(xi.T() == psi.T() && psi.L() == psi.L());
		Base::Weave weave(xi.L(), xi.T());
		//std::vector< Core::Correlator< Dirac::Matrix > > twopoints;
		std::vector<  std::complex<double>  > twopoints;
		//Core::Propagator xi_conj(xi);

		//conjugation
		//xi_conj.conjugate();

		// loop over all operator combinations

		for (size_t iOp=0; iOp<ops.size(); iOp++)
		{
			if (weave.isRoot()) std::cout << "Nop = " << iOp << std::endl;

			if (weave.isRoot()) std::cout << "Not safe function ... exiting ! " << std::endl;
			exit(0);
			Core::Propagator Gamma_psi(psi);

			//apply operator X
			Gamma_psi.rightMultiplyOperator(ops[iOp]);

			Core::Field< complex12 > res(xi.L(),xi.T());
			Core::Field< complex12 >::iterator K(res.begin());

			//			Core::Field< std::complex<double> > field(xi.L(),xi.T());
			//			std::vector< Core::Field< std::complex<double> > > res(12,field);


			//			std::vector< Core::Field< std::complex<double> > >::iterator K(res.begin());
			Core::Propagator::const_iterator J(Gamma_psi.begin());
			Core::Propagator::const_iterator I(xi.begin());

			while(I != xi.end())
			{
				for (size_t i=0; i<12;i++)
				{
					(*K)[i]=innerProduct((*I)[i],(*J)[i]);
				}
				++K;
				++J;
				++I;
			}
			Core::Correlator< complex12 > twopoint(res);

			twopoint.sumOverSpatialVolume(); 
			twopoint.deleteField();

			//		std::cout << "C" <<twopoint[0][0] <<std::endl;

			if (weave.isRoot()) 
			{
				/*	for (size_t i=0; i<12;i++)
					{
					for (size_t t=0;t < xi.T();t++)
					{
					std::cout << "Op = "<< iOp << " " << i <<" " << t  <<" "<<twopoint[t][i].real() << " "<< twopoint[t][i].imag()  <<std::endl;
					}
					}*/

				for (size_t t=0;t < xi.T();t++)
				{
					for (size_t i=0; i<12;i++)
					{
						twopoints.push_back(twopoint[t][i]);
					}
				}
			}

		}
		//contract
		//Core::Correlator< Dirac::Matrix >twopoint( xi_conj * Gamma_psi);
		//sum over space



		//		twopoint.deleteField();	
		//accumulate 
		//		 twopoints.push_back(twopoint);




		return twopoints;
	}	

	std::vector< std::complex<double>  > compute_loop_new(
			Core::Propagator const &xi, Core::Propagator const &psi,
			std::vector< Base::HermitianBilinearOperator > ops)
	{

		clock_t start, finish;
		start = clock();

		assert(xi.T() == psi.T() && psi.L() == psi.L());

		Base::Weave weave(xi.L(), xi.T());
		std::vector<  std::complex<double>  > twopoints;

		//Core::Propagator  psi_conj(psi);
		//Core::Propagator   xi_conj(xi);

		//conjugation
		//psi_conj.conjugate();
		//xi_conj.conjugate();

		Core::Propagator Gamma_psi0(psi);
		//declare iterator of 192 complex. To 
		Core::Field< complex192 > res(xi.L(),xi.T());
		Core::Field< complex192 >::iterator K(res.begin());

		Core::Propagator::const_iterator J0(Gamma_psi0.begin());

		Core::Propagator::const_iterator I(xi.begin());

		/*Base::Weave *d_weave;

		  static Core::Field< size_t > *s_timelabel;

		  s_timelabel =  new Core::Field< size_t >(psi.L() , psi.T());

		  size_t localIndex;
		  for(size_t idx_T = 0; idx_T < psi.T(); idx_T++)
		  {
		  for(size_t idx_Z = 0; idx_Z < psi.L(); idx_Z++)
		  {
		  for(size_t idx_Y = 0; idx_Y < psi.L(); idx_Y++)
		  {
		  for(size_t idx_X = 0; idx_X < psi.L(); idx_X++)
		  {
		  localIndex = d_weave->globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
		  if (localIndex == d_weave->localVolume())
		  continue;
		  (*s_timelabel)[localIndex] = idx_T;
		  }
		  }
		  }
		  }
		 */
//		if (weave.isRoot()) 
//		{
//			std::cout << "before loop iterator ..."<<std::endl;
//		}


		clock_t starti, finishi;
		starti = clock();


		QCD::Tensor tmp0;
		QCD::Tensor tmp1;
		QCD::Tensor tmp2;
		QCD::Tensor tmp3;
		QCD::Tensor tmp4;
		QCD::Tensor tmp5;
		QCD::Tensor tmp6;
		QCD::Tensor tmp7;
		QCD::Tensor tmp8;
		QCD::Tensor tmp9;
		QCD::Tensor tmp10;
		QCD::Tensor tmp11;
		QCD::Tensor tmp12;
		QCD::Tensor tmp13;
		QCD::Tensor tmp14;
		QCD::Tensor tmp15;


		while(I != xi.end())
		{

			tmp0=*J0;
			tmp1 = tmp0;
			tmp2 = tmp0;
			tmp3 = tmp0;
			tmp4 = tmp0;
			tmp5 = tmp0;
			tmp6 = tmp0;
			tmp7 = tmp0;
			tmp8 = tmp0;
			tmp9 = tmp0;
			tmp10= tmp0;
			tmp11= tmp0;
			tmp12= tmp0;
			tmp13= tmp0;
			tmp14= tmp0;
			tmp15= tmp0;

			
			//QCD::Tensor tmp0(*J0);
/*			QCD::Tensor tmp1(*J0);
			QCD::Tensor tmp2(*J0);
			QCD::Tensor tmp3(*J0);
			QCD::Tensor tmp4(*J0);
			QCD::Tensor tmp5(*J0);
			QCD::Tensor tmp6(*J0);
			QCD::Tensor tmp7(*J0);
			QCD::Tensor tmp8(*J0);
			QCD::Tensor tmp9(*J0);
			QCD::Tensor tmp10(*J0);
			QCD::Tensor tmp11(*J0);
			QCD::Tensor tmp12(*J0);
			QCD::Tensor tmp13(*J0);
			QCD::Tensor tmp14(*J0);
			QCD::Tensor tmp15(*J0);*/

			tmp0.rightMultiplyOperator(ops[0]);
			tmp1.rightMultiplyOperator(ops[1]);
			tmp2.rightMultiplyOperator(ops[2]);
			tmp3.rightMultiplyOperator(ops[3]);
			tmp4.rightMultiplyOperator(ops[4]);
			tmp5.rightMultiplyOperator(ops[5]);
			tmp6.rightMultiplyOperator(ops[6]);
			tmp7.rightMultiplyOperator(ops[7]);
			tmp8.rightMultiplyOperator(ops[8]);
			tmp9.rightMultiplyOperator(ops[9]);
			tmp10.rightMultiplyOperator(ops[10]);
			tmp11.rightMultiplyOperator(ops[11]);
			tmp12.rightMultiplyOperator(ops[12]);
			tmp13.rightMultiplyOperator(ops[13]);
			tmp14.rightMultiplyOperator(ops[14]);
			tmp15.rightMultiplyOperator(ops[15]);

			for (size_t i=0; i<12;i++)
			{

				(*K)[0+16*i]=innerProduct((*I)[i],tmp0[i]);
				(*K)[1+16*i]=innerProduct((*I)[i],tmp1[i]);
				(*K)[2+16*i]=innerProduct((*I)[i],tmp2[i]);
				(*K)[3+16*i]=innerProduct((*I)[i],tmp3[i]);
				(*K)[4+16*i]=innerProduct((*I)[i],tmp4[i]);
				(*K)[5+16*i]=innerProduct((*I)[i],tmp5[i]);
				(*K)[6+16*i]=innerProduct((*I)[i],tmp6[i]);
				(*K)[7+16*i]=innerProduct((*I)[i],tmp7[i]);
				(*K)[8+16*i]=innerProduct((*I)[i],tmp8[i]);
				(*K)[9+16*i]=innerProduct((*I)[i],tmp9[i]);
				(*K)[10+16*i]=innerProduct((*I)[i],tmp10[i]);
				(*K)[11+16*i]=innerProduct((*I)[i],tmp11[i]);
				(*K)[12+16*i]=innerProduct((*I)[i],tmp12[i]);
				(*K)[13+16*i]=innerProduct((*I)[i],tmp13[i]);
				(*K)[14+16*i]=innerProduct((*I)[i],tmp14[i]);
				(*K)[15+16*i]=innerProduct((*I)[i],tmp15[i]);
			}
			++K;
			++J0;
			++I;
		}

		finishi = clock();
//		if (weave.isRoot())
//			std::cout << "loop iterator  "<< double(finishi - starti)/CLOCKS_PER_SEC  << "seconds." << std::endl;


		clock_t startv, finishv;
		startv = clock();


		Core::Correlator< complex192 > twopoint(res);

		//sum over space
		twopoint.sumOverSpatialVolume(); 
//		twopoint.deleteField();

		finishv = clock();
//		if (weave.isRoot())
//			std::cout << "loop vector  "<< double(finishv - startv)/CLOCKS_PER_SEC  << "seconds." << std::endl;


		if (weave.isRoot()) 
		{

			//accumulate 
			for (size_t i=0; i < 192;i++)
			{
				for (size_t t=0;t < xi.T();t++)
				{
					twopoints.push_back(twopoint[t][i]);
				}
			}

		}



		finish = clock();
//		if (weave.isRoot())
//			std::cout << "Computation disconnected_loops_test in  "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;

		return twopoints;


	}

}

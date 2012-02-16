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
	std::vector< std::complex<double>  > compute_loop_new(
			Core::Propagator const &xi, Core::Propagator const &psi,
			std::vector< Base::HermitianBilinearOperator > ops,
			int const * const position_offset, std::vector< int* > const &momenta,int const tsrc)
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

		QCD::Tensor tmp0;
		while(I != xi.end())
		{

			tmp0=*J0;
			//QCD::Tensor tmp0(*J0);
			QCD::Tensor tmp1(*J0);
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
			QCD::Tensor tmp15(*J0);

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



		Core::Correlator< complex192 > twopoint(res);




		twopoint.prepareMomentumProjection(position_offset);

		std::vector< Core::Correlator< complex192 > > tmp(twopoint.momentumProjection(momenta));

		twopoint.deleteField();

		if (weave.isRoot()) 
		{

			for (size_t m=0; m <momenta.size();m++)
			{
				tmp[m].setOffset(tsrc);
				for (size_t i=0; i < 192;i++)
				{
					for (size_t t=0;t < xi.T();t++)
					{
						twopoints.push_back(tmp[m][t][i]);
					}
				}
			}
		}



		return twopoints;
	}

}

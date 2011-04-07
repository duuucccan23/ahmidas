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



		assert(xi.T() == psi.T() && psi.L() == psi.L());

		Base::Weave weave(xi.L(), xi.T());
		std::vector<  std::complex<double>  > twopoints;

		//Core::Propagator  psi_conj(psi);
		//Core::Propagator   xi_conj(xi);

		//conjugation
		//psi_conj.conjugate();
		//xi_conj.conjugate();

		// declaration of one propagator for each gamma structure. Memory consuming but efficient ...
		Core::Propagator Gamma_psi0(psi);
		Core::Propagator Gamma_psi1(psi);
		Core::Propagator Gamma_psi2(psi);
		Core::Propagator Gamma_psi3(psi);
		Core::Propagator Gamma_psi4(psi);
		Core::Propagator Gamma_psi5(psi);
		Core::Propagator Gamma_psi6(psi);
		Core::Propagator Gamma_psi7(psi);
		Core::Propagator Gamma_psi8(psi);
		Core::Propagator Gamma_psi9(psi);
		Core::Propagator Gamma_psi10(psi);
		Core::Propagator Gamma_psi11(psi);
		Core::Propagator Gamma_psi12(psi);
		Core::Propagator Gamma_psi13(psi);
		Core::Propagator Gamma_psi14(psi);
		Core::Propagator Gamma_psi15(psi);

		// compute Gamma psi_star 
		Gamma_psi0.rightMultiplyOperator(ops[0]);
		Gamma_psi1.rightMultiplyOperator(ops[1]);
		Gamma_psi2.rightMultiplyOperator(ops[2]);
		Gamma_psi3.rightMultiplyOperator(ops[3]);
		Gamma_psi4.rightMultiplyOperator(ops[4]);
		Gamma_psi5.rightMultiplyOperator(ops[5]);
		Gamma_psi6.rightMultiplyOperator(ops[6]);
		Gamma_psi7.rightMultiplyOperator(ops[7]);
		Gamma_psi8.rightMultiplyOperator(ops[8]);
		Gamma_psi9.rightMultiplyOperator(ops[9]);
		Gamma_psi10.rightMultiplyOperator(ops[10]);
		Gamma_psi11.rightMultiplyOperator(ops[11]);
		Gamma_psi12.rightMultiplyOperator(ops[12]);
		Gamma_psi13.rightMultiplyOperator(ops[13]);
		Gamma_psi14.rightMultiplyOperator(ops[14]);
		Gamma_psi15.rightMultiplyOperator(ops[15]);


		//declare iterator of 192 complex. To 
		Core::Field< complex192 > res(xi.L(),xi.T());
		Core::Field< complex192 >::iterator K(res.begin());

		Core::Propagator::const_iterator J0(Gamma_psi0.begin());
		Core::Propagator::const_iterator J1(Gamma_psi1.begin());
		Core::Propagator::const_iterator J2(Gamma_psi2.begin());
		Core::Propagator::const_iterator J3(Gamma_psi3.begin());
		Core::Propagator::const_iterator J4(Gamma_psi4.begin());
		Core::Propagator::const_iterator J5(Gamma_psi5.begin());
		Core::Propagator::const_iterator J6(Gamma_psi6.begin());
		Core::Propagator::const_iterator J7(Gamma_psi7.begin());
		Core::Propagator::const_iterator J8(Gamma_psi8.begin());
		Core::Propagator::const_iterator J9(Gamma_psi9.begin());
		Core::Propagator::const_iterator J10(Gamma_psi10.begin());
		Core::Propagator::const_iterator J11(Gamma_psi11.begin());
		Core::Propagator::const_iterator J12(Gamma_psi12.begin());
		Core::Propagator::const_iterator J13(Gamma_psi13.begin());
		Core::Propagator::const_iterator J14(Gamma_psi14.begin());
		Core::Propagator::const_iterator J15(Gamma_psi15.begin());

		Core::Propagator::const_iterator I(xi.begin());


		while(I != xi.end())
		{
			for (size_t i=0; i<12;i++)
			{

				(*K)[0+16*i]=innerProduct((*I)[i],(*J0)[i]);
				(*K)[1+16*i]=innerProduct((*I)[i],(*J1)[i]);
				(*K)[2+16*i]=innerProduct((*I)[i],(*J2)[i]);
				(*K)[3+16*i]=innerProduct((*I)[i],(*J3)[i]);
				(*K)[4+16*i]=innerProduct((*I)[i],(*J4)[i]);
				(*K)[5+16*i]=innerProduct((*I)[i],(*J5)[i]);
				(*K)[6+16*i]=innerProduct((*I)[i],(*J6)[i]);
				(*K)[7+16*i]=innerProduct((*I)[i],(*J7)[i]);
				(*K)[8+16*i]=innerProduct((*I)[i],(*J8)[i]);
				(*K)[9+16*i]=innerProduct((*I)[i],(*J9)[i]);
				(*K)[10+16*i]=innerProduct((*I)[i],(*J10)[i]);
				(*K)[11+16*i]=innerProduct((*I)[i],(*J11)[i]);
				(*K)[12+16*i]=innerProduct((*I)[i],(*J12)[i]);
				(*K)[13+16*i]=innerProduct((*I)[i],(*J13)[i]);
				(*K)[14+16*i]=innerProduct((*I)[i],(*J14)[i]);
				(*K)[15+16*i]=innerProduct((*I)[i],(*J15)[i]);
			}
			++K;
			++J0;
			++J1;
			++J2;
			++J3;
			++J4;
			++J5;
			++J6;
			++J7;
			++J8;
			++J9;
			++J10;
			++J11;
			++J12;
			++J13;
			++J14;
			++J15;
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

#include "Meson.ih"
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

class complex16
{
	private :
		std::complex<double> data[16];

	public :


		inline complex16() { std::fill_n(data,16,0.0); }
		inline complex16(std::complex<double> c) { std::fill_n(data,16,c); }
		inline std::complex<double> &operator[](size_t index){ return data[index]; } 
		inline complex16 &operator+=(complex16 const &other){ std::transform (this->data, this->data + 16, other.data, this->data, std::plus< std::complex< double > >()); return *this; }
		inline complex16 &operator=(complex16 const &other) {
			if (&other != this)
				std::copy(other.data, other.data + 16, data); 
			return *this;
		}
		inline complex16 operator*(std::complex<double> const &factor) const {
			complex16 tmp;

			std::transform (this->data, this->data + 16, tmp.data, std::bind2nd(std::multiplies< std::complex< double > >(), factor)); 
			return tmp;
		}

};



namespace Contract
{

	// compute phi_vv Gamma g5 phi^* x g5 Gamma
	std::vector< std::complex<double>  > compute_contract_meson(
			Core::Propagator const &phi, Core::Propagator const &phi_vv,
			std::vector< Base::HermitianBilinearOperator > ops)
	{



		assert(phi.T() == phi_vv.T() && phi.L() == phi_vv.L());

		Base::Weave weave(phi.L(), phi.T());
		std::vector<  std::complex<double>  > twopoints;

		Core::Propagator  g5_phi_conj(phi);
		g5_phi_conj.revert();

 if (weave.isRoot()) std::cout << "\n debut" << std::endl;

		// declaration of one propagator for each gamma structure. Memory consuming but efficient ...
		Core::Propagator Gamma_psi0(g5_phi_conj);
		Core::Propagator Gamma_psi1(g5_phi_conj);
		Core::Propagator Gamma_psi2(g5_phi_conj);
		Core::Propagator Gamma_psi3(g5_phi_conj);
		Core::Propagator Gamma_psi4(g5_phi_conj);
		Core::Propagator Gamma_psi5(g5_phi_conj);
		Core::Propagator Gamma_psi6(g5_phi_conj);
		Core::Propagator Gamma_psi7(g5_phi_conj);
		Core::Propagator Gamma_psi8(g5_phi_conj);
		Core::Propagator Gamma_psi9(g5_phi_conj);
		Core::Propagator Gamma_psi10(g5_phi_conj);
		Core::Propagator Gamma_psi11(g5_phi_conj);
		Core::Propagator Gamma_psi12(g5_phi_conj);
		Core::Propagator Gamma_psi13(g5_phi_conj);
		Core::Propagator Gamma_psi14(g5_phi_conj);
		Core::Propagator Gamma_psi15(g5_phi_conj);

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


		Gamma_psi0.leftMultiplyOperator(ops[0]);
		Gamma_psi1.leftMultiplyOperator(ops[1]);
		Gamma_psi2.leftMultiplyOperator(ops[2]);
		Gamma_psi3.leftMultiplyOperator(ops[3]);
		Gamma_psi4.leftMultiplyOperator(ops[4]);
		Gamma_psi5.leftMultiplyOperator(ops[5]);
		Gamma_psi6.leftMultiplyOperator(ops[6]);
		Gamma_psi7.leftMultiplyOperator(ops[7]);
		Gamma_psi8.leftMultiplyOperator(ops[8]);
		Gamma_psi9.leftMultiplyOperator(ops[9]);
		Gamma_psi10.leftMultiplyOperator(ops[10]);
		Gamma_psi11.leftMultiplyOperator(ops[11]);
		Gamma_psi12.leftMultiplyOperator(ops[12]);
		Gamma_psi13.leftMultiplyOperator(ops[13]);
		Gamma_psi14.leftMultiplyOperator(ops[14]);
		Gamma_psi15.leftMultiplyOperator(ops[15]);

		// compute Gamma psi_star 
		// compute Gamma psi_star 
		//declare iterator of 16 complex. To 
		Core::Field< complex16 > res(phi.L(),phi.T());
		Core::Field< complex16 >::iterator K(res.begin());

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

		Core::Propagator::const_iterator I(phi_vv.begin());

		if (weave.isRoot()) std::cout << "\n bla avant loop space.\n" << std::endl;


		while(I != phi_vv.end())
		{
			for (size_t i=0; i<12;i++)
			{
				if (i==0)
				{
				(*K)[0]=innerProduct((*I)[i],(*J0)[i]);
				(*K)[1]=innerProduct((*I)[i],(*J1)[i]);
				(*K)[2]=innerProduct((*I)[i],(*J2)[i]);
				(*K)[3]=innerProduct((*I)[i],(*J3)[i]);
				(*K)[4]=innerProduct((*I)[i],(*J4)[i]);
				(*K)[5]=innerProduct((*I)[i],(*J5)[i]);
				(*K)[6]=innerProduct((*I)[i],(*J6)[i]);
				(*K)[7]=innerProduct((*I)[i],(*J7)[i]);
				(*K)[8]=innerProduct((*I)[i],(*J8)[i]);
				(*K)[9]=innerProduct((*I)[i],(*J9)[i]);
				(*K)[10]=innerProduct((*I)[i],(*J10)[i]);
				(*K)[11]=innerProduct((*I)[i],(*J11)[i]);
				(*K)[12]=innerProduct((*I)[i],(*J12)[i]);
				(*K)[13]=innerProduct((*I)[i],(*J13)[i]);
				(*K)[14]=innerProduct((*I)[i],(*J14)[i]);
				(*K)[15]=innerProduct((*I)[i],(*J15)[i]);
				}
				else
				{
				(*K)[0]+=innerProduct((*I)[i],(*J0)[i]);
				(*K)[1]+=innerProduct((*I)[i],(*J1)[i]);
				(*K)[2]+=innerProduct((*I)[i],(*J2)[i]);
				(*K)[3]+=innerProduct((*I)[i],(*J3)[i]);
				(*K)[4]+=innerProduct((*I)[i],(*J4)[i]);
				(*K)[5]+=innerProduct((*I)[i],(*J5)[i]);
				(*K)[6]+=innerProduct((*I)[i],(*J6)[i]);
				(*K)[7]+=innerProduct((*I)[i],(*J7)[i]);
				(*K)[8]+=innerProduct((*I)[i],(*J8)[i]);
				(*K)[9]+=innerProduct((*I)[i],(*J9)[i]);
				(*K)[10]+=innerProduct((*I)[i],(*J10)[i]);
				(*K)[11]+=innerProduct((*I)[i],(*J11)[i]);
				(*K)[12]+=innerProduct((*I)[i],(*J12)[i]);
				(*K)[13]+=innerProduct((*I)[i],(*J13)[i]);
				(*K)[14]+=innerProduct((*I)[i],(*J14)[i]);
				(*K)[15]+=innerProduct((*I)[i],(*J15)[i]);
				}
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



		Core::Correlator< complex16 > twopoint(res);

		//sum over space
		twopoint.sumOverSpatialVolume(); 
		twopoint.deleteField();

		if (weave.isRoot()) 
		{

			for (size_t i=0; i < 16;i++)
			{
				for (size_t t=0;t < phi.T();t++)
				{
					twopoints.push_back(twopoint[t][i]);
				}
			}

		}

		return twopoints;
	}



}

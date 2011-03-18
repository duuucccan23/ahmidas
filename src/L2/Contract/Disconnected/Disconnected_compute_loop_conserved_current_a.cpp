#include "Disconnected.ih"


//Compute the disconnected loop of the operator O_{mumu} and O_{mumu}\gamma5
// O_00 and g5 O_00 cross checked independtly using an independent code 
// still to cross check the spatial derivative
class complex96
{
	private :
		std::complex<double> data[96];

	public :


		inline complex96() { std::fill_n(data,96,0.0); }
		inline complex96(std::complex<double> c) { std::fill_n(data,96,c); }
		inline std::complex<double> &operator[](size_t index){ return data[index]; }
		inline complex96 &operator+=(complex96 const &other){ std::transform (this->data, this->data + 96, other.data, this->data, std::plus< std::complex< double > >()); return *this; }

		inline complex96 &operator=(complex96 const &other) {
		 if (&other != this)	
			std::copy(other.data, other.data + 96, data); 
		 return *this;
		}
		inline complex96 operator*(std::complex<double> const &factor) const {
			complex96 tmp;

			std::transform (this->data, this->data + 96, tmp.data, std::bind2nd(std::multiplies< std::complex< double > >(), factor)); 
			return tmp;
		}
};



namespace Contract
{
	// compute xi^* x Gamma x psi
	std::vector<  std::complex<double> > compute_loop_conserved_vector_current(
			Core::Field < QCD::Gauge > &gauge_field,
			Core::Propagator const &xi, Core::Propagator const &psi,int const * const position_offset, std::vector< int* > const &momenta,int const tsrc)
	{

		Dirac::Gamma<1> gamma1;
		Dirac::Gamma<2> gamma2;
		Dirac::Gamma<3> gamma3;
		Dirac::Gamma<4> gamma4;

		size_t const L(psi.L());
		size_t const T(psi.T());


		assert(xi.T() == psi.T() && psi.L() == psi.L());
		Base::Weave weave(xi.L(), xi.T());
		std::vector<  std::complex<double>  > twopoints;

		if (weave.isRoot()) std::cout <<"merde" <<std::endl;

		Core::Propagator Gamma_psi(L,T);

		Core::Propagator xi_shifted_T_UP(xi);
		Core::Propagator xi_shifted_T_DOWN(xi);
		Core::Propagator xi_shifted_X_UP(xi);
		Core::Propagator xi_shifted_X_DOWN(xi);
		Core::Propagator xi_shifted_Y_UP(xi);
		Core::Propagator xi_shifted_Y_DOWN(xi);
		Core::Propagator xi_shifted_Z_UP(xi);
		Core::Propagator xi_shifted_Z_DOWN(xi);

		xi_shifted_T_UP.isolate();
		xi_shifted_T_DOWN.isolate();
		xi_shifted_X_UP.isolate();
		xi_shifted_X_DOWN.isolate();
		xi_shifted_Y_UP.isolate();
		xi_shifted_Y_DOWN.isolate();
		xi_shifted_Z_UP.isolate();
		xi_shifted_Z_DOWN.isolate();

		xi_shifted_T_UP.shift(Base::idx_T, Base::dir_UP);
		xi_shifted_T_DOWN.shift(Base::idx_T, Base::dir_DOWN);
		xi_shifted_X_UP.shift(Base::idx_X, Base::dir_UP);
		xi_shifted_X_DOWN.shift(Base::idx_X, Base::dir_DOWN);
		xi_shifted_Y_UP.shift(Base::idx_Y, Base::dir_UP);
		xi_shifted_Y_DOWN.shift(Base::idx_Y, Base::dir_DOWN);
		xi_shifted_Z_UP.shift(Base::idx_Z, Base::dir_UP);
		xi_shifted_Z_DOWN.shift(Base::idx_Z, Base::dir_DOWN);

				Core::Propagator Gamma_psi0(psi);
		Core::Propagator Gamma_psi1(psi);
		Core::Propagator Gamma_psi2(psi);
		Core::Propagator Gamma_psi3(psi);
		Core::Propagator Gamma_psi4(psi);
		Core::Propagator Gamma_psi5(psi);
		Core::Propagator Gamma_psi6(psi);
		Core::Propagator Gamma_psi7(psi);

		//V_0  


		if (weave.isRoot()) std::cout <<"merde2" <<std::endl;
		// part1 tr (g0 -1) U_0(x) psi(x+0) xi^star(x)
		{	Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma4);
			tmp -= tmp2;
			Gamma_psi0 = tmp;
			Gamma_psi0.isolate();
		}


		{// part2 tr (g0+1) U^dag_0(x) psi(x) xi^star(x+0)
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T); 
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma4);
			tmp +=tmp2;
			Gamma_psi1 = tmp;
			Gamma_psi1.isolate();
		}


		//V_1  

		// part1 tr (g1 -1) U_1(x) psi(x+1) xi^star(x)
		{	Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_X, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_X);
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma1);
			tmp -= tmp2;
			Gamma_psi2 = tmp;
			Gamma_psi2.isolate();
		}


		{// part2 tr (g1+1) U^dag_1(x) psi(x) xi^star(x+1))
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_X); 
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma1);
			tmp +=tmp2;
			Gamma_psi3 = tmp;
			Gamma_psi3.isolate();
		}

		//V2

		// part1 tr (g2 -1) U_2(x) psi(x+2) xi^star(x)
		{	Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_Y, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Y);
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma2);
			tmp -= tmp2;
			Gamma_psi4 = tmp;
			Gamma_psi4.isolate();
		}


		{// part2 tr (g2+1) U^dag_2(x) psi(x) xi^star(x+2)
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Y); 
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma2);
			tmp +=tmp2;
			Gamma_psi5 = tmp;
			Gamma_psi5.isolate();
		}

		//V3

		// part1 tr (g3 -1) U_3(x) psi(x+3) xi^star(x)
		{	Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_Z, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Z);
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma3);
			tmp -= tmp2;
			Gamma_psi6 = tmp;
			Gamma_psi6.isolate();
		}


		{// part2 tr (g3+1) U^dag_3(x) psi(x) xi^star(x+3)
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Z); 
			Core::Propagator tmp2(tmp);
			tmp.rightMultiply(gamma3);
			tmp +=tmp2;
			Gamma_psi7 = tmp;
			Gamma_psi7.isolate();
		}

Core::Propagator::const_iterator I(xi.begin());
		Core::Propagator::const_iterator I0(xi_shifted_T_UP.begin());
		Core::Propagator::const_iterator I1(xi_shifted_T_DOWN.begin());
		Core::Propagator::const_iterator I2(xi_shifted_X_UP.begin());
		Core::Propagator::const_iterator I3(xi_shifted_X_DOWN.begin());
		Core::Propagator::const_iterator I4(xi_shifted_Y_UP.begin());
		Core::Propagator::const_iterator I5(xi_shifted_Y_DOWN.begin());
		Core::Propagator::const_iterator I6(xi_shifted_Z_UP.begin());
		Core::Propagator::const_iterator I7(xi_shifted_Z_DOWN.begin());

		Core::Propagator::const_iterator J0(Gamma_psi0.begin());
		Core::Propagator::const_iterator J1(Gamma_psi1.begin());
		Core::Propagator::const_iterator J2(Gamma_psi2.begin());
		Core::Propagator::const_iterator J3(Gamma_psi3.begin());
		Core::Propagator::const_iterator J4(Gamma_psi4.begin());
		Core::Propagator::const_iterator J5(Gamma_psi5.begin());
		Core::Propagator::const_iterator J6(Gamma_psi6.begin());
		Core::Propagator::const_iterator J7(Gamma_psi7.begin());



		//declare iterator on a field of 96 complex.... 
		Core::Field< complex96 > res(xi.L(),xi.T());
		Core::Field< complex96 >::iterator K(res.begin());

		while(I != xi.end())
		{
			for (size_t i=0; i<12;i++)
			{

				(*K)[0+8*i]=innerProduct((*I)[i],(*J0)[i]);
				(*K)[1+8*i]=innerProduct((*I1)[i],(*J1)[i]);
				(*K)[2+8*i]=innerProduct((*I)[i],(*J2)[i]);
				(*K)[3+8*i]=innerProduct((*I3)[i],(*J3)[i]);
				(*K)[4+8*i]=innerProduct((*I)[i],(*J4)[i]);
				(*K)[5+8*i]=innerProduct((*I5)[i],(*J5)[i]);
				(*K)[6+8*i]=innerProduct((*I)[i],(*J6)[i]);
				(*K)[7+8*i]=innerProduct((*I7)[i],(*J7)[i]);

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

//			++I0;
			++I1;
//			++I2;
			++I3;
//			++I4;
			++I5;
//			++I6;
			++I7;
			++I;
		}


		Core::Correlator< complex96 > twopoint(&res);

		twopoint.prepareMomentumProjection(position_offset);


		std::vector< Core::Correlator< complex96 > > tmp(twopoint.momentumProjection(momenta));



		//sum over space
		//		twopoint.sumOverSpatialVolume(); 
		twopoint.deleteField();

		if (weave.isRoot()) 
		{
			for (size_t m=0; m <momenta.size();m++)
			{
				tmp[m].setOffset(tsrc);

				for (size_t i=0; i < 96;i++)
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

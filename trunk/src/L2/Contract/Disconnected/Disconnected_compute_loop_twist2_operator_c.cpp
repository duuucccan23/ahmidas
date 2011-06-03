#include "Disconnected.ih"

// This function use 2 times less memory by splitting the work into 2 parts.
// 2 loops on space are done

//Compute the disconnected loop of the operator O_{mumu} and O_{mumu}\gamma5
// O_00 and g5 O_00 cross checked independtly using an independent code 
// still to cross check the spatial derivative
class complex192 
{
	private :
		std::complex<double> data[192];

	public :


		inline complex192() { std::fill_n(data,192,0.0); }
		inline complex192(std::complex<double> c) { std::fill_n(data,192,c); }
		//inline ~complex192() { delete data ;} 
		inline std::complex<double> &operator[](size_t index){ return data[index]; }
		inline complex192 &operator+=(complex192 const &other){ std::transform (this->data, this->data + 192, other.data, this->data, std::plus< std::complex< double > >()); return *this; }
};



namespace Contract
{

	// compute xi^* x Gamma x psi
	std::vector<  std::complex<double> > compute_loop_twist2_operator(
			Core::Field < QCD::Gauge > &gauge_field,
			Core::Propagator const &xi, Core::Propagator const &psi, bool pol)
	{
		Dirac::Gamma<5> gamma5;
		Dirac::Gamma<1> gamma1;
		Dirac::Gamma<2> gamma2;
		Dirac::Gamma<3> gamma3;
		Dirac::Gamma<4> gamma4;

		


		size_t const L(psi.L());
		size_t const T(psi.T());

		//declare iterator on a field of 192 complex.... 
		Core::Field< complex192 > res(xi.L(),xi.T());
		Core::Field< complex192 >::iterator K(res.begin());



		assert(xi.T() == psi.T() && psi.L() == psi.L());
		Base::Weave weave(xi.L(), xi.T());
		std::vector<  std::complex<double>  > twopoints;

		{

			Core::Propagator Gamma_psi0(psi);
			Core::Propagator Gamma_psi1(psi);
			Core::Propagator Gamma_psi2(psi);
			Core::Propagator Gamma_psi3(psi);
			Core::Propagator Gamma_psi4(psi);
			Core::Propagator Gamma_psi5(psi);
			Core::Propagator Gamma_psi6(psi);
			Core::Propagator Gamma_psi7(psi);

			// now start the very ugly part ... O_mu mu and gamma_5 O_mu mu 

			//O44  

			// part1 tr g0 U_0(x) psi(x+0) xi^star(x)
			{	Core::Propagator tmp(psi);
				tmp.isolate();
				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_T, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_T);
				tmp.rightMultiply(gamma4);

				Gamma_psi0 = tmp;
				//			Gamma_psi0.isolate();
			}


			{
				Core::Propagator tmp(psi);
				tmp.isolate();
				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_T);  
				tmp.shift(Base::idx_T, Base::dir_UP);
				tmp.rightMultiply(gamma4);

				Gamma_psi1= tmp;
				//			Gamma_psi1.isolate();
			}


			{// part3 tr g0 U^dag_0(x) psi(x) xi^star(x+0)
				Core::Propagator tmp(psi);

				tmp.isolate();
				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_T); 
				tmp.rightMultiply(gamma4);

				Gamma_psi2 = tmp;
				//			Gamma_psi2.isolate();
			}



			{// part4 tr g0 U_0(x-0) psi(x) xi^star(x-0)
				Core::Propagator tmp(psi);

				tmp.isolate();
				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_T, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_T);
				tmp.rightMultiply(gamma4);
				tmp.shift(Base::idx_T, Base::dir_UP);

				Gamma_psi3 = tmp;
				//			Gamma_psi3.isolate();

			}

			//end O44

			//O11
			{
				Core::Propagator tmp(psi);

				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);

				tmp.shift(Base::idx_X, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_X);
				tmp.rightMultiply(gamma1);

				Gamma_psi4 = tmp;
				//			Gamma_psi4.isolate();

			}


			{
				Core::Propagator tmp(psi);

				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_X);  
				tmp.shift(Base::idx_X, Base::dir_UP);
				tmp.rightMultiply(gamma1);

				Gamma_psi5= tmp;
				//			Gamma_psi5.isolate();

			}

			{
				Core::Propagator tmp(psi);
				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_X); 
				tmp.rightMultiply(gamma1);

				Gamma_psi6 = tmp;
				//			Gamma_psi6.isolate();

			}

			{
				Core::Propagator tmp(psi);
				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_X, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_X);
				tmp.rightMultiply(gamma1);
				tmp.shift(Base::idx_X, Base::dir_UP);

				Gamma_psi7 = tmp;
				//			Gamma_psi7.isolate();


				//end O11
			}



			Core::Propagator xi_shifted_T_UP(xi);
			Core::Propagator xi_shifted_T_DOWN(xi);
			Core::Propagator xi_shifted_X_UP(xi);
			Core::Propagator xi_shifted_X_DOWN(xi);

			xi_shifted_T_UP.isolate();
			xi_shifted_T_DOWN.isolate();
			xi_shifted_X_UP.isolate();
			xi_shifted_X_DOWN.isolate();

			xi_shifted_T_UP.shift(Base::idx_T, Base::dir_UP);
			xi_shifted_T_DOWN.shift(Base::idx_T, Base::dir_DOWN);
			xi_shifted_X_UP.shift(Base::idx_X, Base::dir_UP);
			xi_shifted_X_DOWN.shift(Base::idx_X, Base::dir_DOWN);


			Core::Propagator::const_iterator I(xi.begin());
			Core::Propagator::const_iterator I0(xi_shifted_T_UP.begin());
			Core::Propagator::const_iterator I1(xi_shifted_T_DOWN.begin());
			Core::Propagator::const_iterator I2(xi_shifted_X_UP.begin());
			Core::Propagator::const_iterator I3(xi_shifted_X_DOWN.begin());


			Core::Propagator::const_iterator J0(Gamma_psi0.begin());
			Core::Propagator::const_iterator J1(Gamma_psi1.begin());
			Core::Propagator::const_iterator J2(Gamma_psi2.begin());
			Core::Propagator::const_iterator J3(Gamma_psi3.begin());
			Core::Propagator::const_iterator J4(Gamma_psi4.begin());
			Core::Propagator::const_iterator J5(Gamma_psi5.begin());
			Core::Propagator::const_iterator J6(Gamma_psi6.begin());
			Core::Propagator::const_iterator J7(Gamma_psi7.begin());



			while(I != xi.end())
			{
				for (size_t i=0; i<12;i++)
				{

					(*K)[0+16*i]=innerProduct((*I)[i],(*J0)[i]); // correct
					(*K)[1+16*i]=innerProduct((*I)[i],(*J1)[i]);
					(*K)[2+16*i]=innerProduct((*I1)[i],(*J2)[i]);  //should be I1 ?
					(*K)[3+16*i]=innerProduct((*I0)[i],(*J3)[i]);
					(*K)[4+16*i]=innerProduct((*I)[i],(*J4)[i]);
					(*K)[5+16*i]=innerProduct((*I)[i],(*J5)[i]);
					(*K)[6+16*i]=innerProduct((*I3)[i],(*J6)[i]);
					(*K)[7+16*i]=innerProduct((*I2)[i],(*J7)[i]);

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

				++I0;
				++I1;
				++I2;
				++I3;
				++I;
			}


		}



		{
			Core::Propagator Gamma_psi8(psi);
			Core::Propagator Gamma_psi9(psi);
			Core::Propagator Gamma_psi10(psi);
			Core::Propagator Gamma_psi11(psi);
			Core::Propagator Gamma_psi12(psi);
			Core::Propagator Gamma_psi13(psi);
			Core::Propagator Gamma_psi14(psi);
			Core::Propagator Gamma_psi15(psi);




			//O22
			{
				Core::Propagator tmp(psi);


				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_Y, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_Y);
				tmp.rightMultiply(gamma2);

				Gamma_psi8 = tmp;
				//			Gamma_psi8.isolate();

			}


			{
				Core::Propagator tmp(psi);

				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_Y);  
				tmp.shift(Base::idx_Y, Base::dir_UP);
				tmp.rightMultiply(gamma2);

				Gamma_psi9= tmp;
				//			Gamma_psi9.isolate();

			}

			{
				Core::Propagator tmp(psi);
				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_Y); 
				tmp.rightMultiply(gamma2);

				Gamma_psi10 = tmp;
				//			Gamma_psi10.isolate();


			}

			{
				Core::Propagator tmp(psi);

				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_Y, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_Y);
				tmp.rightMultiply(gamma2);
				tmp.shift(Base::idx_Y, Base::dir_UP);

				Gamma_psi11 = tmp;
				//			Gamma_psi11.isolate();


			}
			//end O22

			//O33 
			{
				Core::Propagator tmp(psi);

				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_Z, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_Z);
				tmp.rightMultiply(gamma3);

				Gamma_psi12 = tmp;
				//			Gamma_psi12.isolate();

			}


			{
				Core::Propagator tmp(psi);

				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_Z);  
				tmp.shift(Base::idx_Z, Base::dir_UP);
				tmp.rightMultiply(gamma3);

				Gamma_psi13= tmp;
				//			Gamma_psi13.isolate();


			}
			{
				Core::Propagator tmp(psi);
				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.rightMultiplyDagger(gauge_field,Base::idx_Z); 
				tmp.rightMultiply(gamma3);

				Gamma_psi14 = tmp;
				//			Gamma_psi14.isolate();


			}

			{
				Core::Propagator tmp(psi);
				tmp.isolate();

				if (pol==true) tmp.rightMultiply(gamma5);
				tmp.shift(Base::idx_Z, Base::dir_DOWN);
				tmp.rightMultiply(gauge_field,Base::idx_Z);
				tmp.rightMultiply(gamma3);
				tmp.shift(Base::idx_Z, Base::dir_UP);

				Gamma_psi15 = tmp;
				//			Gamma_psi15.isolate();

			}
			//end O33


			Core::Propagator xi_shifted_Y_UP(xi);
			Core::Propagator xi_shifted_Y_DOWN(xi);
			Core::Propagator xi_shifted_Z_UP(xi);
			Core::Propagator xi_shifted_Z_DOWN(xi);

			xi_shifted_Y_UP.isolate();
			xi_shifted_Y_DOWN.isolate();
			xi_shifted_Z_UP.isolate();
			xi_shifted_Z_DOWN.isolate();

			xi_shifted_Y_UP.shift(Base::idx_Y, Base::dir_UP);
			xi_shifted_Y_DOWN.shift(Base::idx_Y, Base::dir_DOWN);
			xi_shifted_Z_UP.shift(Base::idx_Z, Base::dir_UP);
			xi_shifted_Z_DOWN.shift(Base::idx_Z, Base::dir_DOWN);


			Core::Propagator::const_iterator I(xi.begin());
			Core::Propagator::const_iterator I4(xi_shifted_Y_UP.begin());
			Core::Propagator::const_iterator I5(xi_shifted_Y_DOWN.begin());
			Core::Propagator::const_iterator I6(xi_shifted_Z_UP.begin());
			Core::Propagator::const_iterator I7(xi_shifted_Z_DOWN.begin());


			Core::Propagator::const_iterator J8(Gamma_psi8.begin());
			Core::Propagator::const_iterator J9(Gamma_psi9.begin());
			Core::Propagator::const_iterator J10(Gamma_psi10.begin());
			Core::Propagator::const_iterator J11(Gamma_psi11.begin());
			Core::Propagator::const_iterator J12(Gamma_psi12.begin());
			Core::Propagator::const_iterator J13(Gamma_psi13.begin());
			Core::Propagator::const_iterator J14(Gamma_psi14.begin());
			Core::Propagator::const_iterator J15(Gamma_psi15.begin());



			while(I != xi.end())
			{
				for (size_t i=0; i<12;i++)
				{

					(*K)[8+16*i]=innerProduct((*I)[i],(*J8)[i]);
					(*K)[9+16*i]=innerProduct((*I)[i],(*J9)[i]);
					(*K)[10+16*i]=innerProduct((*I5)[i],(*J10)[i]);
					(*K)[11+16*i]=innerProduct((*I4)[i],(*J11)[i]);
					(*K)[12+16*i]=innerProduct((*I)[i],(*J12)[i]);
					(*K)[13+16*i]=innerProduct((*I)[i],(*J13)[i]);
					(*K)[14+16*i]=innerProduct((*I7)[i],(*J14)[i]);
					(*K)[15+16*i]=innerProduct((*I6)[i],(*J15)[i]);

				}

				++K;
				++J8;
				++J9;
				++J10;
				++J11;
				++J12;
				++J13;
				++J14;
				++J15;

				++I4;
				++I5;
				++I6;
				++I7;
				++I;
			}
		}

		Core::Correlator< complex192 > twopoint(res);

		//sum over space
		twopoint.sumOverSpatialVolume(); 
		twopoint.deleteField();

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


		return twopoints;
	}

}

#include "Baryon.ih"

namespace Contract
{
  // this is another version of sequential source, with fixed operator timeslice.
  // It is basically the forward propagator times the operator
  // and the sequential propagator completes the line to the sink which is therefore
  // variable while one inversion is needed for each operator timeslice.

void create_sequential_source_proton_fixed_insertion_timeslice_derivative_operator(Core::Propagator *seqSrc_u,
										   Core::Propagator *seqSrc_d,
										   Core::Field < QCD::Gauge > &gauge_field,
										   Core::Propagator const &u,
										   Core::Propagator const &d,
										   size_t const t_op, Base::Operator op)


 {

    size_t const L(u.L());
    size_t const T(u.T());
	assert(L == d.L() && T == d.T());
	assert(L == seqSrc_u->L() && T == seqSrc_u->T());
	assert(L == seqSrc_d->L() && T == seqSrc_d->T());


	std::complex<double> const z=exp(  std::complex<double>(0,1) * M_PI/double(T)  );


	std::cout << z << std::endl;
	Base::Weave weave(L, T);

	(*seqSrc_u) *= 0.00;
	(*seqSrc_d) *= 0.00;

	Dirac::Gamma< 4 > gamma4;
	Dirac::Gamma< 1 > gamma1;
	Dirac::Gamma< 2 > gamma2;
	Dirac::Gamma< 3 > gamma3;


	for(int i=0;i<2;i++)
	{

		Core::Propagator seqSrc(L,T);
		Core::Propagator prop(L,T);

		prop *= 0.00;
		seqSrc *= 0.00;

		if (i==0) prop += u;
		if (i==1) prop += d;


		// temporal derivative

		{
			Core::Propagator tmp(prop);
			tmp.isolate();
			//  part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.select_timeslice(t_op);
			tmp.rightMultiply(gamma4);
			seqSrc += tmp;
		}


		{
			Core::Propagator tmp(prop);
			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T);  
			tmp.shift(Base::idx_T, Base::dir_UP);
			tmp.select_timeslice(t_op);
			tmp.rightMultiply(gamma4);
			seqSrc -= tmp;
		}

		{
			Core::Propagator tmp(prop);
			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T); 
			tmp.shift(Base::idx_T, Base::dir_UP);
			tmp.rightMultiply(gamma4);
			tmp.select_timeslice(t_op+1);
			tmp *= std::conj(z);
			seqSrc -= tmp;
		}

		{
			Core::Propagator tmp(prop);
			tmp.isolate();
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.select_timeslice(t_op-1);

			tmp.rightMultiply(gamma4);
			tmp *= z;
			seqSrc += tmp;
		}

		switch (op)
		{
			case Base::op_O44_with_substraction:
				{

					std::cout << "in the O44 with substraction case !" << std::endl; 

					// derivative direction 1
					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.shift(Base::idx_X, Base::dir_DOWN);
						tmp.rightMultiply(gauge_field,Base::idx_X);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma1);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc += tmp;
					}


					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.rightMultiplyDagger(gauge_field,Base::idx_X);  
						tmp.shift(Base::idx_X, Base::dir_UP);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma1);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc -= tmp;
					}

					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.rightMultiplyDagger(gauge_field,Base::idx_X); 
						tmp.shift(Base::idx_X, Base::dir_UP);
						tmp.rightMultiply(gamma1);
						tmp.select_timeslice(t_op);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc -= tmp;
					}

					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.shift(Base::idx_X, Base::dir_DOWN);
						tmp.rightMultiply(gauge_field,Base::idx_X);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma1);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc += tmp;
					}


					// derivative direction 2
					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.shift(Base::idx_Y, Base::dir_DOWN);
						tmp.rightMultiply(gauge_field,Base::idx_Y);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma2);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc += tmp;
					}


					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.rightMultiplyDagger(gauge_field,Base::idx_Y);  
						tmp.shift(Base::idx_Y, Base::dir_UP);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma2);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc -= tmp;
					}

					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.rightMultiplyDagger(gauge_field,Base::idx_Y); 
						tmp.shift(Base::idx_Y, Base::dir_UP);
						tmp.rightMultiply(gamma2);
						tmp.select_timeslice(t_op);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc -= tmp;
					}

					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.shift(Base::idx_Y, Base::dir_DOWN);
						tmp.rightMultiply(gauge_field,Base::idx_Y);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma2);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc += tmp;
					}

					// derivative direction 3
					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.shift(Base::idx_Z, Base::dir_DOWN);
						tmp.rightMultiply(gauge_field,Base::idx_Z);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma3);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc += tmp;
					}


					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.rightMultiplyDagger(gauge_field,Base::idx_Z);  
						tmp.shift(Base::idx_Z, Base::dir_UP);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma3);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc -= tmp;
					}

					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.rightMultiplyDagger(gauge_field,Base::idx_Z); 
						tmp.shift(Base::idx_Z, Base::dir_UP);
						tmp.rightMultiply(gamma3);
						tmp.select_timeslice(t_op);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc -= tmp;
					}

					{
						Core::Propagator tmp(prop);
						tmp.isolate();
						tmp.shift(Base::idx_Z, Base::dir_DOWN);
						tmp.rightMultiply(gauge_field,Base::idx_Z);
						tmp.select_timeslice(t_op);
						tmp.rightMultiply(gamma3);
						tmp *= std::complex<double>(-1./3.,0.);
						seqSrc += tmp;
					}

				}
		}

		if (i==0) (*seqSrc_u) += seqSrc;
		if (i==1) (*seqSrc_d) += seqSrc;
	}

 }
}

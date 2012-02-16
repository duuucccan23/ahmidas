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
	std::vector< std::complex<double>  > compute_loop_gluon_Pmunu(Core::Field < QCD::Gauge > &gauge_field)
	{

		clock_t start, finish;
		start = clock();



		Base::Weave weave(gauge_field.L(), gauge_field.T());
		std::vector<  std::complex<double>  > twopoints;
		size_t L = gauge_field.L();
		size_t T = gauge_field.T();



		Core::Correlator< std::complex<double> > P_XY(L,T);
		Core::Correlator< std::complex<double> > P_XZ(L,T);
		Core::Correlator< std::complex<double> > P_XT(L,T);
		Core::Correlator< std::complex<double> > P_YZ(L,T);
		Core::Correlator< std::complex<double> > P_YT(L,T);
		Core::Correlator< std::complex<double> > P_ZT(L,T);

		{
			P_XY = Tool::Plaquette_timeslice(gauge_field, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP);
			//sum over space
			P_XY.sumOverSpatialVolume();
			P_XY.deleteField();
		}

		{
			P_XZ = Tool::Plaquette_timeslice(gauge_field, Base::idx_X, Base::dir_UP, Base::idx_Z, Base::dir_UP);
			//sum over space
			P_XZ.sumOverSpatialVolume();
			P_XZ.deleteField();
		}

		{
			P_XT = Tool::Plaquette_timeslice(gauge_field, Base::idx_X, Base::dir_UP, Base::idx_T, Base::dir_UP);
			//sum over space
			P_XT.sumOverSpatialVolume();
			P_XT.deleteField();
		}

		{
			P_YZ = Tool::Plaquette_timeslice(gauge_field, Base::idx_Y, Base::dir_UP, Base::idx_Z, Base::dir_UP);
			//sum over space
			P_YZ.sumOverSpatialVolume();
			P_YZ.deleteField();
		}

		{
			P_YT = Tool::Plaquette_timeslice(gauge_field, Base::idx_Y, Base::dir_UP, Base::idx_T, Base::dir_UP);
			//sum over space
			P_YT.sumOverSpatialVolume();
			P_YT.deleteField();
		}

		{
			P_ZT = Tool::Plaquette_timeslice(gauge_field, Base::idx_Z, Base::dir_UP, Base::idx_T, Base::dir_UP);
			//sum over space
			P_ZT.sumOverSpatialVolume();
			P_ZT.deleteField();
		}

		if (weave.isRoot()) 
		{

			for (size_t t=0;t < gauge_field.T();t++) 				twopoints.push_back(P_XY[t]);
			for (size_t t=0;t < gauge_field.T();t++) 				twopoints.push_back(P_XZ[t]);
			for (size_t t=0;t < gauge_field.T();t++) 				twopoints.push_back(P_XT[t]);
			for (size_t t=0;t < gauge_field.T();t++) 				twopoints.push_back(P_YZ[t]);
			for (size_t t=0;t < gauge_field.T();t++) 				twopoints.push_back(P_YT[t]);
			for (size_t t=0;t < gauge_field.T();t++) 				twopoints.push_back(P_ZT[t]);

		}


		finish = clock();

		if (weave.isRoot())
			std::cout << "Computation gluon loop in  "<< double(finish - start)/CLOCKS_PER_SEC  << "seconds." << std::endl;


		return twopoints;

	}

}

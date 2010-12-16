#include "Fuzz.ih"

namespace Smear
{
	void Fuzz::accumDirection(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex idx) const
	{

		Core::Field< QCD::Gauge > shifter(field);
		
		Core::Component< QCD::Gauge,  SU3::Matrix > shiftComponent(shifter, idx);
		Core::Component< QCD::Gauge,  SU3::Matrix > fieldComponent(field, idx);
		
		//Core::Component< QCD::Gauge, L, T, SU3::Matrix > shiftComponent(shifter, idx);
		//Core::Component< QCD::Gauge, L, T, SU3::Matrix > fieldComponent(*field, idx);
		for (size_t ctr = 1; ctr < d_length; ++ctr)
		{
			shifter.shift(idx, Base::dir_UP);
			fieldComponent.rightMultiply(shiftComponent);
		}

	}



	void Fuzz::smear(Core::Field< QCD::Gauge  > &field) const
	{
		if (d_length <= 1) // Or we'll be making 3 pointless copies...
			return;

		accumDirection(field, Base::idx_X);
		accumDirection(field, Base::idx_Y);
		accumDirection(field, Base::idx_Z);
	}
}

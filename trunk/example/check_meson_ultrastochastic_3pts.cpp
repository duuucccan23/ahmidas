//this version will contract a sequential meson propagator with the s0 source

#define PropagatorType Core::StochasticPropagator<1>
#define SourceType     Core::StochasticSource<1>

#define UltraStocCase

#include "check_meson_3pts_common.cpp"

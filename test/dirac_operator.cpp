// test code for dirac operator

#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Ahmidas.h>
#include <L0/Base/Weave.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Matrix.h>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>


int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  double const tolerance(1.e-13);

  const size_t L = 4;
  const size_t T = 8;

  size_t const timeslice_boundary(T - 1);
  // size_t const timeslice_source(0);

  const double kappa(0.1);
  const double mu(0.01);

  const std::string filename_gauge("../../test/conf48.0000");
  const std::string filename_prop_u("../../test/source48_u");
  const std::string filename_prop_d("../../test/source48_d");
  const std::string filename_src("../../test/source48");


  // const std::string filename_src_out("../../test/my_test_output");


  Base::Weave weave(L,T);


  std::vector<std::string> sourcefiles;
  std::vector<std::string> propfilesU;
  std::vector<std::string> propfilesD;
  std::vector<std::string> gaugefiles;

  std::vector<std::string> outputfiles;

  gaugefiles.push_back(filename_gauge);


  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss << ".inverted";
    oss.flush();
    propfilesU.push_back(std::string(filename_prop_u).append(oss.str()));
    propfilesD.push_back(std::string(filename_prop_d).append(oss.str()));
  }
  for (int f=0; f<12; f++)
  {
    std::ostringstream oss;
    oss <<  ".";
    oss.fill('0');
    oss.width(2);
    oss << f;
    oss.flush();
    sourcefiles.push_back(std::string(filename_src).append(oss.str()));
    //outputfiles.push_back(std::string(filename_src_out).append(oss.str()));
  }


  Core::Field< QCD::Gauge > gauge_field(L, T);
  Tool::IO::load(&gauge_field, gaugefiles[0], Tool::IO::fileILDG);
  if (weave.isRoot())
    std::cout << "gauge field successfully loaded\n";


  Core::Propagator source(L, T);
  Tool::IO::load(&source, sourcefiles, Tool::IO::fileSCIDAC, 64);
  if (weave.isRoot())
    std::cout << "source successfully loaded\n";

  Core::Propagator uProp(L, T);
  Tool::IO::load(&uProp, propfilesU, Tool::IO::fileSCIDAC, 64);
  if (weave.isRoot())
    std::cout << "u quark propagator successfully loaded\n";

  Core::Propagator dProp(L, T);
  Tool::IO::load(&dProp, propfilesD, Tool::IO::fileSCIDAC, 64);
  if (weave.isRoot())
    std::cout << "d quark propagator successfully loaded\n";


  double const temporal_plaq = Tool::temporalPlaquette(gauge_field);
  double const spatial_plaq  = Tool::spatialPlaquette(gauge_field);


  // note: the boundary conditions are already the desired ones for the sources read
  // therefore we don't need the following two function calls
  // uProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  // dProp.changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);



  std::cout.precision(16);
  if (weave.isRoot())
  {
    std::cout << "temporal plaquette: " << temporal_plaq << std::endl;
    std::cout << "spatial plaquette:  " << spatial_plaq << std::endl;
  }


  // norm of source

  double const srcNorm = source.norm();
  if (weave.isRoot())
    std::cout << "norm of source: " << srcNorm << std::endl;


  double const uPropNorm = uProp.norm();
  double const dPropNorm = dProp.norm();
  if (weave.isRoot())
  {
    std::cout << "norm of u propagator: " << uPropNorm << std::endl;
    std::cout << "norm of d propagator: " << dPropNorm << std::endl;
  }

  // apply Dirac operator to u and d propagator



  Core::Propagator solution_u(uProp.applyDiracOperator(gauge_field, kappa, mu, timeslice_boundary));
  double const solNorm_u = solution_u.norm();
  if (weave.isRoot())
    std::cout << "norm of source calculated from u propagator: " << solNorm_u << std::endl;



  // note that we have a negative mu for the d quark
  Core::Propagator solution_d(dProp.applyDiracOperator(gauge_field, kappa, -mu, timeslice_boundary));
  double const solNorm_d = solution_d.norm();
  if (weave.isRoot())
    std::cout << "norm of source calculated from d propagator: " << solNorm_d << std::endl;


  solution_u -= source;
  solution_d -= source;



  double const diffNorm_u = solution_u.norm();
  if (weave.isRoot())
    std::cout << "norm of difference between real source and source calculated from u propagator:\n" << diffNorm_u << std::endl;


  double const diffNorm_d = solution_d.norm();
  if (weave.isRoot())
    std::cout << "norm of difference between real source and source calculated from d propagator:\n" << diffNorm_d << std::endl;


  double const temporal_plaq_new = Tool::temporalPlaquette(gauge_field);
  double const spatial_plaq_new  = Tool::spatialPlaquette(gauge_field);


  if (fabs(1.0 - temporal_plaq/temporal_plaq_new) > tolerance)
  {
    if (weave.isRoot())
      std::cerr << "ERROR: temporal plaquette has changed from " << temporal_plaq << " to " << temporal_plaq_new << std::endl;
    return EXIT_FAILURE;
  }
  if (fabs(1.0 - spatial_plaq/spatial_plaq_new) > tolerance)
  {
    if (weave.isRoot())
      std::cerr << "ERROR: spatial plaquette has changed from " << spatial_plaq << " to " << spatial_plaq_new << std::endl;
    return EXIT_FAILURE;
  }

  if (diffNorm_u*diffNorm_u > tolerance)
  {
    if (weave.isRoot())
    {
      std::cerr << "ERROR: significant difference between true source and source calculated from u propagator: ";
      std::cerr << diffNorm_u << " > " << tolerance << std::endl;
    }
    return EXIT_FAILURE;
  }

  if (diffNorm_d*diffNorm_d > tolerance)
  {
    if (weave.isRoot())
    {
      std::cerr << "ERROR: significant difference between true source and source calculated from d propagator: ";
      std::cerr << diffNorm_d << " > " << tolerance << std::endl;
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

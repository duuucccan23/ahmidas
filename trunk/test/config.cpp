#include <L0/Ahmidas.h>
#include <L0/Print.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L1/Tool.h>
#include <L1/Tool/IO.h>

int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);
  Core::Field< QCD::Gauge > myfield(4,8);
  Tool::IO::load(&myfield, "../../test/conf.48", Tool::IO::fileILDG);

  double stored = 0.63500557222887;
  double prec = 1e-14;
  double plaqs = Tool::spatialPlaquette(myfield);
  double plaqds = Tool::spatialDownPlaquette(myfield);
  double plaqt = Tool::temporalPlaquette(myfield);
  double plaqdt = Tool::temporalDownPlaquette(myfield);
  double plaq = 0.5 * (plaqt + plaqs);

  std::ostringstream ostr("Spatial plaquette value using UP plaquettes: ", std::ios::ate);
  ostr << std::setprecision(14) << plaqs << std::endl
       << "Spatial plaquette value using DOWN plaquettes: " << plaqds << std::endl
       << "Temporal plaquette value using UP plaquettes: " << plaqt << std::endl
       << "Temporal plaquette value using DOWN plaquettes: " << plaqdt << std::endl
       << "Summed plaquette value: "  << plaq << std::endl
       << "Stored plaquette value: "  << stored;
  Print(ostr.str());
  bool plaqerr = (fabs(plaq/stored - 1) > prec);
  ostr.str("This differs ");
  ostr << (plaqerr ? "more" : "less") << " than " << prec << " from the stored value.";
  Print(ostr.str());

  return (plaqerr);
}

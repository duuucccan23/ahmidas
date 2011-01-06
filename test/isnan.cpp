#include <L0/Base/Base.h>
#include <L0/Print.h>
#include <L0/Ahmidas.h>

int main(int argc, char **argv)
{
  Ahmidas start(&argc, &argv);
  double a = 0;
  bool stdflag = Base::isNaN(a);
  bool infflag = Base::isNaN(1.0/a);
  bool nanflag = Base::isNaN(a/a);
  std::ostringstream ostr("Base::isnan(", std::ios::ate);
  ostr << a << ") = " << stdflag << std::endl
       << "Base::isnan(" << 1/a << ") = " << infflag << std::endl
       << "Base::isnan(" << a/a << ") = " << nanflag;
  Print(ostr.str());
  return !(!stdflag & !infflag & nanflag);
}

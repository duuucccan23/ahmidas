#include <L0/Base/Base.h>
#include <iostream>

int main(int argc, char **argv)
{
  double a = 0;
  bool stdflag = Base::isnan(a);
  bool infflag = Base::isnan(1/a);
  bool nanflag = Base::isnan(a/a);
  std::cout << "Base::isnan(" << a << ") = " << stdflag << std::endl;
  std::cout << "Base::isnan(" << 1/a << ") = " << infflag << std::endl;
  std::cout << "Base::isnan(" << a/a << ") = " << nanflag << std::endl;
  return !(!stdflag & !infflag & nanflag);
}

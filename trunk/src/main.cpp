#include <iostream>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/Base/IO.h>
#include <L1/Path.h>
#include <L0/Core/Component.h>
#include <L0/SU3/Matrix.h>

int main(int argc, char **argv)
{
  Core::Field< QCD::Gauge, 8, 8 > myfield;
  Base::IO::loadILDG(&myfield, "../test/conf.88");
  
  /*
  Core::Field< SU3::Matrix, 8, 8 > result(myfield.component< SU3::Matrix >(Base::idx_X).dagger()); //Store first link hc'ed as result (down direction)
  std::cout << result[0];
  */
  Core::Field< SU3::Matrix, 8, 8 > mystaple = Path::staple(myfield, Base::idx_X, Base::dir_UP, Base::idx_Y, Base::dir_UP);
  SU3::Matrix temp0(myfield[0][Base::idx_X]); //link in +x direction
  SU3::Matrix temp1(myfield[1][Base::idx_Y]); //link in +y direction at +x
  SU3::Matrix temp2((myfield[8][Base::idx_X]).dagger()); //link in -x direction at +x,+y

  SU3::Matrix shift0(myfield[0][Base::idx_X]); //link in +x direction
  myfield.shift(Base::idx_X, Base::dir_DOWN);
  SU3::Matrix shift1(myfield[0][Base::idx_Y]);
  myfield.shift(Base::idx_Y, Base::dir_DOWN);
  myfield.shift(Base::idx_X, Base::dir_UP);
  SU3::Matrix shift2((myfield[0][Base::idx_X]).dagger());

  std::cout << "Staple check, before first mult?\n";
  std::cout << temp0;
  temp0.rightMultiply(temp1);
  std::cout << "Staple check, first mult?\n";
  std::cout << temp0;
  temp0.rightMultiply(temp2.dagger());
  std::cout << "Staple check, second mult?\n";
  std::cout << temp0;

  std::cout << "Staple check, daggered\n";
  std::cout << temp2;
  
  std::cout << "Staple\n";
  std::cout << mystaple[0];

  return 0;
}


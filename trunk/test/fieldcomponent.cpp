#include <iostream>
#include <L0/Core/Field.h>
#include <L0/Core/Component.h>

int main(int argc, char **argv)
{
  Core::Field< double, 4, 8 > field;
  Core::Component< double, 4, 8, double > component = field.component< double >(Base::idx_T);
  return 0;
}

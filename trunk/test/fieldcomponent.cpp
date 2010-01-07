#include <iostream>
#include <L0/Core/Field.h>
#include <L0/Core/Component.h>

int main(int argc, char **argv)
{
  Core::Field< double > field(4, 8);
  Core::Component< double, double > component = field.component< double >(Base::idx_T);
  return 0;
}

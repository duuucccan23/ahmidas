#include <iostream>
#include <L0/Core/Propagator.h>
#include <L1/Tool.h>

int main(int argc, char **argv)
{
  Core::Propagator< 12, 8 > myprop;
  Tool::randomize(&myprop);
  
  return 0;
}

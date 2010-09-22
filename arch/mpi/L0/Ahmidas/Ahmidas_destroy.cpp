#include "Ahmidas.ih"

Ahmidas::~Ahmidas()
{
  MPI::Finalize();
}
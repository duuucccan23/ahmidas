#include "Ahmidas.ih"

Ahmidas::~Ahmidas()
{
  MPI_Finalize();
}
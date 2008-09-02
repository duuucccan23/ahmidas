/* Privately accessible constructor */

#include "../Com.h"
#include <cstdlib>
#include <mpi.h>

namespace Core
{
  Com::Com()
  {
    if (!MPI::Is_initialized())
      exit(1);
    MPI::Datatype su3mat = MPI::DOUBLE_COMPLEX.Create_contiguous(9); //just for now!
    su3mat.Commit();
  }
}

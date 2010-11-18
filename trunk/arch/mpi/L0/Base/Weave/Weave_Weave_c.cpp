#include "Weave.ih"

Base::Weave::~Weave()
{
//   s_grid->grid().Barrier();
//   std::cout << "rank = " << rank() << ", this = " << this << ", s_grid_ref = " << s_grid_ref << std::endl;
//   s_grid->grid().Barrier();
  d_grid = NULL;
  s_grid_ref--;
  if(s_grid_ref == 0)
  {
    delete reinterpret_cast< Base::Grid *>(s_grid);
    s_grid = NULL;
//     MPI::Finalize();
//     std::cout << "MPI finalized" << std::endl;
  }
}
#ifndef COM_H
#define COM_H

/* MPI Communication handling */

namespace Core
{
  template< typename Element >
  class Com
  {
    Grid< L, T>   &d_grid;
    MPI::Datatype  d_surfaces[4];

    public:
      Com(Com const &other);
      Com(Core::Grid< L, T > &grid, Element const &value);

      Com< Element > &operator=(Com< Element > const &other);

      ~Com();

  };

}
#endif

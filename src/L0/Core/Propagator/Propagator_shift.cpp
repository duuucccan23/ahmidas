#include "Propagator.ih"

namespace Core
{
  Propagator &Propagator::shift(Base::SpaceTimeIndex const idx, Base::Direction const dir, size_t const times)
  {
    if (times == 0)
      return *this;

    Base::Direction my_dir(dir);
    if(times < 0)
    {
      switch (dir)
      {
        case Base::dir_UP:
          my_dir = Base::dir_DOWN;
          break;
        case Base::dir_DOWN:
          my_dir = Base::dir_UP;
          break;
        default:
          std::cerr << "Base::Direction not known in "
                    << "Propagator &Propagator::shift(Base::SpaceTimeIndex, Base::Direction , int)\n"
                    << "Aborting..."
                    << std::endl;
          exit(1);
       }
    }

    size_t my_times = size_t(abs(times)) % T();
    // make use of periodic boundary conditions to save time
    if (my_times > T()/2)
    {
      my_times = T()-my_times;
      switch (dir)
      {
        case Base::dir_UP:
          my_dir = Base::dir_DOWN;
          break;
        case Base::dir_DOWN:
          my_dir = Base::dir_UP;
          break;
        default:
          std::cerr << "Base::Direction not known in "
                    << "Propagator &Propagator::shift(Base::SpaceTimeIndex, Base::Direction , int)\n"
                    << "Aborting..."
                    << std::endl;
          exit(1);
       }
    }
    for (size_t I=0; I<my_times; I++)
      d_components->shift(idx, my_dir);
    return *this;
  }
}


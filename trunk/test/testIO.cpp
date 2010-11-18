#include <iomanip>
#include <iostream>

#include <L0/Ahmidas.h>
#include <L0/Base/Base.h>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L1/Tool/IO.h>

// #define __SAVE_FIELD__
#define __LOAD_FIELD__


int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);


  size_t const L(6);
  size_t const T(2);

  std::string file("../../test/my_test_conf");

  Base::Weave weave(L, T);

#ifdef __SAVE_FIELD__
  {
    Core::Field< double > myfield(L, T);

    size_t localIndex;
    long unsigned int count(0);

    for(size_t idx_T = 0; idx_T < T; idx_T++)
    {
      for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
      {
        for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
        {
          for(size_t idx_X = 0; idx_X < L; idx_X++)
          {
            localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
            /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
            if (localIndex == weave.localVolume())
            {
              count++;
              continue;
            }
            myfield[localIndex] = double(count++);
          }
        }
      }
    }

    Tool::IO::saveScidac(myfield, file);
    if (weave.isRoot())
      std::cout << "field created and saved." << std::endl;
    weave.barrier();


  }
#endif
#ifdef __LOAD_FIELD__
  {

    if (weave.isRoot())
      std::cout << "loading field ..." << std::flush;
    weave.barrier();

    Core::Field< double > yourfield(L, T);

    Tool::IO::loadScidac(&yourfield, file);

    if (weave.isRoot())
      std::cout << "done." << std::endl;
    weave.barrier();


    size_t localIndex;
    long unsigned int count(0);

    for(size_t idx_T = 0; idx_T < T; idx_T++)
    {
      for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
      {
        for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
        {
          for(size_t idx_X = 0; idx_X < L; idx_X++)
          {
            localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
            /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
            if (localIndex == weave.localVolume())
            {
              count++;
              continue;
            }
            std::cout << yourfield[localIndex] << " <-> " << count << std::endl;
            assert (yourfield[localIndex] == double(count));
            count++;
          }
        }
      }
    }

    weave.barrier();

    if (weave.isRoot())
      std::cout << "loading works!" << std::endl;
    weave.barrier();

  }
#endif

    weave.barrier();
    if (weave.isRoot())
      std::cout << "everything ok!" << std::endl;
    weave.barrier();

  return EXIT_SUCCESS;
}

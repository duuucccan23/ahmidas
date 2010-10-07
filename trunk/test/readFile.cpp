
#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>

#include <L2/Input/FileReader.h>

int main(int argc, char **argv)
{

  size_t L = 0;
  size_t T = 0;

  Input::FileReader reader("../../test/read_input_file_input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L, T, files, floats, positions, operators);

  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;
  std::cout << "(should be 4x4x4x8)" << std::endl;

  if (L != 4)
  {
    std::cerr << "wrong value for L!";
    return EXIT_FAILURE;
  }
  if (T != 8)
  {
    std::cerr << "wrong value for T !";
    return EXIT_FAILURE;
  }

  double kappa = floats["kappa"];
  double mu = floats["mu"];

  std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;
  std::cout << "(should be 0.15 and 0.01)" << std::endl;

  if (kappa != 0.15)
  {
    std::cerr << "wrong kappa value!" << std::endl;
    return EXIT_FAILURE;
  }
  if (mu != 0.01)
  {
    std::cerr << "wrong mu value!" << std::endl;
    return EXIT_FAILURE;
  }

  size_t const *position = positions[0];
  std::cout << "test source position: (" << position[0] << ", " << position[1] << ", "
                                         << position[2] << ", " << position[3] << ")" << std::endl;
  std::cout << "(should be (2, 0, 3, 1))" << std::endl;

  if (position[0] != 2 || position[1] != 0 || position[2] != 3 || position[3] != 1)
  {
    std::cerr << "source position not correct!" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "files.size() is " << files.size() << std::endl;
  std::cout << "(should be 3)" << std::endl;

  if (files.size() != 3)
  {
    std::cerr << "wrong number of file containers" << std::endl;
    return EXIT_FAILURE;
  }

  if ((files[0])[ 0] != "../../test/source4x4_d.00.inverted" ||
      (files[0])[ 1] != "../../test/source4x4_d.01.inverted" ||
      (files[0])[ 2] != "../../test/source4x4_d.02.inverted" ||
      (files[0])[ 3] != "../../test/source4x4_d.03.inverted" ||
      (files[0])[ 4] != "../../test/source4x4_d.04.inverted" ||
      (files[0])[ 5] != "../../test/source4x4_d.05.inverted" ||
      (files[0])[ 6] != "../../test/source4x4_d.06.inverted" ||
      (files[0])[ 7] != "../../test/source4x4_d.07.inverted" ||
      (files[0])[ 8] != "../../test/source4x4_d.08.inverted" ||
      (files[0])[ 9] != "../../test/source4x4_d.09.inverted" ||
      (files[0])[10] != "../../test/source4x4_d.10.inverted" ||
      (files[0])[11] != "../../test/source4x4_d.11.inverted")
  {
    std::cerr << "last file name in first container is " << (files[0])[11] << std::endl;
    std::cerr << "(should be ../../test/source4x4_d.11.inverted)" << std::endl;
    std::cerr << "there could be another file name that is not correct" << std::endl;
    return EXIT_FAILURE;
  }

  if ((files[1])[ 0] != "../../test/source4x4_u.00.inverted" ||
      (files[1])[ 1] != "../../test/source4x4_u.01.inverted" ||
      (files[1])[ 2] != "../../test/source4x4_u.02.inverted" ||
      (files[1])[ 3] != "../../test/source4x4_u.03.inverted" ||
      (files[1])[ 4] != "../../test/source4x4_u.04.inverted" ||
      (files[1])[ 5] != "../../test/source4x4_u.05.inverted" ||
      (files[1])[ 6] != "../../test/source4x4_u.06.inverted" ||
      (files[1])[ 7] != "../../test/source4x4_u.07.inverted" ||
      (files[1])[ 8] != "../../test/source4x4_u.08.inverted" ||
      (files[1])[ 9] != "../../test/source4x4_u.09.inverted" ||
      (files[1])[10] != "../../test/source4x4_u.10.inverted" ||
      (files[1])[11] != "../../test/source4x4_u.11.inverted")
  {
    std::cerr << "last file name in second container is " << (files[1])[11] << std::endl;
    std::cerr << "(should be ../../test/source4x4_u.11.inverted)" << std::endl;
    std::cerr << "there could be another file name that is not correct" << std::endl;
    return EXIT_FAILURE;
  }

  if ((files[2])[ 0] != "../../test/source4x4.00" ||
      (files[2])[ 1] != "../../test/source4x4.01" ||
      (files[2])[ 2] != "../../test/source4x4.02" ||
      (files[2])[ 3] != "../../test/source4x4.03" ||
      (files[2])[ 4] != "../../test/source4x4.04" ||
      (files[2])[ 5] != "../../test/source4x4.05" ||
      (files[2])[ 6] != "../../test/source4x4.06" ||
      (files[2])[ 7] != "../../test/source4x4.07" ||
      (files[2])[ 8] != "../../test/source4x4.08" ||
      (files[2])[ 9] != "../../test/source4x4.09" ||
      (files[2])[10] != "../../test/source4x4.10" ||
      (files[2])[11] != "../../test/source4x4.11")
  {
    std::cerr << "last file name in third container is " << (files[2])[11] << std::endl;
    std::cerr << "(should be ../../test/source4x4.11)" << std::endl;
    std::cerr << "there could be another file name that is not correct" << std::endl;
    return EXIT_FAILURE;
  }

  /* some operator test should be performed here as well */


  return EXIT_SUCCESS;
}


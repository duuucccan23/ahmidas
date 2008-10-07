#include <iostream>
#include <L0/Core/Field.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Field< double, 4, 4 > field;
  size_t idx[4] = {0, 0, 0, 0};
  for (size_t ctr = 0; ctr < 4; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 4 increments: " << idx[Base::idx_X] << "  " << idx[Base::idx_Y] << "  "
                                      << idx[Base::idx_Z] << "  " << idx[Base::idx_T] << std::endl;
  if (!(idx[Base::idx_X] == 0 && idx[Base::idx_Y] == 1 && idx[Base::idx_Z] == 0 && idx[Base::idx_T] == 0))
    return 1;
  for (size_t ctr = 0; ctr < 12; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 16 increments: " << idx[Base::idx_X] << "  " << idx[Base::idx_Y] << "  "
                                       << idx[Base::idx_Z] << "  " << idx[Base::idx_T] << std::endl;
  if (!(idx[Base::idx_X] == 0 && idx[Base::idx_Y] == 0 && idx[Base::idx_Z] == 1 && idx[Base::idx_T] == 0))
    return 1;
  for (size_t ctr = 0; ctr < 48; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 64 increments: " << idx[Base::idx_X] << "  " << idx[Base::idx_Y] << "  "
                                       << idx[Base::idx_Z] << "  " << idx[Base::idx_T] << std::endl;
  if (!(idx[Base::idx_X] == 0 && idx[Base::idx_Y] == 0 && idx[Base::idx_Z] == 0 && idx[Base::idx_T] == 1))
    return 1;
  for (size_t ctr = 0; ctr < 192; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 256 increments: " << idx[Base::idx_X] << "  " << idx[Base::idx_Y] << "  "
                                        << idx[Base::idx_Z] << "  " << idx[Base::idx_T] << std::endl;
  if (!(idx[Base::idx_X] == 0 && idx[Base::idx_Y] == 0 && idx[Base::idx_Z] == 0 && idx[Base::idx_T] == 4))
    return 1;
  return 0;
}

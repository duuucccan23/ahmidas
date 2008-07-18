#include <iostream>
#include <l0/Core/Field/Field.h>

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);
  Core::Grid< 4, 4 > grid;
  Core::Field< double, 4, 4 > field(grid);
  size_t idx[4] = {0, 0, 0, 0};
  for (size_t ctr = 0; ctr < 4; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 4 increments: " << idx[Core::idx_X] << "  " << idx[Core::idx_Y] << "  "
				      << idx[Core::idx_Z] << "  " << idx[Core::idx_T] << std::endl;
  if (!(idx[Core::idx_X] == 0 && idx[Core::idx_Y] == 1 && idx[Core::idx_Z] == 0 && idx[Core::idx_T] == 0))
    return 1;
  for (size_t ctr = 0; ctr < 12; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 16 increments: " << idx[Core::idx_X] << "  " << idx[Core::idx_Y] << "  "
				       << idx[Core::idx_Z] << "  " << idx[Core::idx_T] << std::endl;
  if (!(idx[Core::idx_X] == 0 && idx[Core::idx_Y] == 0 && idx[Core::idx_Z] == 1 && idx[Core::idx_T] == 0))
    return 1;
  for (size_t ctr = 0; ctr < 48; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 64 increments: " << idx[Core::idx_X] << "  " << idx[Core::idx_Y] << "  "
				        << idx[Core::idx_Z] << "  " << idx[Core::idx_T] << std::endl;
  if (!(idx[Core::idx_X] == 0 && idx[Core::idx_Y] == 0 && idx[Core::idx_Z] == 0 && idx[Core::idx_T] == 1))
    return 1;
  for (size_t ctr = 0; ctr < 192; ++ctr)
    field.increaseIdx(idx);
  std::cout << "After 256 increments: " << idx[Core::idx_X] << "  " << idx[Core::idx_Y] << "  "
				        << idx[Core::idx_Z] << "  " << idx[Core::idx_T] << std::endl;
  if (!(idx[Core::idx_X] == 0 && idx[Core::idx_Y] == 0 && idx[Core::idx_Z] == 0 && idx[Core::idx_T] == 4))
    return 1;
  return 0;
}

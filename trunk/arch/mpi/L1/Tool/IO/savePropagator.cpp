#include <L1/Tool/IO.h>

void Tool::IO::savePropagator(Core::Field< QCD::Spinor > const &field, std::string const &filename)
{
  MPI_File fp;
  Base::Weave weave(field.L(), field.T());
  char *file = new char[filename.length() + 1];
  std::copy(filename.begin(), filename.end() + filename.length(), file);
  file[filename.length()] = '\0';
  MPI_File_open(weave.grid(), file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);
  delete[] file;
  LemonWriter* writer = lemonCreateWriter(&fp, weave.grid());

  savePropagatorType(writer);
  savePropagatorInfo(writer);
  saveScidacBinary(writer, field);

  lemonDestroyWriter(writer);
  MPI_File_close(&fp);
}

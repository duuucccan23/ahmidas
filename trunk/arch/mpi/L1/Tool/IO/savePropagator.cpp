#include <L1/Tool/IO.h>

void Tool::IO::savePropagator(Core::Propagator const &prop, std::vector< std::string > const &filename)
{
  if (filename.size() < 12)
    exit(1);
  MPI_File fp;
  Core::Field< QCD::Tensor > const &field = prop.components(); // Hack to get rid of yet another layer of cr*p...
  Base::Weave weave(field.L(), field.T());
  char *file = new char[256];
  file[255] = '\0'; // Ensure proper termination
  int dims[] = {field.L(), field.L(), field.L(), field.T()};
  for (size_t idx = 0; idx < 12; ++idx)
  {
    strncpy(file, filename[idx].c_str(), 255);

    MPI_File_open(weave.grid(), file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);

    LemonWriter* writer = lemonCreateWriter(&fp, weave.grid());

    savePropagatorType(writer);
    savePropagatorInfo(writer, field.L(), field.T());

    QCD::Spinor *buffer = new QCD::Spinor[field.localVolume()];
    Core::Component< QCD::Tensor, QCD::Spinor > comp = field.component< QCD::Spinor >(idx);

    for (size_t ctr = 0; ctr < field.localVolume(); ++ctr)
      buffer[ctr] = comp[ctr];

    saveScidacBinary(writer, buffer, field.volume(), field.localVolume(), dims); // Takes care of BE issues

    lemonDestroyWriter(writer);
    MPI_File_close(&fp);
  }
  delete[] file;
}

void Tool::IO::savePropagator(Core::StochasticPropagator< 1 > const &prop, std::vector< std::string > const &filename)
{
  if (filename.size() < 1)
    exit(1);
  MPI_File fp;
  Core::Field< QCD::Tensor > const &field = prop.components(); // Hack to get rid of yet another layer of cr*p...
  Base::Weave weave(field.L(), field.T());
  char *file = new char[256];
  file[255] = '\0'; // Ensure proper termination
  int dims[] = {field.L(), field.L(), field.L(), field.T()};
  strncpy(file, filename[0].c_str(), 255);

  MPI_File_open(weave.grid(), file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);

  LemonWriter* writer = lemonCreateWriter(&fp, weave.grid());

  savePropagatorType(writer);
  savePropagatorInfo(writer, field.L(), field.T());

  QCD::Spinor *buffer = new QCD::Spinor[field.localVolume()];
  Core::Component< QCD::Tensor, QCD::Spinor > comp = field.component< QCD::Spinor >(0);

  for (size_t ctr = 0; ctr < field.localVolume(); ++ctr)
    buffer[ctr] = comp[ctr];

  saveScidacBinary(writer, buffer, field.volume(), field.localVolume(), dims); // Takes care of BE issues

  lemonDestroyWriter(writer);
  MPI_File_close(&fp);
  delete[] file;
}

void Tool::IO::savePropagator(Core::StochasticPropagator< 4 > const &prop, std::vector< std::string > const &filename)
{
  if (filename.size() < 4)
    exit(1);
  MPI_File fp;
  Core::Field< QCD::Tensor > const &field = prop.components(); // Hack to get rid of yet another layer of cr*p...
  Base::Weave weave(field.L(), field.T());
  char *file = new char[256];
  file[255] = '\0'; // Ensure proper termination
  int dims[] = {field.L(), field.L(), field.L(), field.T()};
  for (size_t idx = 0; idx < 12; idx += 3)
  {
    strncpy(file, filename[idx / 3].c_str(), 255);

    MPI_File_open(weave.grid(), file, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);

    LemonWriter* writer = lemonCreateWriter(&fp, weave.grid());

    savePropagatorType(writer);
    savePropagatorInfo(writer, field.L(), field.T());

    QCD::Spinor *buffer = new QCD::Spinor[field.localVolume()];
    Core::Component< QCD::Tensor, QCD::Spinor > comp = field.component< QCD::Spinor >(idx);

    for (size_t ctr = 0; ctr < field.localVolume(); ++ctr)
      buffer[ctr] = comp[ctr];

    saveScidacBinary(writer, buffer, field.volume(), field.localVolume(), dims); // Takes care of BE issues

    lemonDestroyWriter(writer);
    MPI_File_close(&fp);
  }
  delete[] file;
}

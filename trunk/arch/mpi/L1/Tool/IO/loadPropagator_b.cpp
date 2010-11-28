#include <L1/Tool/IO.h>

void Tool::IO::loadPropagator(Core::Field< QCD::Spinor > *field, size_t amount, std::string const &filename)
{
  MPI_File fp;
  Base::Weave weave(field->L(), field->T());
  char *file = new char[filename.length() + 1];
  std::copy(filename.begin(), filename.end() + filename.length(), file);
  file[filename.length()] = '\0';
  MPI_File_open(weave.grid(), file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
  delete[] file;
  LemonReader* reader = lemonCreateReader(&fp, weave.grid());

  int status;
  ScidacInfo info;
  info.precision = "64";

  size_t fieldsRead = 0;

  while ((status =lemonReaderNextRecord(reader)) != LEMON_EOF)
  {
    if (status != LEMON_SUCCESS)
    {
      fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
      break;
    }
    char const *headerType = lemonReaderType(reader);

    if (strcmp("etmc-propagator-format", headerType) == 0)
      loadScidacInfo(reader, info);
    if (strcmp("scidac-binary-data", headerType) == 0)
    {
      loadScidacBinary(reader, field[fieldsRead], info);
      ++fieldsRead;
      if (fieldsRead == amount)
	break;
    }

    lemonReaderCloseRecord(reader);
  }
  lemonDestroyReader(reader);
  MPI_File_close(&fp);
}

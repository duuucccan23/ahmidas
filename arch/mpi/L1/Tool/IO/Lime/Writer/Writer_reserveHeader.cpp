#include "Writer.ih"

void Tool::IO::Lime::Writer::reserveHeader()
{

  // this routine is more or less redundant for parallel I/O, but still the metadata is written serially


  if (!d_writeHeader)
  {
    return;
  }

  // d_record.offset should be zero, and position in file should be d_record.recOffset
  assert(d_MPI_FILE.Get_position() == d_record.recOffset);
  assert(d_record.offset == MPI::Offset(0));

  //static int cnt(0);

  // std::cout << "now I am going to write a header of zeros! ("  << cnt ++ << ")" << std::endl;

  for (size_t ctr = 0; ctr < s_headerSize / 8; ++ctr)
    d_MPI_FILE.Write(s_padding, 8, MPI::BYTE);
  //d_MPI_FILE.flush();

  d_hasWritten = false;
}

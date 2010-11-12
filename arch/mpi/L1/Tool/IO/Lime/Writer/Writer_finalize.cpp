#include "Writer.ih"

void Tool::IO::Lime::Writer::finalize()
{

  if(!d_writeHeader || d_hasWritten)
  {
    //d_MPI_FILE.flush();
    // this process has done its duty
    return;
  }

  assert(d_MPI_FILE.Get_position() >= d_record.recOffset);
  // size can be zero if only one process has written the current record
  // in that case we can ask the stream how many bytes have been written
  uint64_t written = 0;
  if (d_MPI_FILE.Get_position() > 0)
    written = d_MPI_FILE.Get_position() - d_record.recOffset - MPI::Offset(s_headerSize);
  //std::cout << "d_MPI_FILE.tellp() = " << d_MPI_FILE.Get_position() << ", d_record.recOffset = " << d_record.recOffset << std::endl;
  //std::cout << "d_record.size = " << d_record.size << ", written = " << written << std::endl;
  // in the other case the actual record size had to be passed and is stored in d_record.size
  written = d_record.size > written ? d_record.size : written;

  if ( written % 8 != 0)
    d_MPI_FILE.Write(s_padding, (8 - (written % 8)) % 8, MPI::BYTE);

  d_startOfNextRecord = d_MPI_FILE.Get_position();

  uHeader header;
  std::fill(header.as8, header.as8 + s_headerSize, 0x00);

  header.as32[0] = s_limeMagic;
  if (!Base::bigEndian)
    Base::swapEndian(header.as32[0]);

  header.as16[2] = d_record.version;
  if (!Base::bigEndian)
    Base::swapEndian(header.as16[2]);

  header.as8[6] = d_record.mesBeg ? s_mesBeginMask : 0x00;
  header.as8[6] = header.as8[6] | (d_record.mesEnd ? s_mesEndMask : 0x00);

  header.as64[1] = written;
  if (!Base::bigEndian)
    Base::swapEndian(header.as64[1]);

  std::copy(d_record.type, d_record.type + 128, header.as8 + 16);

  // Go back and do the actual writing

  // offset is supposed to be zero or s_headerSize,
  // depending on whether Writer::seekp(streampos const) has been called
  // not the case anymore for writing of non-contiguous data
  // assert (d_record.offset == std::streampos(0) || d_record.offset == std::streampos(s_headerSize));

  //std::cout << "d_record.recOffset = " << d_record.recOffset << std::endl;
  //std::cout<<"position before writing the record end: "<<d_MPI_FILE.Get_position()<<std::endl;
  d_MPI_FILE.Seek(d_record.recOffset, MPI_SEEK_SET);
  //std::cout<<"position after seeking record end: "<<d_MPI_FILE.Get_position()<<std::endl;
  d_MPI_FILE.Write(header.as8, s_headerSize,MPI::BYTE);
  //std::cout<<"position after writing record end: "<<d_MPI_FILE.Get_position()<<std::endl;
  assert(good());

  //d_MPI_FILE.flush();

  //static int cnt(0);

//   std::cout << "now I am going to write the actual header! ("  << cnt ++ << ")" << std::endl;
//   std::cout << "---" << std::endl;
//   std::cout << std::string(header.as8, s_headerSize) << std::endl;
//   std::cout << "---" << std::endl;

  d_MPI_FILE.Seek(d_startOfNextRecord, MPI_SEEK_SET);
  d_hasWritten = true;
}

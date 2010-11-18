#include "Reader.ih"

Tool::IO::Lime::Reader::Reader(Base::Weave * weave, std::string const &filename)
  : d_weave(weave), d_filename(filename), d_currentRecord(0),
    d_fail(false), d_messagesCorrect(true)
{
  uHeader header;
  int32_t message = -1;
  bool messageOpened = true;

  d_MPI_FILE = MPI::File::Open(d_weave->d_grid->grid(), filename.c_str(),
                               MPI::MODE_RDONLY, MPI::INFO_NULL);

  if (d_weave->rank() == 0)
  {
    while (good())
    {
      d_MPI_FILE.Read(header.as8, s_headerSize, MPI_BYTE, d_MPI_STATUS);

      if (!Base::bigEndian)
        Base::swapEndian(header.as32[0]);

      if (header.as32[0] != s_limeMagic)
        break; // We're done, apparently.

      if (!Base::bigEndian)
        Base::swapEndian(header.as16[2]);
      d_versions.push_back(header.as16[2]);

      if (!messageOpened)
      {
        if (header.as8[6] & s_mesBeginMask)
        {
          messageOpened = true;
          ++message;
        }
        else
          d_messagesCorrect = false;
      }

      if (messageOpened)
      {
        if (header.as8[6] & s_mesEndMask)
          messageOpened = false;
        else if (header.as8[6] & s_mesBeginMask)
          d_messagesCorrect = false;
      }

      d_messages.push_back(message);

      if (!Base::bigEndian)
        Base::swapEndian(header.as64[1]);
      d_sizes.push_back(header.as64[1]);

      d_types.push_back(std::string(&header.as8[16]));

      d_offsets.push_back(d_MPI_FILE.Get_position());
      d_MPI_FILE.Seek(((d_sizes.back() + 7) / 8) * 8 , MPI_SEEK_CUR);
    }
    if (!d_offsets.empty())
    {
      // here we need to set some error status or exit
      d_fail = false;
      d_MPI_FILE.Seek(d_offsets.front(), MPI_SEEK_SET);
    }
  }
  d_weave->barrier();
}

#include "Reader.ih"

IO::Lime::Reader::Reader(std::string const &filename)
  : d_filename(filename), d_stream(0), d_reader(0), d_currentRecord(0), d_fail(0)
{
  {
    struct stat statusBuffer;
    int readStatus = stat(filename.c_str(), &statusBuffer);

    if (readStatus)
    {
      std::cerr << "An error occurred while attempting to query " << filename << '.' << std::endl;
      MPI::COMM_WORLD.Abort(readStatus);
    }

    d_stream = fopen(filename.c_str(), "r");
    d_reader = limeCreateReader(d_stream);

    if (!d_stream || !d_reader)
      MPI::COMM_WORLD.Abort(EIO);
  }

  
  int32_t messageCtr = -1; // This would indicate a missing message start
  while (true)
  {
    d_recordOffsets.push_back(limeGetReaderPointer(d_reader)); 
    std::cerr << "[DEBUG] offset: " << d_recordOffsets.back() << std::endl;

    d_fail = limeReaderNextRecord(d_reader);
    if (d_fail != LIME_SUCCESS)
      break;
    if (limeReaderMBFlag(d_reader))
      ++messageCtr;

    // We now build up a catalogue of meta data.
    // This data is only available serially in the file,
    // but we'd like random access to it.
    d_messageIndices.push_back(messageCtr);
    std::cerr << "[DEBUG] message: " << d_messageIndices.back() << std::endl;
    d_limeTypes.push_back(limeReaderType(d_reader));
    std::cerr << "[DEBUG] type: " << d_limeTypes.back() << std::endl;
    d_recordSizes.push_back(limeReaderBytes(d_reader));
    std::cerr << "[DEBUG] size: " << d_recordSizes.back() << std::endl;
    std::cerr << "[DEBUG] ***-----------------*** " << std::endl;
  }

  if (d_fail != LIME_EOF)
  {
    std::cerr << "File " << filename << " gave an I/O error." << std::endl;
    MPI::COMM_WORLD.Abort(EIO);
  }

  d_fail = limeSetReaderPointer(d_reader, d_recordOffsets[0]);
}

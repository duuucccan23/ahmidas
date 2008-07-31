namespace Core
{
  template< typename Element, size_t L, size_t T >
  template< >
  void Field< Element, L, T >::readDataFromIO< IO::ILDG::Generic >(IO::ILDG::Generic &inputIO)
  {
    size_t written = 0; //Number of locally written elements
    size_t const elementSize = sizeof(Element) / sizeof(double); // We assume Element is natively a collection of doubles
    size_t const precIO = inputIO.precision() / 8;
  
    // We process all input as chars, to avoid any unpleasantness with floating point registers.
    if (d_grid.rank()) // Prepare to receive
    {
      char *fileBuffer = new char[elementSize * precIO * d_grid.contiguousBlock()];
      for (size_t ctr = 0; ctr < (d_grid.localVolume() / d_grid.contiguousBlock()); ++ctr)
      {
        d_grid.grid().Recv(fileBuffer, d_grid.contiguousBlock() * elementSize * precIO, 
                          MPI::BYTE, 0, TAG_FILE_DISTRIBUTION);
        written = moveBufferToData(fileBuffer, written, precIO);
      }
      delete[] fileBuffer;
      return;
    }
  
      // At this point, we know we're node 0
      // Check the content of this file for a simple sanity condition -- is there enough data to begin with.
    if (d_grid.totalVolume() * elementSize * precIO != confFile.size())
    {
      std::cerr << "Content is of inappropriate size.\nAborting." << std::endl;
      MPI::COMM_WORLD.Abort(EIO);
    }
      
    char *fileBuffer = new char[elementSize * precIO * d_grid.contiguousBlock()];
  
    size_t nBlocks = d_grid.totalVolume() / d_grid.contiguousBlock();
    for (size_t ctr = 0; ctr < nBlocks; ++ctr)
    {
      confFile.read(fileBuffer, d_grid.contiguousBlock() * elementSize * precIO);
      if (confFile.fail())
      {
        std::cerr << "Unexpected error while reading file " << fileName << ".\nAborting." << std::endl;
        MPI::COMM_WORLD.Abort(EIO);
      }
  
      size_t destination = d_grid.rank(ctr * d_grid.contiguousBlock());
      if (!destination) // This block of data should be stored locally
      {
        written = moveBufferToData(fileBuffer, written, precIO);
        continue;
      }
      d_grid.grid().Send(fileBuffer, d_grid.contiguousBlock() * elementSize * precIO, 
                        MPI::BYTE, destination, TAG_FILE_DISTRIBUTION);
    }
    delete[] fileBuffer;
  }
}

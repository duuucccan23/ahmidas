#include <iostream>
#include <mpi.h>
#include <Lime/Reader/Reader.h>
#include <QCD/Gauge/Gauge.h>
#include <Core/Core.h>
#include <Core/Grid/Grid.h>
#include <Core/Field/Field.h>

using namespace std;

int main(int argc, char **argv)
{
  MPI::Init(argc, argv);

  Core::Grid<24, 48> grid;
  size_t myrank = grid.rank();
 
  if (myrank == 0)
  {
    size_t mylocvol = grid.localVolume();
    size_t totvol = grid.totalVolume();
    cout << "I am node: " << myrank << ", my local volume is: " << mylocvol << ", total volume is: " << totvol << endl;
    
    size_t const *sizes = grid.sizes();
    cout << "Local lattice size is: " << sizes[0] << " " << sizes[1] << " " << sizes[2] << " "  << sizes[3] << endl;
    
    size_t const *dims = grid.dims();
    cout << "Node dimensions are: " << dims[0] << " " << dims[1] << " " << dims[2] << " "  << dims[3] << endl;
    
    size_t const *surfaces = grid.surfaces();
    cout << "Surfaces are: " << surfaces[0] << " " << surfaces[1] << " " << surfaces[2] << " "  << surfaces[3] << endl;
    
    size_t const *dimSizes = grid.dimSizes();
    cout << "Dimsizes are: " << dimSizes[0] << " " << dimSizes[1] << " " << dimSizes[2] << " "  << dimSizes[3] << endl;
    
    cout << "The communication buffer volume is: " << grid.bufferVolume() << endl;
    cout << "The largest contiguous block is: " << grid.contiguousBlock() << endl;
  }
  Core::Field<QCD::Gauge, 24, 48> field(grid);
  field.readFromFile("../test/conf.save");  
  
  MPI::Finalize();
  return 0;
}


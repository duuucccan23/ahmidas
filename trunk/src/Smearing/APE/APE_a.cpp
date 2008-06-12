#include "APE.ih"

// Performs an APE smearing step.

void Smearing::APE::smear(Fields::GaugeField &field)
{
  // Create memory for storing the links
  // NOTE This is a busload of extra memory we need here.
  // More memory efficient implementations are available,

  Fields::GaugeField shifter(field);
  
  // We create buffers as we need them...
  Buffer< SU3::Matrix > stapleXplusY(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  Buffer< SU3::Matrix > stapleXplusZ(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  Buffer< SU3::Matrix > stapleXminZ(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  Buffer< SU3::Matrix > stapleYplusX(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  Buffer< SU3::Matrix > stapleZplusY(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  Buffer< SU3::Matrix > stapleZplusX(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  Buffer< SU3::Matrix > stapleYminZ(GaugeField.grid(), d_alpha * SU3::Matrix::identity());
  
  // We start from (0, 0, 0), where we will end as well...
    stapleYplusX.rightMultiply(shifter.component(idx_X)); // 1
    stapleZplusY.rightMultiply(shifter.component(idx_Y)); // 1
    stapleXplusZ.rightMultiply(shifter.component(idx_Z)); // 1
    stapleZplusX.rightMultiply(shifter.component(idx_X)); // 1

  shifter.shift(idx_X, dir_UP);   // ( 1, 0, 0)   1
  
    stapleXplusY.leftMultiply(shifter.component(idx_Y)); // 1
    stapleYplusX.rightMultiply(shifter.component(idx_Y)); // 2 
    stapleZplusX.rightMultiply(shifter.component(idx_Z)); // 2
  
  shifter.shift(idx_Z, dir_DOWN); // ( 1, 0,-1)   2
  
    stapleXminZ.leftMultiply(shifter.component(idx_Z)); // 1
  
  shifter.shift(idx_X, dir_DOWN); // ( 0, 0,-1)   3
  
    stapleXminZ.leftMultiply(shifter.component(idx_X)); // 2
    stapleXminZ.leftMultiply(shifter.component(idx_Z).dagger()); // 3
    stapleYminZ.rightMultiply(shifter.component(idx_Z).dagger()); // 1       
    stapleYminZ.rightMultiply(shifter.component(idx_Y)); // 2
    
    // Finished a staple! Add and reassign the memory
    std_foreach(field.component(idx_X).begin(), field.component(idx_X).end(), 
                Fields::add(stapleXminZ));   
    Buffer< SU3::Matrix > stapleYplusZ(stapleXminZ, d_alpha * SU3::Matrix::identity());
    
  shifter.shift(idx_Y, dir_UP);   // ( 0, 1,-1)   4
  
    stapleYminZ.rightMultiply(fieldComponent(idx_Z)); // 3
    
    // Finished a staple! Add and reassign the memory
    std_foreach(field.component(idx_Y).begin(), field.component(idx_Y).end(), 
                Fields::add(stapleYminZ));   
    Buffer< SU3::Matrix > stapleZminX(stapleYminZ, d_alpha * SU3::Matrix::identity());
    
  shifter.shift(idx_Z, dir_UP);    // ( 0, 1, 0)   5
  
    stapleXplusY.leftMultiply(shifter.component(idx_X));
    stapleYplusX.rightMultiply(shifter.component(idx_X).dagger());
    stapleYplusZ.leftMultiply(shifter.component(idx_Z).dagger());
    stapleZplusy.rightMultiply(shifter.component(idx_z));
    
    std_foreach(field.component(idx_Y).begin(), field.component(idx_Y).end(), 
                Fields::add(stapleYplusX));
    Buffer< SU3::Matrix > stapleYminX(stapleYplusX, d_alpha * SU3::Matrix::identity());
  
  shifter.shift(idx_X, dir_DOWN); // (-1, 1, 0)   6
  
    stapleYminX.leftMultiply(shifter.component(idx_X));
    
  shifter.shift(idx_Y, dir_DOWN); // (-1, 0, 0)   7
    
    stapleYminX.leftMultiply(shifter.component(idx_Y));
    stapleYminX.leftMultiply(shifter.component(idx_X).dagger());
    stapleZminX.rightMultiply(shifter.component(idx_X).dagger());
    stapleZminX.rightMultiply(shifter.component(idx_Z));
    
    std_foreach(field.component(idx_Y).begin(), field.component(idx_Y).end(), 
                Fields::add(stapleYminX));
    Buffer< SU3::Matrix > stapleXminY(stapleYminX, d_alpha * SU3::Matrix::identity());
  
  shifter.shift(idx_Z, dir_UP);   // (-1, 0, 1)   8
  
    stapleZminX.rightMultiply(shifter.component(idx_X));
    
    std_foreach(field.component(idx_Z).begin(), field.component(idx_Z).end(), 
                Fields::add(stapleZminX));
    // NOTE We're finished with this buffer - we may be able to compress things further.
      
  shifter.shift(idx_X, dir_UP);   // ( 0, 0, 1)   9
  
    stapleZminY.leftMultiply(shifter.component(idx_Y).dagger());
    stapleXplusZ.rightMultiply(shifter.component(idx_X)); 
    stapleYplusZ.leftMultiply(shifter.component(idx_Y));
    stapleZplusX.rightMultiply(shifter.component(idx_X).dagger());
    
    std_foreach(field.component(idx_Z).begin(), field.component(idx_Z).end(), 
                Fields::add(stapleZplusX));
    Buffer< SU3::Matrix > stapleZminY(stapleZplusX, d_alpha * SU3::Matrix::identity());
    std_foreach(field.component(idx_Z).begin(), field.component(idx_Z).end(), 
                Fields::add(stapleZplusY));
    Buffer< SU3::Matrix > stapleZminY(stapleZplusY, d_alpha * SU3::Matrix::identity());    
   
  shifter.shift(idx_Y, dir_DOWN); // ( 0,-1, 1)  10
  
    stapleZminY.leftMultiply(shifter.component(idx_Y));
    
  shifter.shift(idx_Z, dir_DOWN); // ( 0,-1, 0)  11
  
    stapleXminY.rightMultiply(shifter.component(idx_Y).dagger());
    stapleXminY.rightMultiply(shifter.component(idx_X));        
    stapleZminY.leftMultiply(shifter.component(idx_Z));
    stapleZminY.leftMultiply(shifter.component(idx_Y).dagger());    
    
    std_foreach(field.component(idx_Z).begin(), field.component(idx_Z).end(), 
                Fields::add(stapleZminY));
  
  shifter.shift(idx_X, dir_UP);   // ( 1,-1, 0)  12
  
    stapleXminY.rightMultiply(shifter.component(idx_Y));
    
    std_foreach(field.component(idx_X).begin(), field.component(idx_X).end(), 
                Fields::add(stapleXminY));
  
  shifter.shift(idx_Y, dir_UP);   // ( 1, 0, 0)  13
  
    stapleXplusZ.rightMultiply(shifter.component(idx_Z).dagger());
    
    std_foreach(field.component(idx_X).begin(), field.component(idx_X).end(), 
                Fields::add(stapleXminY));    
  
  shifter.shift(idx_X, dir_DOWN); // ( 0, 0, 0)  14   
    stapleXplusY.leftMultiply(shifter.component(idx_Y));
    stapleYplusZ.leftMultiply(shifter.component(idx_Z));
    
    std_foreach(field.component(idx_X).begin(), field.component(idx_X).end(), 
                Fields::add(stapleXplusY));    
    std_foreach(field.component(idx_Y).begin(), field.component(idx_Y).end(), 
                Fields::add(stapleYplusZ));    

  // That should be it! We made a full round and summed 12 staples over the field in total.
  // All that remains is making sure we get an SU(3) matrix back.
  field.reunitarize();
}

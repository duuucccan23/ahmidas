#include "Matrix.ih"

/***********************************************************
 *   The polar decomposition is implemented 
 *   following Fast polar decomposition 
 *   N.J. Higham /R.S. Schreiber (1998)
 * *********************************************************/

void SU3::Matrix::reunitarize() 
{ double check=10;
  int j=0;
  do
  {
    Matrix X_j = *this;
    Matrix Xbuff = X_j.inverse();
    double gamma = sqrt(Xbuff.norm()/Xbuff.norm());
    operator*=(0.5 * gamma);  operator+=((0.5/gamma) * Xbuff.dagger());
    X_j -= *this;
    check = sqrt(X_j.norm());  j++;
  }
  while(check > 1e-15);//std::cout << "j" << j << std::endl;
}





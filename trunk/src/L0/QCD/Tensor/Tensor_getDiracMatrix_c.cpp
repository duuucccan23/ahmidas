#include "Tensor.ih"

namespace QCD
{
  void getSpinDilutedDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, QCD::Tensor const &B)
  {
      /*  For a white source propagator we can use a faster constructor
      since we know which entries of the tensors are zero.
      This constructor is in particular used for of Dirac-diluted
      stochastic propagators, and therefore important in appication
      of the one-end-trick.
    */

    Dirac::Matrix xrA;
    Dirac::Matrix xgA;
    Dirac::Matrix xbA;
    Dirac::Matrix xrB;
    Dirac::Matrix xgB;
    Dirac::Matrix xbB;
    
    A.getDiracMatrix(xrA, Base::col_RED,   Base::col_RED);
    A.getDiracMatrix(xgA, Base::col_RED,   Base::col_GREEN);
    A.getDiracMatrix(xbA, Base::col_RED,   Base::col_BLUE);
    
    B.getDiracMatrix(xrB, Base::col_RED,   Base::col_RED);
    B.getDiracMatrix(xgB, Base::col_RED,   Base::col_GREEN);
    B.getDiracMatrix(xbB, Base::col_RED,   Base::col_BLUE);
    
    xrA *= xrB;
    xgA *= xgB;
    xbA *= xbB;
    
    xrA += xgA;
    xrA += xbA;
    std::copy(&(xrA[0]), &(xrA[0]) + 16, &(dMatrix[0]));
  }  
}

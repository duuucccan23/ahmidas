#include "Tensor.ih"

namespace QCD
{
  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, QCD::Tensor const &B, bool const colourDilutedSource)
  {
      /*  For a white source propagator we can use a faster constructor
      since we know which entries of the tensors are zero.
      This constructor is in particular used for of Dirac-diluted
      stochastic propagators, and therefore important in appication
      of the one-end-trick.
    */

    if (!colourDilutedSource)
    {
      // note that "RED" as the "source" index actually addresses the only non-zero entries
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
      B.getDiracMatrix(xgB, Base::col_GREEN, Base::col_RED);
      B.getDiracMatrix(xbB, Base::col_BLUE,  Base::col_RED);

      // note that second QCD::Tensor is supposed to appear as a
      // hermitian conjugate one, so source and sink indices are swapped
      xrA *= xrB;
      xgA *= xgB;
      xbA *= xbB;

      xrA += xgA;
      xrA += xbA;
      std::copy(&(xrA[0]), &(xrA[0]) + 16, &(dMatrix[0]));
    }
    else if (colourDilutedSource)
      getDiracMatrix(dMatrix, A, B);
  }
}

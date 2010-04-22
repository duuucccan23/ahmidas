#include "Tensor.ih"

namespace QCD
{

  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B)
  {
    Dirac::Matrix rrA;
    Dirac::Matrix rgA;
    Dirac::Matrix rbA;
    Dirac::Matrix grA;
    Dirac::Matrix gbA;
    Dirac::Matrix ggA;
    Dirac::Matrix brA;
    Dirac::Matrix bbA;
    Dirac::Matrix bgA;

    Dirac::Matrix rrB;
    Dirac::Matrix rgB;
    Dirac::Matrix rbB;
    Dirac::Matrix grB;
    Dirac::Matrix gbB;
    Dirac::Matrix ggB;
    Dirac::Matrix brB;
    Dirac::Matrix bbB;
    Dirac::Matrix bgB;

    A.getDiracMatrix(rrA, Base::col_RED,   Base::col_RED);
    A.getDiracMatrix(rgA, Base::col_RED,   Base::col_GREEN);
    A.getDiracMatrix(rbA, Base::col_RED,   Base::col_BLUE);
    A.getDiracMatrix(grA, Base::col_GREEN, Base::col_RED);
    A.getDiracMatrix(gbA, Base::col_GREEN, Base::col_BLUE);
    A.getDiracMatrix(ggA, Base::col_GREEN, Base::col_GREEN);
    A.getDiracMatrix(brA, Base::col_BLUE,  Base::col_RED);
    A.getDiracMatrix(bbA, Base::col_BLUE,  Base::col_BLUE);
    A.getDiracMatrix(bgA, Base::col_BLUE,  Base::col_GREEN);

    B.getDiracMatrix(rrB, Base::col_RED,   Base::col_RED);
    B.getDiracMatrix(rgB, Base::col_GREEN, Base::col_RED);
    B.getDiracMatrix(rbB, Base::col_BLUE,  Base::col_RED);
    B.getDiracMatrix(grB, Base::col_RED,   Base::col_GREEN);
    B.getDiracMatrix(gbB, Base::col_BLUE,  Base::col_GREEN);
    B.getDiracMatrix(ggB, Base::col_GREEN, Base::col_GREEN);
    B.getDiracMatrix(brB, Base::col_RED,   Base::col_BLUE);
    B.getDiracMatrix(bbB, Base::col_BLUE,  Base::col_BLUE);
    B.getDiracMatrix(bgB, Base::col_GREEN, Base::col_BLUE);

    rrA *= rrB;
    rgA *= rgB;
    rbA *= rbB;
    grA *= grB;
    gbA *= gbB;
    ggA *= ggB;
    brA *= brB;
    bbA *= bbB;
    bgA *= bgB;

    rrA += rgA;
    rrA += rbA;
    rrA += grA;
    rrA += gbA;
    rrA += ggA;
    rrA += brA;
    rrA += bbA;
    rrA += bgA;

    std::copy(&(rrA[0]), &(rrA[0]) + 16, &(dMatrix[0]));
  }
}

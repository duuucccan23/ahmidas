#include "Tensor.ih"

namespace QCD
{
  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol)
  {
  /*
      this is the underlying routine for a baryon contraction
  */
    // pre-contraction using interpolating field
    Dirac::Matrix result(std::complex< double >(0, 0));

    Tensor B_(B);

    switch(iPol)
    {
      case Base::bar_PROTON:
      {
        B_.left_multiply_proton();
        B_.right_multiply_proton();
        B_.transposeFull();
        break;
      }
      case Base::bar_PROTON_VAR:
      {
        Dirac::Gamma< 24 > g2g0;
        std::complex< double > const I(0.0, 1.0);
        B_*=g2g0;
        B_*= I*I;
        B_ = g2g0 * B_;
        B_.transposeFull();
        break;
      }
      default:
      std::cerr << "unknown interpolating field in Matrix::Matrix(3 x QCD::Tensor, Base::BaryonInterpolatingField)!"
                << std::endl;
      std::cerr << "Aborting..." << std::endl;
      exit(1);
    }

    // For each of the nine colour combinations for QCD::Tensor A there are only
    // four combinations of colour for the other tensors B and C
    // that yield a non-zero epsilon.
    // Those are, in terms of the Matrixs:
    //  +(rrA, ggB, bbC), -(rrA, gbB, bgC), -(rrA, bgB, gbC), +(rrA, bbB, ggC)
    //  -(rgA, grB, bbC), +(rgA, gbB, brC), +(rgA, brB, gbC), -(rgA, bbB, grC)
    //  +(rbA, grB, bgC), -(rbA, ggB, brC), -(rbA, brB, ggC), +(rbA, bgB, grC)
    //  ---
    //  +(grA, bgB, rbC), -(grA, bbB, rgC), -(grA, rgB, bbC), +(grA, rbB, bgC)
    //  +(ggA, bbB, rrC), -(ggA, brB, rbC), -(ggA, rbB, brC), +(ggA, rrB, bbC)
    //  -(gbA, bgB, rrC), +(gbA, brB, rgC), +(gbA, rgB, brC), -(gbA, rrB, bgC)
    //  ---
    //  +(brA, rgB, gbC), -(brA, rbB, ggC), -(brA, ggB, rbC), +(brA, gbB, rgC)
    //  +(bgA, rbB, grC), -(bgA, rrB, gbC), -(bgA, gbB, rrC), +(bgA, grB, rbC)
    //  -(bbA, rgB, grC), +(bbA, rrB, ggC), +(bbA, ggB, rrC), -(bbA, grB, rgC)

    Dirac::Matrix rrA;
    Dirac::Matrix rgA;
    Dirac::Matrix rbA;
    Dirac::Matrix grA;
    Dirac::Matrix ggA;
    Dirac::Matrix gbA;
    Dirac::Matrix brA;
    Dirac::Matrix bgA;
    Dirac::Matrix bbA;

    Dirac::Matrix rrB;
    Dirac::Matrix rgB;
    Dirac::Matrix rbB;
    Dirac::Matrix grB;
    Dirac::Matrix ggB;
    Dirac::Matrix gbB;
    Dirac::Matrix brB;
    Dirac::Matrix bgB;
    Dirac::Matrix bbB;

    Dirac::Matrix rrC;
    Dirac::Matrix rgC;
    Dirac::Matrix rbC;
    Dirac::Matrix grC;
    Dirac::Matrix ggC;
    Dirac::Matrix gbC;
    Dirac::Matrix brC;
    Dirac::Matrix bgC;
    Dirac::Matrix bbC;

    A.getDiracMatrix(rrA, Base::col_RED,   Base::col_RED);
    A.getDiracMatrix(rgA, Base::col_RED,   Base::col_GREEN);
    A.getDiracMatrix(rbA, Base::col_RED,   Base::col_BLUE);
    A.getDiracMatrix(grA, Base::col_GREEN, Base::col_RED);
    A.getDiracMatrix(ggA, Base::col_GREEN, Base::col_GREEN);
    A.getDiracMatrix(gbA, Base::col_GREEN, Base::col_BLUE);
    A.getDiracMatrix(brA, Base::col_BLUE,  Base::col_RED);
    A.getDiracMatrix(bgA, Base::col_BLUE,  Base::col_GREEN);
    A.getDiracMatrix(bbA, Base::col_BLUE,  Base::col_BLUE);

    // the source and sink colours are interchanged due to transposing
    B_.getDiracMatrix(rrB, Base::col_RED,   Base::col_RED);
    B_.getDiracMatrix(rgB, Base::col_GREEN, Base::col_RED);
    B_.getDiracMatrix(rbB, Base::col_BLUE,  Base::col_RED);
    B_.getDiracMatrix(grB, Base::col_RED,   Base::col_GREEN);
    B_.getDiracMatrix(ggB, Base::col_GREEN, Base::col_GREEN);
    B_.getDiracMatrix(gbB, Base::col_BLUE,  Base::col_GREEN);
    B_.getDiracMatrix(brB, Base::col_RED,   Base::col_BLUE);
    B_.getDiracMatrix(bgB, Base::col_GREEN, Base::col_BLUE);
    B_.getDiracMatrix(bbB, Base::col_BLUE,  Base::col_BLUE);

    C.getDiracMatrix(rrC, Base::col_RED,   Base::col_RED);
    C.getDiracMatrix(rgC, Base::col_RED,   Base::col_GREEN);
    C.getDiracMatrix(rbC, Base::col_RED,   Base::col_BLUE);
    C.getDiracMatrix(grC, Base::col_GREEN, Base::col_RED);
    C.getDiracMatrix(ggC, Base::col_GREEN, Base::col_GREEN);
    C.getDiracMatrix(gbC, Base::col_GREEN, Base::col_BLUE);
    C.getDiracMatrix(brC, Base::col_BLUE,  Base::col_RED);
    C.getDiracMatrix(bgC, Base::col_BLUE,  Base::col_GREEN);
    C.getDiracMatrix(bbC, Base::col_BLUE,  Base::col_BLUE);

    result += rrC*(ggB*bbA);
    result -= rrC*(gbB*bgA);
    result -= rrC*(bgB*gbA);
    result += rrC*(bbB*ggA);
    //  ---
    result -= rgC*(grB*bbA);
    result += rgC*(gbB*brA);
    result += rgC*(brB*gbA);
    result -= rgC*(bbB*grA);
    //  ---
    result += rbC*(grB*bgA);
    result -= rbC*(ggB*brA);
    result -= rbC*(brB*ggA);
    result += rbC*(bgB*grA);
    //  --- ---  ---
    result += grC*(bgB*rbA);
    result -= grC*(bbB*rgA);
    result -= grC*(rgB*bbA);
    result += grC*(rbB*bgA);
    //  ---
    result += ggC*(bbB*rrA);
    result -= ggC*(brB*rbA);
    result -= ggC*(rbB*brA);
    result += ggC*(rrB*bbA);
    //  ---
    result -= gbC*(bgB*rrA);
    result += gbC*(brB*rgA);
    result += gbC*(rgB*brA);
    result -= gbC*(rrB*bgA);
    //  --- ---  ---
    result += brC*(rgB*gbA);
    result -= brC*(rbB*ggA);
    result -= brC*(ggB*rbA);
    result += brC*(gbB*rgA);
    //  ---
    result += bgC*(rbB*grA);
    result -= bgC*(rrB*gbA);
    result -= bgC*(gbB*rrA);
    result += bgC*(grB*rbA);
    //  ---
    result -= bbC*(rgB*grA);
    result += bbC*(rrB*ggA);
    result += bbC*(ggB*rrA);
    result -= bbC*(grB*rgA);

    switch(iPol)
    {
      case Base::bar_PROTON_VAR:
//       {
//         Dirac::Gamma< 5 > gamma5;
//         result = gamma5 * result;
//         result *= gamma5;
//         
//       }
      case Base::bar_PROTON:

        result += rrA*((ggB*bbC).trace());
        result -= rrA*((gbB*bgC).trace());
        result -= rrA*((bgB*gbC).trace());
        result += rrA*((bbB*ggC).trace());
        //  ---
        result -= rgA*((grB*bbC).trace());
        result += rgA*((gbB*brC).trace());
        result += rgA*((brB*gbC).trace());
        result -= rgA*((bbB*grC).trace());
        //  ---
        result += rbA*((grB*bgC).trace());
        result -= rbA*((ggB*brC).trace());
        result -= rbA*((brB*ggC).trace());
        result += rbA*((bgB*grC).trace());
        //  --- --- ---
        result += grA*((bgB*rbC).trace());
        result -= grA*((bbB*rgC).trace());
        result -= grA*((rgB*bbC).trace());
        result += grA*((rbB*bgC).trace());
        //  ---
        result += ggA*((bbB*rrC).trace());
        result -= ggA*((brB*rbC).trace());
        result -= ggA*((rbB*brC).trace());
        result += ggA*((rrB*bbC).trace());
        //  ---
        result -= gbA*((bgB*rrC).trace());
        result += gbA*((brB*rgC).trace());
        result += gbA*((rgB*brC).trace());
        result -= gbA*((rrB*bgC).trace());
        //  --- --- ---
        result += brA*((rgB*gbC).trace());
        result -= brA*((rbB*ggC).trace());
        result -= brA*((ggB*rbC).trace());
        result += brA*((gbB*rgC).trace());
        //  ---
        result += bgA*((rbB*grC).trace());
        result -= bgA*((rrB*gbC).trace());
        result -= bgA*((gbB*rrC).trace());
        result += bgA*((grB*rbC).trace());
        //  ---
        result -= bbA*((rgB*grC).trace());
        result += bbA*((rrB*ggC).trace());
        result += bbA*((ggB*rrC).trace());
        result -= bbA*((grB*rgC).trace());

      break;
    }
//     std::cout << result << std::endl;
//     std::cout << "trace: " << result.trace() << std::endl;
    std::copy(&(result[0]), &(result[0]) + 16, &(dMatrix[0]));
  }
}

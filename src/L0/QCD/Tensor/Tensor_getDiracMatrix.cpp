#include "Tensor.ih"

namespace QCD
{

  void Tensor::getDiracMatrix(Dirac::Matrix &dMatrix, Base::ColourIndex const colour_src, Base::ColourIndex const colour_snk) const
  {
    size_t index = 12*colour_src + colour_snk;
    dMatrix[ 0] = d_data[index     ];
    dMatrix[ 1] = d_data[index +   3];
    dMatrix[ 2] = d_data[index +   6];
    dMatrix[ 3] = d_data[index +   9];
    dMatrix[ 4] = d_data[index +  36];
    dMatrix[ 5] = d_data[index +  39];
    dMatrix[ 6] = d_data[index +  42];
    dMatrix[ 7] = d_data[index +  45];
    dMatrix[ 8] = d_data[index +  72];
    dMatrix[ 9] = d_data[index +  75];
    dMatrix[10] = d_data[index +  78];
    dMatrix[11] = d_data[index +  81];
    dMatrix[12] = d_data[index + 108];
    dMatrix[13] = d_data[index + 111];
    dMatrix[14] = d_data[index + 114];
    dMatrix[15] = d_data[index + 117];
  }



  /* the following functions are used to perform contractions */

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

  /*  For a white source propagator we can use a faster constructor
      since we know which entries of the tensors are zero.
      This constructor is in particular used for of Dirac-diluted
      stochastic propagators, and therefore important in appication
      of the one-end-trick.
  */
  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, QCD::Tensor const &B, bool const colourDilutedSource)
  {
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


  /*
      this is the underlying routine for a baryon contraction
  */
  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol)
  {
    // pre-contraction using interpolating field
    Dirac::Matrix result(std::complex< double >(0, 0));

    Tensor B_(B);

    switch(iPol)
    {
      case Base::bar_PROTON:
        B_.left_multiply_proton();
        B_.right_multiply_proton();
        B_.transposeFull();
      break;
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

  // alternative version, gives the same result as function above if Tensors A and C are equal
  void getDiracMatrix_alternative(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol)
  {
    // pre-contraction using interpolating field
    Dirac::Matrix result(std::complex< double >(0, 0));

    Tensor B_(B);

    switch(iPol)
    {
      case Base::bar_PROTON:
        B_.left_multiply_proton();
        B_.right_multiply_proton();
        B_.transposeFull();
      break;
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

    result += rrA*(ggB*bbC);
    result -= rrA*(gbB*bgC);
    result -= rrA*(bgB*gbC);
    result += rrA*(bbB*ggC);
    //  ---                
    result -= rgA*(grB*bbC);
    result += rgA*(gbB*brC);
    result += rgA*(brB*gbC);
    result -= rgA*(bbB*grC);
    //  ---                
    result += rbA*(grB*bgC);
    result -= rbA*(ggB*brC);
    result -= rbA*(brB*ggC);
    result += rbA*(bgB*grC);
    //  --- --- ---        
    result += grA*(bgB*rbC);
    result -= grA*(bbB*rgC);
    result -= grA*(rgB*bbC);
    result += grA*(rbB*bgC);
    //  ---                
    result += ggA*(bbB*rrC);
    result -= ggA*(brB*rbC);
    result -= ggA*(rbB*brC);
    result += ggA*(rrB*bbC);
    //  ---                
    result -= gbA*(bgB*rrC);
    result += gbA*(brB*rgC);
    result += gbA*(rgB*brC);
    result -= gbA*(rrB*bgC);
    //  --- --- ---        
    result += brA*(rgB*gbC);
    result -= brA*(rbB*ggC);
    result -= brA*(ggB*rbC);
    result += brA*(gbB*rgC);
    //  ---                
    result += bgA*(rbB*grC);
    result -= bgA*(rrB*gbC);
    result -= bgA*(gbB*rrC);
    result += bgA*(grB*rbC);
    //  ---                
    result -= bbA*(rgB*grC);
    result += bbA*(rrB*ggC);
    result += bbA*(ggB*rrC);
    result -= bbA*(grB*rgC);

    switch(iPol)
    {
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

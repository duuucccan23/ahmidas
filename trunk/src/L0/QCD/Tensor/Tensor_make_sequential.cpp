#include "Tensor.ih"

namespace QCD
{

  // this has to be reviewed concerning efficiency
  // this routine contracts the color indices according to the antisymmetric epsilon tensor
  // i.e. the result with color index c (source or sink index, respectively)
  // contains the combinations epsilon_abc - epsilon_bac
  // where the first index is an index of Tensor A and the second of Tensor B
  void Tensor::make_sequential(Tensor const &A, Tensor const &B)
  {
    reducedTensor rrA(A, Base::col_RED,   Base::col_RED);
    reducedTensor rgA(A, Base::col_RED,   Base::col_GREEN);
    reducedTensor rbA(A, Base::col_RED,   Base::col_BLUE);
    reducedTensor grA(A, Base::col_GREEN, Base::col_RED);
    reducedTensor ggA(A, Base::col_GREEN, Base::col_GREEN);
    reducedTensor gbA(A, Base::col_GREEN, Base::col_BLUE);
    reducedTensor brA(A, Base::col_BLUE,  Base::col_RED);
    reducedTensor bgA(A, Base::col_BLUE,  Base::col_GREEN);
    reducedTensor bbA(A, Base::col_BLUE,  Base::col_BLUE);

    reducedTensor rrB(B, Base::col_RED,   Base::col_RED);
    reducedTensor rgB(B, Base::col_RED,   Base::col_GREEN);
    reducedTensor rbB(B, Base::col_RED,   Base::col_BLUE);
    reducedTensor grB(B, Base::col_GREEN, Base::col_RED);
    reducedTensor gbB(B, Base::col_GREEN, Base::col_BLUE);
    reducedTensor ggB(B, Base::col_GREEN, Base::col_GREEN);
    reducedTensor brB(B, Base::col_BLUE,  Base::col_RED);
    reducedTensor bbB(B, Base::col_BLUE,  Base::col_BLUE);
    reducedTensor bgB(B, Base::col_BLUE,  Base::col_GREEN);

    // Non-zero combinations
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


    reducedTensor rrC  = ggA.elementwise_product(bbB);
    rrC               -= gbA.elementwise_product(bgB);
    rrC               -= bgA.elementwise_product(gbB);
    rrC               += bbA.elementwise_product(ggB);

    reducedTensor rgC  = gbA.elementwise_product(brB);
    rgC               -= grA.elementwise_product(bbB);
    rgC               -= bbA.elementwise_product(grB);
    rgC               += brA.elementwise_product(gbB);

    reducedTensor rbC  = grA.elementwise_product(bgB);
    rbC               -= ggA.elementwise_product(brB);
    rbC               -= brA.elementwise_product(ggB);
    rbC               += bgA.elementwise_product(grB);

    /* ----------------------- */

    reducedTensor grC  = rbA.elementwise_product(bgB);
    grC               -= rgA.elementwise_product(bbB);
    grC               -= bbA.elementwise_product(rgB);
    grC               += bgA.elementwise_product(rbB);

    reducedTensor ggC  = rrA.elementwise_product(bbB);
    ggC               -= rbA.elementwise_product(brB);
    ggC               -= brA.elementwise_product(rbB);
    ggC               += bbA.elementwise_product(rrB);

    reducedTensor gbC  = rgA.elementwise_product(brB);
    gbC               -= rrA.elementwise_product(bgB);
    gbC               -= bgA.elementwise_product(rrB);
    gbC               += brA.elementwise_product(rgB);

    /* ----------------------- */

    reducedTensor brC  = rgA.elementwise_product(gbB);
    brC               -= rbA.elementwise_product(ggB);
    brC               -= ggA.elementwise_product(rbB);
    brC               += gbA.elementwise_product(rgB);

    reducedTensor bgC  = rbA.elementwise_product(grB);
    bgC               -= rrA.elementwise_product(gbB);
    bgC               -= gbA.elementwise_product(rrB);
    bgC               += grA.elementwise_product(rbB);

    reducedTensor bbC  = rrA.elementwise_product(ggB);
    bbC               -= rgA.elementwise_product(grB);
    bbC               -= grA.elementwise_product(rgB);
    bbC               += ggA.elementwise_product(rrB);

    /* --------------------------------------------- */

    size_t index = 12*Base::col_RED  + Base::col_RED;
    d_data[index      ] =  rrC.d_data[ 0];
    d_data[index +   3] =  rrC.d_data[ 1];
    d_data[index +   6] =  rrC.d_data[ 2];
    d_data[index +   9] =  rrC.d_data[ 3];
    d_data[index +  36] =  rrC.d_data[ 4];
    d_data[index +  39] =  rrC.d_data[ 5];
    d_data[index +  42] =  rrC.d_data[ 6];
    d_data[index +  45] =  rrC.d_data[ 7];
    d_data[index +  72] =  rrC.d_data[ 8];
    d_data[index +  75] =  rrC.d_data[ 9];
    d_data[index +  78] =  rrC.d_data[10];
    d_data[index +  81] =  rrC.d_data[11];
    d_data[index + 108] =  rrC.d_data[12];
    d_data[index + 111] =  rrC.d_data[13];
    d_data[index + 114] =  rrC.d_data[14];
    d_data[index + 117] =  rrC.d_data[15];

    index = 12*Base::col_RED  + Base::col_GREEN;
    d_data[index      ] =  rgC.d_data[ 0];
    d_data[index +   3] =  rgC.d_data[ 1];
    d_data[index +   6] =  rgC.d_data[ 2];
    d_data[index +   9] =  rgC.d_data[ 3];
    d_data[index +  36] =  rgC.d_data[ 4];
    d_data[index +  39] =  rgC.d_data[ 5];
    d_data[index +  42] =  rgC.d_data[ 6];
    d_data[index +  45] =  rgC.d_data[ 7];
    d_data[index +  72] =  rgC.d_data[ 8];
    d_data[index +  75] =  rgC.d_data[ 9];
    d_data[index +  78] =  rgC.d_data[10];
    d_data[index +  81] =  rgC.d_data[11];
    d_data[index + 108] =  rgC.d_data[12];
    d_data[index + 111] =  rgC.d_data[13];
    d_data[index + 114] =  rgC.d_data[14];
    d_data[index + 117] =  rgC.d_data[15];

    index = 12*Base::col_RED  + Base::col_BLUE;
    d_data[index      ] =  rbC.d_data[ 0];
    d_data[index +   3] =  rbC.d_data[ 1];
    d_data[index +   6] =  rbC.d_data[ 2];
    d_data[index +   9] =  rbC.d_data[ 3];
    d_data[index +  36] =  rbC.d_data[ 4];
    d_data[index +  39] =  rbC.d_data[ 5];
    d_data[index +  42] =  rbC.d_data[ 6];
    d_data[index +  45] =  rbC.d_data[ 7];
    d_data[index +  72] =  rbC.d_data[ 8];
    d_data[index +  75] =  rbC.d_data[ 9];
    d_data[index +  78] =  rbC.d_data[10];
    d_data[index +  81] =  rbC.d_data[11];
    d_data[index + 108] =  rbC.d_data[12];
    d_data[index + 111] =  rbC.d_data[13];
    d_data[index + 114] =  rbC.d_data[14];
    d_data[index + 117] =  rbC.d_data[15];

    /* --------------------------------------------- */

    index = 12*Base::col_GREEN  + Base::col_RED;
    d_data[index      ] =  grC.d_data[ 0];
    d_data[index +   3] =  grC.d_data[ 1];
    d_data[index +   6] =  grC.d_data[ 2];
    d_data[index +   9] =  grC.d_data[ 3];
    d_data[index +  36] =  grC.d_data[ 4];
    d_data[index +  39] =  grC.d_data[ 5];
    d_data[index +  42] =  grC.d_data[ 6];
    d_data[index +  45] =  grC.d_data[ 7];
    d_data[index +  72] =  grC.d_data[ 8];
    d_data[index +  75] =  grC.d_data[ 9];
    d_data[index +  78] =  grC.d_data[10];
    d_data[index +  81] =  grC.d_data[11];
    d_data[index + 108] =  grC.d_data[12];
    d_data[index + 111] =  grC.d_data[13];
    d_data[index + 114] =  grC.d_data[14];
    d_data[index + 117] =  grC.d_data[15];

    index = 12*Base::col_GREEN  + Base::col_GREEN;
    d_data[index      ] =  ggC.d_data[ 0];
    d_data[index +   3] =  ggC.d_data[ 1];
    d_data[index +   6] =  ggC.d_data[ 2];
    d_data[index +   9] =  ggC.d_data[ 3];
    d_data[index +  36] =  ggC.d_data[ 4];
    d_data[index +  39] =  ggC.d_data[ 5];
    d_data[index +  42] =  ggC.d_data[ 6];
    d_data[index +  45] =  ggC.d_data[ 7];
    d_data[index +  72] =  ggC.d_data[ 8];
    d_data[index +  75] =  ggC.d_data[ 9];
    d_data[index +  78] =  ggC.d_data[10];
    d_data[index +  81] =  ggC.d_data[11];
    d_data[index + 108] =  ggC.d_data[12];
    d_data[index + 111] =  ggC.d_data[13];
    d_data[index + 114] =  ggC.d_data[14];
    d_data[index + 117] =  ggC.d_data[15];

    index = 12*Base::col_GREEN  + Base::col_BLUE;
    d_data[index      ] =  gbC.d_data[ 0];
    d_data[index +   3] =  gbC.d_data[ 1];
    d_data[index +   6] =  gbC.d_data[ 2];
    d_data[index +   9] =  gbC.d_data[ 3];
    d_data[index +  36] =  gbC.d_data[ 4];
    d_data[index +  39] =  gbC.d_data[ 5];
    d_data[index +  42] =  gbC.d_data[ 6];
    d_data[index +  45] =  gbC.d_data[ 7];
    d_data[index +  72] =  gbC.d_data[ 8];
    d_data[index +  75] =  gbC.d_data[ 9];
    d_data[index +  78] =  gbC.d_data[10];
    d_data[index +  81] =  gbC.d_data[11];
    d_data[index + 108] =  gbC.d_data[12];
    d_data[index + 111] =  gbC.d_data[13];
    d_data[index + 114] =  gbC.d_data[14];
    d_data[index + 117] =  gbC.d_data[15];

    /* --------------------------------------------- */

    index = 12*Base::col_BLUE  + Base::col_RED;
    d_data[index      ] =  brC.d_data[ 0];
    d_data[index +   3] =  brC.d_data[ 1];
    d_data[index +   6] =  brC.d_data[ 2];
    d_data[index +   9] =  brC.d_data[ 3];
    d_data[index +  36] =  brC.d_data[ 4];
    d_data[index +  39] =  brC.d_data[ 5];
    d_data[index +  42] =  brC.d_data[ 6];
    d_data[index +  45] =  brC.d_data[ 7];
    d_data[index +  72] =  brC.d_data[ 8];
    d_data[index +  75] =  brC.d_data[ 9];
    d_data[index +  78] =  brC.d_data[10];
    d_data[index +  81] =  brC.d_data[11];
    d_data[index + 108] =  brC.d_data[12];
    d_data[index + 111] =  brC.d_data[13];
    d_data[index + 114] =  brC.d_data[14];
    d_data[index + 117] =  brC.d_data[15];

    index = 12*Base::col_BLUE  + Base::col_GREEN;
    d_data[index      ] =  bgC.d_data[ 0];
    d_data[index +   3] =  bgC.d_data[ 1];
    d_data[index +   6] =  bgC.d_data[ 2];
    d_data[index +   9] =  bgC.d_data[ 3];
    d_data[index +  36] =  bgC.d_data[ 4];
    d_data[index +  39] =  bgC.d_data[ 5];
    d_data[index +  42] =  bgC.d_data[ 6];
    d_data[index +  45] =  bgC.d_data[ 7];
    d_data[index +  72] =  bgC.d_data[ 8];
    d_data[index +  75] =  bgC.d_data[ 9];
    d_data[index +  78] =  bgC.d_data[10];
    d_data[index +  81] =  bgC.d_data[11];
    d_data[index + 108] =  bgC.d_data[12];
    d_data[index + 111] =  bgC.d_data[13];
    d_data[index + 114] =  bgC.d_data[14];
    d_data[index + 117] =  bgC.d_data[15];

    index = 12*Base::col_BLUE  + Base::col_BLUE;
    d_data[index      ] =  bbC.d_data[ 0];
    d_data[index +   3] =  bbC.d_data[ 1];
    d_data[index +   6] =  bbC.d_data[ 2];
    d_data[index +   9] =  bbC.d_data[ 3];
    d_data[index +  36] =  bbC.d_data[ 4];
    d_data[index +  39] =  bbC.d_data[ 5];
    d_data[index +  42] =  bbC.d_data[ 6];
    d_data[index +  45] =  bbC.d_data[ 7];
    d_data[index +  72] =  bbC.d_data[ 8];
    d_data[index +  75] =  bbC.d_data[ 9];
    d_data[index +  78] =  bbC.d_data[10];
    d_data[index +  81] =  bbC.d_data[11];
    d_data[index + 108] =  bbC.d_data[12];
    d_data[index + 111] =  bbC.d_data[13];
    d_data[index + 114] =  bbC.d_data[14];
    d_data[index + 117] =  bbC.d_data[15];

//     std::cout << *this << std::endl;
//     exit(1);
  }




  void make_sequential(Tensor result[16], Tensor const &A, Tensor const &B)
  {
    reducedTensor rrA(A, Base::col_RED,   Base::col_RED);
    reducedTensor rgA(A, Base::col_RED,   Base::col_GREEN);
    reducedTensor rbA(A, Base::col_RED,   Base::col_BLUE);
    reducedTensor grA(A, Base::col_GREEN, Base::col_RED);
    reducedTensor ggA(A, Base::col_GREEN, Base::col_GREEN);
    reducedTensor gbA(A, Base::col_GREEN, Base::col_BLUE);
    reducedTensor brA(A, Base::col_BLUE,  Base::col_RED);
    reducedTensor bgA(A, Base::col_BLUE,  Base::col_GREEN);
    reducedTensor bbA(A, Base::col_BLUE,  Base::col_BLUE);

    reducedTensor rrB(B, Base::col_RED,   Base::col_RED);
    reducedTensor rgB(B, Base::col_RED,   Base::col_GREEN);
    reducedTensor rbB(B, Base::col_RED,   Base::col_BLUE);
    reducedTensor grB(B, Base::col_GREEN, Base::col_RED);
    reducedTensor gbB(B, Base::col_GREEN, Base::col_BLUE);
    reducedTensor ggB(B, Base::col_GREEN, Base::col_GREEN);
    reducedTensor brB(B, Base::col_BLUE,  Base::col_RED);
    reducedTensor bbB(B, Base::col_BLUE,  Base::col_BLUE);
    reducedTensor bgB(B, Base::col_BLUE,  Base::col_GREEN);

    size_t const n_complex = 16*16;

    std::complex< double > tmp[n_complex];

    std::complex< double > zero(0, 0);

    //reducedTensor *tmpTensors = reinterpret_cast< reducedTensor* >(tmp);

    std::complex< double > *tmpTensors = tmp;

    std::complex< double > rrC[n_complex];
    std::fill_n(rrC, n_complex, zero);
    ggA.outer_product(bbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    gbA.outer_product(bgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bgA.outer_product(gbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bbA.outer_product(ggB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());


    std::complex< double > rgC[n_complex];
    std::fill_n(rgC, n_complex, zero);
    gbA.outer_product(brB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    grA.outer_product(bbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    bbA.outer_product(grB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    brA.outer_product(gbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());

    std::complex< double > rbC[n_complex];
    std::fill_n(rbC, n_complex, zero);
    grA.outer_product(bgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    ggA.outer_product(brB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    brA.outer_product(ggB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    bgA.outer_product(grB, tmpTensors);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > grC[n_complex];
    std::fill_n(grC, n_complex, zero);
    rbA.outer_product(bgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rgA.outer_product(bbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bbA.outer_product(rgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bgA.outer_product(rbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());

    std::complex< double > ggC[n_complex];
    std::fill_n(ggC, n_complex, zero);
    rrA.outer_product(bbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rbA.outer_product(brB, tmpTensors);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    brA.outer_product(rbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    bbA.outer_product(rrB, tmpTensors);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());

    std::complex< double > gbC[n_complex];
    std::fill_n(gbC, n_complex, zero);
    rgA.outer_product(brB, tmpTensors);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rrA.outer_product(bgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    bgA.outer_product(rrB, tmpTensors);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    brA.outer_product(rgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > brC[n_complex];
    std::fill_n(brC, n_complex, zero);
    rgA.outer_product(gbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rbA.outer_product(ggB, tmpTensors);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    ggA.outer_product(rbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    gbA.outer_product(rgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());

    std::complex< double > bgC[n_complex];
    std::fill_n(bgC, n_complex, zero);
    rbA.outer_product(grB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rrA.outer_product(gbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    gbA.outer_product(rrB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    grA.outer_product(rbB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());

    std::complex< double > bbC[n_complex];
    std::fill_n(bbC, n_complex, zero);
    rrA.outer_product(ggB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rgA.outer_product(grB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    grA.outer_product(rgB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    ggA.outer_product(rrB, tmpTensors);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());


//     for (size_t iDirac=0; iDirac<16*16; iDirac++)
//     {
//       if (iDirac % 16 == 0)
//         std::cout << std::endl;
//       std::cout << tmp[iDirac] << " ";
//     }
//     std::cout << std::endl;

    /* --------------------------------------------- */

    std::complex< double >* data_ptr = NULL;

    size_t index(-1);

    for (size_t iDirac=0; iDirac<16; iDirac++)
    {
      index = 12*Base::col_RED  + Base::col_RED;
      data_ptr = (result[iDirac]).d_data + index;

      size_t const offset = iDirac * 16;

      data_ptr[  0] =  rrC[offset ];
      data_ptr[  3] =  rrC[offset + 1];
      data_ptr[  6] =  rrC[offset + 2];
      data_ptr[  9] =  rrC[offset + 3];
      data_ptr[ 36] =  rrC[offset + 4];
      data_ptr[ 39] =  rrC[offset + 5];
      data_ptr[ 42] =  rrC[offset + 6];
      data_ptr[ 45] =  rrC[offset + 7];
      data_ptr[ 72] =  rrC[offset + 8];
      data_ptr[ 75] =  rrC[offset + 9];
      data_ptr[ 78] =  rrC[offset +10];
      data_ptr[ 81] =  rrC[offset +11];
      data_ptr[108] =  rrC[offset +12];
      data_ptr[111] =  rrC[offset +13];
      data_ptr[114] =  rrC[offset +14];
      data_ptr[117] =  rrC[offset +15];

      index = 12*Base::col_RED  + Base::col_GREEN;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  rgC[offset + 0];
      data_ptr[  3] =  rgC[offset + 1];
      data_ptr[  6] =  rgC[offset + 2];
      data_ptr[  9] =  rgC[offset + 3];
      data_ptr[ 36] =  rgC[offset + 4];
      data_ptr[ 39] =  rgC[offset + 5];
      data_ptr[ 42] =  rgC[offset + 6];
      data_ptr[ 45] =  rgC[offset + 7];
      data_ptr[ 72] =  rgC[offset + 8];
      data_ptr[ 75] =  rgC[offset + 9];
      data_ptr[ 78] =  rgC[offset +10];
      data_ptr[ 81] =  rgC[offset +11];
      data_ptr[108] =  rgC[offset +12];
      data_ptr[111] =  rgC[offset +13];
      data_ptr[114] =  rgC[offset +14];
      data_ptr[117] =  rgC[offset +15];

      index = 12*Base::col_RED  + Base::col_BLUE;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  rbC[offset + 0];
      data_ptr[  3] =  rbC[offset + 1];
      data_ptr[  6] =  rbC[offset + 2];
      data_ptr[  9] =  rbC[offset + 3];
      data_ptr[ 36] =  rbC[offset + 4];
      data_ptr[ 39] =  rbC[offset + 5];
      data_ptr[ 42] =  rbC[offset + 6];
      data_ptr[ 45] =  rbC[offset + 7];
      data_ptr[ 72] =  rbC[offset + 8];
      data_ptr[ 75] =  rbC[offset + 9];
      data_ptr[ 78] =  rbC[offset +10];
      data_ptr[ 81] =  rbC[offset +11];
      data_ptr[108] =  rbC[offset +12];
      data_ptr[111] =  rbC[offset +13];
      data_ptr[114] =  rbC[offset +14];
      data_ptr[117] =  rbC[offset +15];

      /* --------------------------------------------- */

      index = 12*Base::col_GREEN  + Base::col_RED;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  grC[offset + 0];
      data_ptr[  3] =  grC[offset + 1];
      data_ptr[  6] =  grC[offset + 2];
      data_ptr[  9] =  grC[offset + 3];
      data_ptr[ 36] =  grC[offset + 4];
      data_ptr[ 39] =  grC[offset + 5];
      data_ptr[ 42] =  grC[offset + 6];
      data_ptr[ 45] =  grC[offset + 7];
      data_ptr[ 72] =  grC[offset + 8];
      data_ptr[ 75] =  grC[offset + 9];
      data_ptr[ 78] =  grC[offset +10];
      data_ptr[ 81] =  grC[offset +11];
      data_ptr[108] =  grC[offset +12];
      data_ptr[111] =  grC[offset +13];
      data_ptr[114] =  grC[offset +14];
      data_ptr[117] =  grC[offset +15];

      index = 12*Base::col_GREEN  + Base::col_GREEN;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  ggC[offset + 0];
      data_ptr[  3] =  ggC[offset + 1];
      data_ptr[  6] =  ggC[offset + 2];
      data_ptr[  9] =  ggC[offset + 3];
      data_ptr[ 36] =  ggC[offset + 4];
      data_ptr[ 39] =  ggC[offset + 5];
      data_ptr[ 42] =  ggC[offset + 6];
      data_ptr[ 45] =  ggC[offset + 7];
      data_ptr[ 72] =  ggC[offset + 8];
      data_ptr[ 75] =  ggC[offset + 9];
      data_ptr[ 78] =  ggC[offset +10];
      data_ptr[ 81] =  ggC[offset +11];
      data_ptr[108] =  ggC[offset +12];
      data_ptr[111] =  ggC[offset +13];
      data_ptr[114] =  ggC[offset +14];
      data_ptr[117] =  ggC[offset +15];

      index = 12*Base::col_GREEN  + Base::col_BLUE;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  gbC[offset + 0];
      data_ptr[  3] =  gbC[offset + 1];
      data_ptr[  6] =  gbC[offset + 2];
      data_ptr[  9] =  gbC[offset + 3];
      data_ptr[ 36] =  gbC[offset + 4];
      data_ptr[ 39] =  gbC[offset + 5];
      data_ptr[ 42] =  gbC[offset + 6];
      data_ptr[ 45] =  gbC[offset + 7];
      data_ptr[ 72] =  gbC[offset + 8];
      data_ptr[ 75] =  gbC[offset + 9];
      data_ptr[ 78] =  gbC[offset +10];
      data_ptr[ 81] =  gbC[offset +11];
      data_ptr[108] =  gbC[offset +12];
      data_ptr[111] =  gbC[offset +13];
      data_ptr[114] =  gbC[offset +14];
      data_ptr[117] =  gbC[offset +15];

      /* --------------------------------------------- */

      index = 12*Base::col_BLUE  + Base::col_RED;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  brC[offset + 0];
      data_ptr[  3] =  brC[offset + 1];
      data_ptr[  6] =  brC[offset + 2];
      data_ptr[  9] =  brC[offset + 3];
      data_ptr[ 36] =  brC[offset + 4];
      data_ptr[ 39] =  brC[offset + 5];
      data_ptr[ 42] =  brC[offset + 6];
      data_ptr[ 45] =  brC[offset + 7];
      data_ptr[ 72] =  brC[offset + 8];
      data_ptr[ 75] =  brC[offset + 9];
      data_ptr[ 78] =  brC[offset +10];
      data_ptr[ 81] =  brC[offset +11];
      data_ptr[108] =  brC[offset +12];
      data_ptr[111] =  brC[offset +13];
      data_ptr[114] =  brC[offset +14];
      data_ptr[117] =  brC[offset +15];

      index = 12*Base::col_BLUE  + Base::col_GREEN;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  bgC[offset + 0];
      data_ptr[  3] =  bgC[offset + 1];
      data_ptr[  6] =  bgC[offset + 2];
      data_ptr[  9] =  bgC[offset + 3];
      data_ptr[ 36] =  bgC[offset + 4];
      data_ptr[ 39] =  bgC[offset + 5];
      data_ptr[ 42] =  bgC[offset + 6];
      data_ptr[ 45] =  bgC[offset + 7];
      data_ptr[ 72] =  bgC[offset + 8];
      data_ptr[ 75] =  bgC[offset + 9];
      data_ptr[ 78] =  bgC[offset +10];
      data_ptr[ 81] =  bgC[offset +11];
      data_ptr[108] =  bgC[offset +12];
      data_ptr[111] =  bgC[offset +13];
      data_ptr[114] =  bgC[offset +14];
      data_ptr[117] =  bgC[offset +15];

      index = 12*Base::col_BLUE  + Base::col_BLUE;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  bbC[offset + 0];
      data_ptr[  3] =  bbC[offset + 1];
      data_ptr[  6] =  bbC[offset + 2];
      data_ptr[  9] =  bbC[offset + 3];
      data_ptr[ 36] =  bbC[offset + 4];
      data_ptr[ 39] =  bbC[offset + 5];
      data_ptr[ 42] =  bbC[offset + 6];
      data_ptr[ 45] =  bbC[offset + 7];
      data_ptr[ 72] =  bbC[offset + 8];
      data_ptr[ 75] =  bbC[offset + 9];
      data_ptr[ 78] =  bbC[offset +10];
      data_ptr[ 81] =  bbC[offset +11];
      data_ptr[108] =  bbC[offset +12];
      data_ptr[111] =  bbC[offset +13];
      data_ptr[114] =  bbC[offset +14];
      data_ptr[117] =  bbC[offset +15];

      data_ptr = NULL;
    }

//     for (size_t iDirac=0; iDirac<16; iDirac++)
//     {
//       std::cout << "index: ";
//       std::cout.width(3);
//       std::cout << iDirac << std::endl;
//       std::cout << result[iDirac] << std::endl;
//     }
//     exit(1);
  }

//   getDiracTr(Tensor &result, Tensor const * const doubleDirac)
//   {
//     reducedTensor *tmp[9];
//     result = Tensor(tmp);
//   }


}

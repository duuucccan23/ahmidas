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


    Dirac::Matrix rrC  = ggA.elementwise_product(bbB);
    rrC               -= gbA.elementwise_product(bgB);
    rrC               -= bgA.elementwise_product(gbB);
    rrC               += bbA.elementwise_product(ggB);

    Dirac::Matrix rgC  = gbA.elementwise_product(brB);
    rgC               -= grA.elementwise_product(bbB);
    rgC               -= bbA.elementwise_product(grB);
    rgC               += brA.elementwise_product(gbB);

    Dirac::Matrix rbC  = grA.elementwise_product(bgB);
    rbC               -= ggA.elementwise_product(brB);
    rbC               -= brA.elementwise_product(ggB);
    rbC               += bgA.elementwise_product(grB);

    /* ----------------------- */

    Dirac::Matrix grC  = rbA.elementwise_product(bgB);
    grC               -= rgA.elementwise_product(bbB);
    grC               -= bbA.elementwise_product(rgB);
    grC               += bgA.elementwise_product(rbB);

    Dirac::Matrix ggC  = rrA.elementwise_product(bbB);
    ggC               -= rbA.elementwise_product(brB);
    ggC               -= brA.elementwise_product(rbB);
    ggC               += bbA.elementwise_product(rrB);

    Dirac::Matrix gbC  = rgA.elementwise_product(brB);
    gbC               -= rrA.elementwise_product(bgB);
    gbC               -= bgA.elementwise_product(rrB);
    gbC               += brA.elementwise_product(rgB);

    /* ----------------------- */

    Dirac::Matrix brC  = rgA.elementwise_product(gbB);
    brC               -= rbA.elementwise_product(ggB);
    brC               -= ggA.elementwise_product(rbB);
    brC               += gbA.elementwise_product(rgB);

    Dirac::Matrix bgC  = rbA.elementwise_product(grB);
    bgC               -= rrA.elementwise_product(gbB);
    bgC               -= gbA.elementwise_product(rrB);
    bgC               += grA.elementwise_product(rbB);

    Dirac::Matrix bbC  = rrA.elementwise_product(ggB);
    bbC               -= rgA.elementwise_product(grB);
    bbC               -= grA.elementwise_product(rgB);
    bbC               += ggA.elementwise_product(rrB);

    /* --------------------------------------------- */

    size_t index = 12*Base::col_RED  + Base::col_RED;
    d_data[index      ] =  rrC[ 0];
    d_data[index +   3] =  rrC[ 1];
    d_data[index +   6] =  rrC[ 2];
    d_data[index +   9] =  rrC[ 3];
    d_data[index +  36] =  rrC[ 4];
    d_data[index +  39] =  rrC[ 5];
    d_data[index +  42] =  rrC[ 6];
    d_data[index +  45] =  rrC[ 7];
    d_data[index +  72] =  rrC[ 8];
    d_data[index +  75] =  rrC[ 9];
    d_data[index +  78] =  rrC[10];
    d_data[index +  81] =  rrC[11];
    d_data[index + 108] =  rrC[12];
    d_data[index + 111] =  rrC[13];
    d_data[index + 114] =  rrC[14];
    d_data[index + 117] =  rrC[15];

    index = 12*Base::col_RED  + Base::col_GREEN;
    d_data[index      ] =  rgC[ 0];
    d_data[index +   3] =  rgC[ 1];
    d_data[index +   6] =  rgC[ 2];
    d_data[index +   9] =  rgC[ 3];
    d_data[index +  36] =  rgC[ 4];
    d_data[index +  39] =  rgC[ 5];
    d_data[index +  42] =  rgC[ 6];
    d_data[index +  45] =  rgC[ 7];
    d_data[index +  72] =  rgC[ 8];
    d_data[index +  75] =  rgC[ 9];
    d_data[index +  78] =  rgC[10];
    d_data[index +  81] =  rgC[11];
    d_data[index + 108] =  rgC[12];
    d_data[index + 111] =  rgC[13];
    d_data[index + 114] =  rgC[14];
    d_data[index + 117] =  rgC[15];

    index = 12*Base::col_RED  + Base::col_BLUE;
    d_data[index      ] =  rbC[ 0];
    d_data[index +   3] =  rbC[ 1];
    d_data[index +   6] =  rbC[ 2];
    d_data[index +   9] =  rbC[ 3];
    d_data[index +  36] =  rbC[ 4];
    d_data[index +  39] =  rbC[ 5];
    d_data[index +  42] =  rbC[ 6];
    d_data[index +  45] =  rbC[ 7];
    d_data[index +  72] =  rbC[ 8];
    d_data[index +  75] =  rbC[ 9];
    d_data[index +  78] =  rbC[10];
    d_data[index +  81] =  rbC[11];
    d_data[index + 108] =  rbC[12];
    d_data[index + 111] =  rbC[13];
    d_data[index + 114] =  rbC[14];
    d_data[index + 117] =  rbC[15];

    /* --------------------------------------------- */

    index = 12*Base::col_GREEN  + Base::col_RED;
    d_data[index      ] =  grC[ 0];
    d_data[index +   3] =  grC[ 1];
    d_data[index +   6] =  grC[ 2];
    d_data[index +   9] =  grC[ 3];
    d_data[index +  36] =  grC[ 4];
    d_data[index +  39] =  grC[ 5];
    d_data[index +  42] =  grC[ 6];
    d_data[index +  45] =  grC[ 7];
    d_data[index +  72] =  grC[ 8];
    d_data[index +  75] =  grC[ 9];
    d_data[index +  78] =  grC[10];
    d_data[index +  81] =  grC[11];
    d_data[index + 108] =  grC[12];
    d_data[index + 111] =  grC[13];
    d_data[index + 114] =  grC[14];
    d_data[index + 117] =  grC[15];

    index = 12*Base::col_GREEN  + Base::col_GREEN;
    d_data[index      ] =  ggC[ 0];
    d_data[index +   3] =  ggC[ 1];
    d_data[index +   6] =  ggC[ 2];
    d_data[index +   9] =  ggC[ 3];
    d_data[index +  36] =  ggC[ 4];
    d_data[index +  39] =  ggC[ 5];
    d_data[index +  42] =  ggC[ 6];
    d_data[index +  45] =  ggC[ 7];
    d_data[index +  72] =  ggC[ 8];
    d_data[index +  75] =  ggC[ 9];
    d_data[index +  78] =  ggC[10];
    d_data[index +  81] =  ggC[11];
    d_data[index + 108] =  ggC[12];
    d_data[index + 111] =  ggC[13];
    d_data[index + 114] =  ggC[14];
    d_data[index + 117] =  ggC[15];

    index = 12*Base::col_GREEN  + Base::col_BLUE;
    d_data[index      ] =  gbC[ 0];
    d_data[index +   3] =  gbC[ 1];
    d_data[index +   6] =  gbC[ 2];
    d_data[index +   9] =  gbC[ 3];
    d_data[index +  36] =  gbC[ 4];
    d_data[index +  39] =  gbC[ 5];
    d_data[index +  42] =  gbC[ 6];
    d_data[index +  45] =  gbC[ 7];
    d_data[index +  72] =  gbC[ 8];
    d_data[index +  75] =  gbC[ 9];
    d_data[index +  78] =  gbC[10];
    d_data[index +  81] =  gbC[11];
    d_data[index + 108] =  gbC[12];
    d_data[index + 111] =  gbC[13];
    d_data[index + 114] =  gbC[14];
    d_data[index + 117] =  gbC[15];

    /* --------------------------------------------- */

    index = 12*Base::col_BLUE  + Base::col_RED;
    d_data[index      ] =  brC[ 0];
    d_data[index +   3] =  brC[ 1];
    d_data[index +   6] =  brC[ 2];
    d_data[index +   9] =  brC[ 3];
    d_data[index +  36] =  brC[ 4];
    d_data[index +  39] =  brC[ 5];
    d_data[index +  42] =  brC[ 6];
    d_data[index +  45] =  brC[ 7];
    d_data[index +  72] =  brC[ 8];
    d_data[index +  75] =  brC[ 9];
    d_data[index +  78] =  brC[10];
    d_data[index +  81] =  brC[11];
    d_data[index + 108] =  brC[12];
    d_data[index + 111] =  brC[13];
    d_data[index + 114] =  brC[14];
    d_data[index + 117] =  brC[15];

    index = 12*Base::col_BLUE  + Base::col_GREEN;
    d_data[index      ] =  bgC[ 0];
    d_data[index +   3] =  bgC[ 1];
    d_data[index +   6] =  bgC[ 2];
    d_data[index +   9] =  bgC[ 3];
    d_data[index +  36] =  bgC[ 4];
    d_data[index +  39] =  bgC[ 5];
    d_data[index +  42] =  bgC[ 6];
    d_data[index +  45] =  bgC[ 7];
    d_data[index +  72] =  bgC[ 8];
    d_data[index +  75] =  bgC[ 9];
    d_data[index +  78] =  bgC[10];
    d_data[index +  81] =  bgC[11];
    d_data[index + 108] =  bgC[12];
    d_data[index + 111] =  bgC[13];
    d_data[index + 114] =  bgC[14];
    d_data[index + 117] =  bgC[15];

    index = 12*Base::col_BLUE  + Base::col_BLUE;
    d_data[index      ] =  bbC[ 0];
    d_data[index +   3] =  bbC[ 1];
    d_data[index +   6] =  bbC[ 2];
    d_data[index +   9] =  bbC[ 3];
    d_data[index +  36] =  bbC[ 4];
    d_data[index +  39] =  bbC[ 5];
    d_data[index +  42] =  bbC[ 6];
    d_data[index +  45] =  bbC[ 7];
    d_data[index +  72] =  bbC[ 8];
    d_data[index +  75] =  bbC[ 9];
    d_data[index +  78] =  bbC[10];
    d_data[index +  81] =  bbC[11];
    d_data[index + 108] =  bbC[12];
    d_data[index + 111] =  bbC[13];
    d_data[index + 114] =  bbC[14];
    d_data[index + 117] =  bbC[15];

//     std::cout << *this << std::endl;
//     exit(1);
  }




  void make_sequential(Tensor result[16], Tensor const &A, Tensor const &B)
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

    size_t const n_complex = 16*16;

    std::complex< double > tmp[n_complex];

    std::complex< double > zero(0, 0);

    //Dirac::Matrix *tmpTensors = reinterpret_cast< Dirac::Matrix* >(tmp);

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
//     Dirac::Matrix *tmp[9];
//     result = Tensor(tmp);
//   }


}

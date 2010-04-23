#include "Tensor.ih"

namespace QCD
{

  // this is for the d_bar * Op * d current (and it keeps all indices free)
  void make_sequential_d(Tensor result[16], Tensor const &A, Tensor const &B)
  {
    /* under construction */
    // assert(false);
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
    ggA.outer_product(bbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    gbA.outer_product(bgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bgA.outer_product(gbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bbA.outer_product(ggB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    ggA.outer_product(bbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    gbA.outer_product(bgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bgA.outer_product(gbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bbA.outer_product(ggB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());


    std::complex< double > rgC[n_complex];
    std::fill_n(rgC, n_complex, zero);
    gbA.outer_product(brB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    grA.outer_product(bbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    bbA.outer_product(grB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    brA.outer_product(gbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    gbA.outer_product(brB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    grA.outer_product(bbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    bbA.outer_product(grB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    brA.outer_product(gbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());

    std::complex< double > rbC[n_complex];
    std::fill_n(rbC, n_complex, zero);
    grA.outer_product(bgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    ggA.outer_product(brB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    brA.outer_product(ggB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    bgA.outer_product(grB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    grA.outer_product(bgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    ggA.outer_product(brB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    brA.outer_product(ggB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    bgA.outer_product(grB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > grC[n_complex];
    std::fill_n(grC, n_complex, zero);
    rbA.outer_product(bgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rgA.outer_product(bbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bbA.outer_product(rgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bgA.outer_product(rbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rbA.outer_product(bgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rgA.outer_product(bbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bbA.outer_product(rgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bgA.outer_product(rbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());

    std::complex< double > ggC[n_complex];
    std::fill_n(ggC, n_complex, zero);
    rrA.outer_product(bbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rbA.outer_product(brB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    brA.outer_product(rbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    bbA.outer_product(rrB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rrA.outer_product(bbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rbA.outer_product(brB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    brA.outer_product(rbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    bbA.outer_product(rrB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());

    std::complex< double > gbC[n_complex];
    std::fill_n(gbC, n_complex, zero);
    rgA.outer_product(brB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rrA.outer_product(bgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    bgA.outer_product(rrB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    brA.outer_product(rgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rgA.outer_product(brB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rrA.outer_product(bgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    bgA.outer_product(rrB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    brA.outer_product(rgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > brC[n_complex];
    std::fill_n(brC, n_complex, zero);
    rgA.outer_product(gbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rbA.outer_product(ggB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    ggA.outer_product(rbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    gbA.outer_product(rgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rgA.outer_product(gbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rbA.outer_product(ggB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    ggA.outer_product(rbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    gbA.outer_product(rgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());

    std::complex< double > bgC[n_complex];
    std::fill_n(bgC, n_complex, zero);
    rbA.outer_product(grB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rrA.outer_product(gbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    gbA.outer_product(rrB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    grA.outer_product(rbB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rbA.outer_product(grB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rrA.outer_product(gbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    gbA.outer_product(rrB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    grA.outer_product(rbB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());

    std::complex< double > bbC[n_complex];
    std::fill_n(bbC, n_complex, zero);
    rrA.outer_product(ggB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rgA.outer_product(grB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    grA.outer_product(rgB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    ggA.outer_product(rrB, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rrA.outer_product(ggB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rgA.outer_product(grB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    grA.outer_product(rgB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    ggA.outer_product(rrB, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());

//     for (size_t iDirac=0; iDirac<16*16; iDirac++)
//     {
//       if (iDirac % 16 == 0)
//         std::cout << std::endl;
//       std::cout << tmp[iDirac] << " ";
//     }
//     std::cout << std::endl;

    /* --------------------------------------------- */

    // the following part takes 9 Dirac::Matrixes and translates them to a Tensor for each od the 16 Dirac combinations

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

      data_ptr[  0] =  rbC[offset];
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

      data_ptr[  0] =  grC[offset];
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

      data_ptr[  0] =  ggC[offset];
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

      data_ptr[  0] =  gbC[offset];
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

      data_ptr[  0] =  brC[offset];
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

      data_ptr[  0] =  bgC[offset];
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

      data_ptr[  0] =  bbC[offset];
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
  }

}

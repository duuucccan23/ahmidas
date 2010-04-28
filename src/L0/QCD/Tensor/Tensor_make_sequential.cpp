#include "Tensor.ih"

namespace QCD
{

  // this is for the d_bar * Op * d current (and it keeps all indices free)
  void make_sequential_d(Tensor result[16], Tensor const &U1, Tensor const &U2)
  {
    /* under construction */
    // assert(false);
    Dirac::Matrix rrU1;
    Dirac::Matrix rgU1;
    Dirac::Matrix rbU1;
    Dirac::Matrix grU1;
    Dirac::Matrix gbU1;
    Dirac::Matrix ggU1;
    Dirac::Matrix brU1;
    Dirac::Matrix bbU1;
    Dirac::Matrix bgU1;

    Dirac::Matrix rrU2;
    Dirac::Matrix rgU2;
    Dirac::Matrix rbU2;
    Dirac::Matrix grU2;
    Dirac::Matrix gbU2;
    Dirac::Matrix ggU2;
    Dirac::Matrix brU2;
    Dirac::Matrix bbU2;
    Dirac::Matrix bgU2;

    U1.getDiracMatrix(rrU1, Base::col_RED,   Base::col_RED);
    U1.getDiracMatrix(rgU1, Base::col_RED,   Base::col_GREEN);
    U1.getDiracMatrix(rbU1, Base::col_RED,   Base::col_BLUE);
    U1.getDiracMatrix(grU1, Base::col_GREEN, Base::col_RED);
    U1.getDiracMatrix(gbU1, Base::col_GREEN, Base::col_BLUE);
    U1.getDiracMatrix(ggU1, Base::col_GREEN, Base::col_GREEN);
    U1.getDiracMatrix(brU1, Base::col_BLUE,  Base::col_RED);
    U1.getDiracMatrix(bbU1, Base::col_BLUE,  Base::col_BLUE);
    U1.getDiracMatrix(bgU1, Base::col_BLUE,  Base::col_GREEN);

    U2.getDiracMatrix(rrU2, Base::col_RED,   Base::col_RED);
    U2.getDiracMatrix(rgU2, Base::col_RED,   Base::col_GREEN);
    U2.getDiracMatrix(rbU2, Base::col_RED,   Base::col_BLUE);
    U2.getDiracMatrix(grU2, Base::col_GREEN, Base::col_RED);
    U2.getDiracMatrix(ggU2, Base::col_GREEN, Base::col_GREEN);
    U2.getDiracMatrix(gbU2, Base::col_GREEN, Base::col_BLUE);
    U2.getDiracMatrix(brU2, Base::col_BLUE,  Base::col_RED);
    U2.getDiracMatrix(bgU2, Base::col_BLUE,  Base::col_GREEN);
    U2.getDiracMatrix(bbU2, Base::col_BLUE,  Base::col_BLUE);

    size_t const n_complex = 16*16;

    std::complex< double > tmp[n_complex];

    std::complex< double > zero(0, 0);

    //Dirac::Matrix *tmpTensors = reinterpret_cast< Dirac::Matrix* >(tmp);

    std::complex< double > *tmpTensors = tmp;

    std::complex< double > rrC[n_complex];
    std::fill_n(rrC, n_complex, zero);
    ggU1.outer_product(bbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    gbU1.outer_product(bgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bgU1.outer_product(gbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bbU1.outer_product(ggU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    ggU1.outer_product(bbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());
    gbU1.outer_product(bgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bgU1.outer_product(gbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::minus< std::complex< double> >());
    bbU1.outer_product(ggU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rrC, rrC, std::plus< std::complex< double> >());


    std::complex< double > rgC[n_complex];
    std::fill_n(rgC, n_complex, zero);
    gbU1.outer_product(brU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    grU1.outer_product(bbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    bbU1.outer_product(grU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    brU1.outer_product(gbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    gbU1.outer_product(brU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());
    grU1.outer_product(bbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    bbU1.outer_product(grU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::minus< std::complex< double> >());
    brU1.outer_product(gbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rgC, rgC, std::plus< std::complex< double> >());

    std::complex< double > rbC[n_complex];
    std::fill_n(rbC, n_complex, zero);
    grU1.outer_product(bgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    ggU1.outer_product(brU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    brU1.outer_product(ggU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    bgU1.outer_product(grU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    grU1.outer_product(bgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());
    ggU1.outer_product(brU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    brU1.outer_product(ggU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::minus< std::complex< double> >());
    bgU1.outer_product(grU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, rbC, rbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > grC[n_complex];
    std::fill_n(grC, n_complex, zero);
    rbU1.outer_product(bgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rgU1.outer_product(bbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bbU1.outer_product(rgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bgU1.outer_product(rbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rbU1.outer_product(bgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());
    rgU1.outer_product(bbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bbU1.outer_product(rgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::minus< std::complex< double> >());
    bgU1.outer_product(rbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, grC, grC, std::plus< std::complex< double> >());

    std::complex< double > ggC[n_complex];
    std::fill_n(ggC, n_complex, zero);
    rrU1.outer_product(bbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rbU1.outer_product(brU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    brU1.outer_product(rbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    bbU1.outer_product(rrU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rrU1.outer_product(bbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());
    rbU1.outer_product(brU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    brU1.outer_product(rbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::minus< std::complex< double> >());
    bbU1.outer_product(rrU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, ggC, ggC, std::plus< std::complex< double> >());

    std::complex< double > gbC[n_complex];
    std::fill_n(gbC, n_complex, zero);
    rgU1.outer_product(brU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rrU1.outer_product(bgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    bgU1.outer_product(rrU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    brU1.outer_product(rgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rgU1.outer_product(brU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());
    rrU1.outer_product(bgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    bgU1.outer_product(rrU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::minus< std::complex< double> >());
    brU1.outer_product(rgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, gbC, gbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > brC[n_complex];
    std::fill_n(brC, n_complex, zero);
    rgU1.outer_product(gbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rbU1.outer_product(ggU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    ggU1.outer_product(rbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    gbU1.outer_product(rgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rgU1.outer_product(gbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());
    rbU1.outer_product(ggU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    ggU1.outer_product(rbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::minus< std::complex< double> >());
    gbU1.outer_product(rgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, brC, brC, std::plus< std::complex< double> >());

    std::complex< double > bgC[n_complex];
    std::fill_n(bgC, n_complex, zero);
    rbU1.outer_product(grU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rrU1.outer_product(gbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    gbU1.outer_product(rrU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    grU1.outer_product(rbU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rbU1.outer_product(grU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());
    rrU1.outer_product(gbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    gbU1.outer_product(rrU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::minus< std::complex< double> >());
    grU1.outer_product(rbU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bgC, bgC, std::plus< std::complex< double> >());

    std::complex< double > bbC[n_complex];
    std::fill_n(bbC, n_complex, zero);
    rrU1.outer_product(ggU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rgU1.outer_product(grU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    grU1.outer_product(rgU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    ggU1.outer_product(rrU2, tmpTensors, Dirac::order_FIRST_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rrU1.outer_product(ggU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::plus< std::complex< double> >());
    rgU1.outer_product(grU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    grU1.outer_product(rgU2, tmpTensors, Dirac::order_OUTER_FIXED);
    std::transform(tmp, tmp + n_complex, bbC, bbC, std::minus< std::complex< double> >());
    ggU1.outer_product(rrU2, tmpTensors, Dirac::order_OUTER_FIXED);
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

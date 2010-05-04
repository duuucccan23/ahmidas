#include "Tensor.ih"

#include <vector>

namespace QCD
{

  // this is for the u_bar * Op * u current
  // it is assuming that D comes as (Gamma*D*Gamma_bar)^T, where Gamma is chosen according to interpolating field
  void make_sequential_u(Tensor result[16], Tensor const &D, Tensor const &U)
  {

    // note that D and U have different order when order is order_FIRST_OUTER_DELTA or order_BOTH_OUTER_DELTA than
    // for the other Dirac::OuterIndexOrder(s)

    Dirac::Matrix rrU;
    Dirac::Matrix rgU;
    Dirac::Matrix rbU;
    Dirac::Matrix grU;
    Dirac::Matrix gbU;
    Dirac::Matrix ggU;
    Dirac::Matrix brU;
    Dirac::Matrix bbU;
    Dirac::Matrix bgU;

    Dirac::Matrix rrD;
    Dirac::Matrix rgD;
    Dirac::Matrix rbD;
    Dirac::Matrix grD;
    Dirac::Matrix gbD;
    Dirac::Matrix ggD;
    Dirac::Matrix brD;
    Dirac::Matrix bbD;
    Dirac::Matrix bgD;

    U.getDiracMatrix(rrU, Base::col_RED,   Base::col_RED);
    U.getDiracMatrix(rgU, Base::col_RED,   Base::col_GREEN);
    U.getDiracMatrix(rbU, Base::col_RED,   Base::col_BLUE);
    U.getDiracMatrix(grU, Base::col_GREEN, Base::col_RED);
    U.getDiracMatrix(ggU, Base::col_GREEN, Base::col_GREEN);
    U.getDiracMatrix(gbU, Base::col_GREEN, Base::col_BLUE);
    U.getDiracMatrix(brU, Base::col_BLUE,  Base::col_RED);
    U.getDiracMatrix(bgU, Base::col_BLUE,  Base::col_GREEN);
    U.getDiracMatrix(bbU, Base::col_BLUE,  Base::col_BLUE);

    // implicit color transpose
    D.getDiracMatrix(rrD, Base::col_RED,   Base::col_RED);
    D.getDiracMatrix(rgD, Base::col_GREEN, Base::col_RED);
    D.getDiracMatrix(rbD, Base::col_BLUE,  Base::col_RED);
    D.getDiracMatrix(grD, Base::col_RED,   Base::col_GREEN);
    D.getDiracMatrix(ggD, Base::col_GREEN, Base::col_GREEN);
    D.getDiracMatrix(gbD, Base::col_BLUE,  Base::col_GREEN);
    D.getDiracMatrix(brD, Base::col_RED,   Base::col_BLUE);
    D.getDiracMatrix(bgD, Base::col_GREEN, Base::col_BLUE);
    D.getDiracMatrix(bbD, Base::col_BLUE,  Base::col_BLUE);


    size_t const n_complex = 16*16;

    std::complex< double > tmp[n_complex];

    std::complex< double > const zero(0, 0);

    std::complex< double > rrC[n_complex];
    std::fill_n(rrC, n_complex, zero);
    ggU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    gbU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    bgU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    bbU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    bbD.outer_product(ggU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    bgD.outer_product(gbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    gbD.outer_product(bgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    ggD.outer_product(bbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    bbD.outer_product(ggU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    bgD.outer_product(gbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    gbD.outer_product(bgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    ggD.outer_product(bbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    ggU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());
    gbU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    bgU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::minus< std::complex< double> >());
    bbU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rrC, rrC+ n_complex, tmp, rrC, std::plus< std::complex< double> >());

    std::complex< double > rgC[n_complex];
    std::fill_n(rgC, n_complex, zero);
    gbU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    grU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    bbU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    brU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    brD.outer_product(gbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    bbD.outer_product(grU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    grD.outer_product(bbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    gbD.outer_product(brU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    brD.outer_product(gbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    bbD.outer_product(grU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    grD.outer_product(bbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    gbD.outer_product(brU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    gbU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());
    grU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    bbU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::minus< std::complex< double> >());
    brU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rgC, rgC+ n_complex, tmp, rgC, std::plus< std::complex< double> >());

    std::complex< double > rbC[n_complex];
    std::fill_n(rbC, n_complex, zero);
    grU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    ggU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    brU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    bgU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    bgD.outer_product(grU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    brD.outer_product(ggU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    ggD.outer_product(brU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    grD.outer_product(bgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    bgD.outer_product(grU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    brD.outer_product(ggU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    ggD.outer_product(brU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    grD.outer_product(bgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    grU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());
    ggU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    brU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::minus< std::complex< double> >());
    bgU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(rbC, rbC+ n_complex, tmp, rbC, std::plus< std::complex< double> >());

    /* ----------------------- */


    std::complex< double > grC[n_complex];
    std::fill_n(grC, n_complex, zero);
    rbU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    rgU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    bbU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    bgU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    bgD.outer_product(rbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    bbD.outer_product(rgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    rgD.outer_product(bbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    rbD.outer_product(bgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    bgD.outer_product(rbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    bbD.outer_product(rgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    rgD.outer_product(bbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    rbD.outer_product(bgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    rbU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());
    rgU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    bbU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::minus< std::complex< double> >());
    bgU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(grC, grC+ n_complex, tmp, grC, std::plus< std::complex< double> >());

    std::complex< double > ggC[n_complex];
    std::fill_n(ggC, n_complex, zero);
    rrU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    rbU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    brU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    bbU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    bbD.outer_product(rrU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    brD.outer_product(rbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    rbD.outer_product(brU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    rrD.outer_product(bbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    bbD.outer_product(rrU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    brD.outer_product(rbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    rbD.outer_product(brU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    rrD.outer_product(bbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    rrU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());
    rbU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    brU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::minus< std::complex< double> >());
    bbU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(ggC, ggC+ n_complex, tmp, ggC, std::plus< std::complex< double> >());

    std::complex< double > gbC[n_complex];
    std::fill_n(gbC, n_complex, zero);
    rgU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    rrU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    bgU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    brU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    brD.outer_product(rgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    bgD.outer_product(rrU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    rrD.outer_product(bgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    rgD.outer_product(brU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    brD.outer_product(rgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    bgD.outer_product(rrU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    rrD.outer_product(bgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    rgD.outer_product(brU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    rgU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());
    rrU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    bgU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::minus< std::complex< double> >());
    brU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(gbC, gbC+ n_complex, tmp, gbC, std::plus< std::complex< double> >());

    /* ----------------------- */

    std::complex< double > brC[n_complex];
    std::fill_n(brC, n_complex, zero);
    rgU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    rbU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    ggU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    gbU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    gbD.outer_product(rgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    ggD.outer_product(rbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    rbD.outer_product(ggU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    rgD.outer_product(gbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    gbD.outer_product(rgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    ggD.outer_product(rbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    rbD.outer_product(ggU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    rgD.outer_product(gbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    rgU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());
    rbU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    ggU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::minus< std::complex< double> >());
    gbU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(brC, brC+ n_complex, tmp, brC, std::plus< std::complex< double> >());

    std::complex< double > bgC[n_complex];
    std::fill_n(bgC, n_complex, zero);
    rbU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    rrU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    gbU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    grU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    grD.outer_product(rbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    gbD.outer_product(rrU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    rrD.outer_product(gbU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    rbD.outer_product(grU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    grD.outer_product(rbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    gbD.outer_product(rrU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    rrD.outer_product(gbU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    rbD.outer_product(grU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    rbU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());
    rrU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    gbU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::minus< std::complex< double> >());
    grU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bgC, bgC+ n_complex, tmp, bgC, std::plus< std::complex< double> >());

    std::complex< double > bbC[n_complex];
    std::fill_n(bbC, n_complex, zero);
    rrU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    rgU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    grU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    ggU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    ggD.outer_product(rrU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    grD.outer_product(rgU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    rgD.outer_product(grU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    rrD.outer_product(ggU, tmp, Dirac::order_BOTH_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    ggD.outer_product(rrU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    grD.outer_product(rgU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    rgD.outer_product(grU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    rrD.outer_product(ggU, tmp, Dirac::order_FIRST_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    rrU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());
    rgU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    grU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::minus< std::complex< double> >());
    ggU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA);
    std::transform(bbC, bbC+ n_complex, tmp, bbC, std::plus< std::complex< double> >());

    /* --------------------------------------------- */


    // the following part takes 9 Dirac::Matrixes and translates them to a Tensor for each of the 16 Dirac combinations

    std::complex< double >* data_ptr = NULL;

    size_t index(-1);

    for (size_t iDirac=0; iDirac<16; iDirac++)
    {
      size_t const offset = iDirac * 16;

      index = 12*Base::col_RED  + Base::col_RED;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  rrC[offset +  0];
      data_ptr[  3] =  rrC[offset +  1];
      data_ptr[  6] =  rrC[offset +  2];
      data_ptr[  9] =  rrC[offset +  3];
      data_ptr[ 36] =  rrC[offset +  4];
      data_ptr[ 39] =  rrC[offset +  5];
      data_ptr[ 42] =  rrC[offset +  6];
      data_ptr[ 45] =  rrC[offset +  7];
      data_ptr[ 72] =  rrC[offset +  8];
      data_ptr[ 75] =  rrC[offset +  9];
      data_ptr[ 78] =  rrC[offset + 10];
      data_ptr[ 81] =  rrC[offset + 11];
      data_ptr[108] =  rrC[offset + 12];
      data_ptr[111] =  rrC[offset + 13];
      data_ptr[114] =  rrC[offset + 14];
      data_ptr[117] =  rrC[offset + 15];

      index = 12*Base::col_GREEN + Base::col_RED;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  rgC[offset +  0];
      data_ptr[  3] =  rgC[offset +  1];
      data_ptr[  6] =  rgC[offset +  2];
      data_ptr[  9] =  rgC[offset +  3];
      data_ptr[ 36] =  rgC[offset +  4];
      data_ptr[ 39] =  rgC[offset +  5];
      data_ptr[ 42] =  rgC[offset +  6];
      data_ptr[ 45] =  rgC[offset +  7];
      data_ptr[ 72] =  rgC[offset +  8];
      data_ptr[ 75] =  rgC[offset +  9];
      data_ptr[ 78] =  rgC[offset + 10];
      data_ptr[ 81] =  rgC[offset + 11];
      data_ptr[108] =  rgC[offset + 12];
      data_ptr[111] =  rgC[offset + 13];
      data_ptr[114] =  rgC[offset + 14];
      data_ptr[117] =  rgC[offset + 15];

      index = 12*Base::col_BLUE + Base::col_RED;
      data_ptr = (result[iDirac]).d_data + index;;

      data_ptr[  0] =  rbC[offset +  0];
      data_ptr[  3] =  rbC[offset +  1];
      data_ptr[  6] =  rbC[offset +  2];
      data_ptr[  9] =  rbC[offset +  3];
      data_ptr[ 36] =  rbC[offset +  4];
      data_ptr[ 39] =  rbC[offset +  5];
      data_ptr[ 42] =  rbC[offset +  6];
      data_ptr[ 45] =  rbC[offset +  7];
      data_ptr[ 72] =  rbC[offset +  8];
      data_ptr[ 75] =  rbC[offset +  9];
      data_ptr[ 78] =  rbC[offset + 10];
      data_ptr[ 81] =  rbC[offset + 11];
      data_ptr[108] =  rbC[offset + 12];
      data_ptr[111] =  rbC[offset + 13];
      data_ptr[114] =  rbC[offset + 14];
      data_ptr[117] =  rbC[offset + 15];

      /* --------------------------------------------- */

      index = 12* Base::col_RED + Base::col_GREEN;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  grC[offset +  0];
      data_ptr[  3] =  grC[offset +  1];
      data_ptr[  6] =  grC[offset +  2];
      data_ptr[  9] =  grC[offset +  3];
      data_ptr[ 36] =  grC[offset +  4];
      data_ptr[ 39] =  grC[offset +  5];
      data_ptr[ 42] =  grC[offset +  6];
      data_ptr[ 45] =  grC[offset +  7];
      data_ptr[ 72] =  grC[offset +  8];
      data_ptr[ 75] =  grC[offset +  9];
      data_ptr[ 78] =  grC[offset + 10];
      data_ptr[ 81] =  grC[offset + 11];
      data_ptr[108] =  grC[offset + 12];
      data_ptr[111] =  grC[offset + 13];
      data_ptr[114] =  grC[offset + 14];
      data_ptr[117] =  grC[offset + 15];

      index = 12*Base::col_GREEN  + Base::col_GREEN;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  ggC[offset +  0];
      data_ptr[  3] =  ggC[offset +  1];
      data_ptr[  6] =  ggC[offset +  2];
      data_ptr[  9] =  ggC[offset +  3];
      data_ptr[ 36] =  ggC[offset +  4];
      data_ptr[ 39] =  ggC[offset +  5];
      data_ptr[ 42] =  ggC[offset +  6];
      data_ptr[ 45] =  ggC[offset +  7];
      data_ptr[ 72] =  ggC[offset +  8];
      data_ptr[ 75] =  ggC[offset +  9];
      data_ptr[ 78] =  ggC[offset + 10];
      data_ptr[ 81] =  ggC[offset + 11];
      data_ptr[108] =  ggC[offset + 12];
      data_ptr[111] =  ggC[offset + 13];
      data_ptr[114] =  ggC[offset + 14];
      data_ptr[117] =  ggC[offset + 15];

      index = 12*Base::col_BLUE + Base::col_GREEN;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  gbC[offset +  0];
      data_ptr[  3] =  gbC[offset +  1];
      data_ptr[  6] =  gbC[offset +  2];
      data_ptr[  9] =  gbC[offset +  3];
      data_ptr[ 36] =  gbC[offset +  4];
      data_ptr[ 39] =  gbC[offset +  5];
      data_ptr[ 42] =  gbC[offset +  6];
      data_ptr[ 45] =  gbC[offset +  7];
      data_ptr[ 72] =  gbC[offset +  8];
      data_ptr[ 75] =  gbC[offset +  9];
      data_ptr[ 78] =  gbC[offset + 10];
      data_ptr[ 81] =  gbC[offset + 11];
      data_ptr[108] =  gbC[offset + 12];
      data_ptr[111] =  gbC[offset + 13];
      data_ptr[114] =  gbC[offset + 14];
      data_ptr[117] =  gbC[offset + 15];

      /* --------------------------------------------- */

      index = 12*Base::col_RED + Base::col_BLUE;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  brC[offset +  0];
      data_ptr[  3] =  brC[offset +  1];
      data_ptr[  6] =  brC[offset +  2];
      data_ptr[  9] =  brC[offset +  3];
      data_ptr[ 36] =  brC[offset +  4];
      data_ptr[ 39] =  brC[offset +  5];
      data_ptr[ 42] =  brC[offset +  6];
      data_ptr[ 45] =  brC[offset +  7];
      data_ptr[ 72] =  brC[offset +  8];
      data_ptr[ 75] =  brC[offset +  9];
      data_ptr[ 78] =  brC[offset + 10];
      data_ptr[ 81] =  brC[offset + 11];
      data_ptr[108] =  brC[offset + 12];
      data_ptr[111] =  brC[offset + 13];
      data_ptr[114] =  brC[offset + 14];
      data_ptr[117] =  brC[offset + 15];

      index = 12*Base::col_GREEN + Base::col_BLUE;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  bgC[offset +  0];
      data_ptr[  3] =  bgC[offset +  1];
      data_ptr[  6] =  bgC[offset +  2];
      data_ptr[  9] =  bgC[offset +  3];
      data_ptr[ 36] =  bgC[offset +  4];
      data_ptr[ 39] =  bgC[offset +  5];
      data_ptr[ 42] =  bgC[offset +  6];
      data_ptr[ 45] =  bgC[offset +  7];
      data_ptr[ 72] =  bgC[offset +  8];
      data_ptr[ 75] =  bgC[offset +  9];
      data_ptr[ 78] =  bgC[offset + 10];
      data_ptr[ 81] =  bgC[offset + 11];
      data_ptr[108] =  bgC[offset + 12];
      data_ptr[111] =  bgC[offset + 13];
      data_ptr[114] =  bgC[offset + 14];
      data_ptr[117] =  bgC[offset + 15];

      index = 12*Base::col_BLUE  + Base::col_BLUE;
      data_ptr = (result[iDirac]).d_data + index;

      data_ptr[  0] =  bbC[offset +  0];
      data_ptr[  3] =  bbC[offset +  1];
      data_ptr[  6] =  bbC[offset +  2];
      data_ptr[  9] =  bbC[offset +  3];
      data_ptr[ 36] =  bbC[offset +  4];
      data_ptr[ 39] =  bbC[offset +  5];
      data_ptr[ 42] =  bbC[offset +  6];
      data_ptr[ 45] =  bbC[offset +  7];
      data_ptr[ 72] =  bbC[offset +  8];
      data_ptr[ 75] =  bbC[offset +  9];
      data_ptr[ 78] =  bbC[offset + 10];
      data_ptr[ 81] =  bbC[offset + 11];
      data_ptr[108] =  bbC[offset + 12];
      data_ptr[111] =  bbC[offset + 13];
      data_ptr[114] =  bbC[offset + 14];
      data_ptr[117] =  bbC[offset + 15];

      data_ptr = NULL;
    }
  }

}

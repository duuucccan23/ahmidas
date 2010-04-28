#include "Tensor.ih"

#include <vector>

namespace QCD
{

  // this is for the u_bar * Op * u current
  // it is assuming that D comes as (Gamma*D*Gamma_bar)^T, where Gamma is chosen according to interpolating field
  void make_sequential_u(Tensor &result, Tensor const &D, Tensor const &U, Base::BaryonPropagatorProjector const projector)
  {

    std::vector< Base::DiracIndex > alpha_f_std;
    alpha_f_std.push_back(Base::gam_1);
    alpha_f_std.push_back(Base::gam_2);
    alpha_f_std.push_back(Base::gam_3);
    alpha_f_std.push_back(Base::gam_4);
    std::vector< Base::DiracIndex > alpha_f_perms;
    alpha_f_perms.reserve(8);
    std::vector< std::complex< double > > alpha_f_factors;
    alpha_f_factors.reserve(8);

    static std::complex< double > const COMPLEX_P_I( 0,  1);
    static std::complex< double > const COMPLEX_M_I( 0, -1);
    static std::complex< double > const COMPLEX_P_1( 1,  0);
    static std::complex< double > const COMPLEX_M_1(-1,  0);
    static std::complex< double > const COMPLEX_0(0, 0);

    switch (projector)
    {
      // note that the projector has to be daggered
      case Base::proj_PARITY_PLUS_TM:
        alpha_f_perms.push_back(Base::gam_3);
        alpha_f_factors.push_back(0.5*COMPLEX_M_1);
        alpha_f_perms.push_back(Base::gam_4);
        alpha_f_factors.push_back(0.5*COMPLEX_M_1);
        alpha_f_perms.push_back(Base::gam_1);
        alpha_f_factors.push_back(0.5*COMPLEX_M_1);
        alpha_f_perms.push_back(Base::gam_2);
        alpha_f_factors.push_back(0.5*COMPLEX_M_1);
        alpha_f_perms.push_back(Base::gam_1);
        alpha_f_factors.push_back(0.5*COMPLEX_P_I);
        alpha_f_perms.push_back(Base::gam_2);
        alpha_f_factors.push_back(0.5*COMPLEX_P_I);
        alpha_f_perms.push_back(Base::gam_3);
        alpha_f_factors.push_back(0.5*COMPLEX_M_I);
        alpha_f_perms.push_back(Base::gam_4);
        alpha_f_factors.push_back(0.5*COMPLEX_M_I);
      break;
      default:
      std::cerr << "Unknown projector in QCD::make_sequential_d() " << std::endl;
    }

    std::fill_n(result.d_data, result.size(), COMPLEX_0);

    /* under construction */
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
    U.getDiracMatrix(gbU, Base::col_GREEN, Base::col_BLUE);
    U.getDiracMatrix(ggU, Base::col_GREEN, Base::col_GREEN);
    U.getDiracMatrix(brU, Base::col_BLUE,  Base::col_RED);
    U.getDiracMatrix(bbU, Base::col_BLUE,  Base::col_BLUE);
    U.getDiracMatrix(bgU, Base::col_BLUE,  Base::col_GREEN);

    D.getDiracMatrix(rrD, Base::col_RED,   Base::col_RED);
    D.getDiracMatrix(rgD, Base::col_GREEN, Base::col_RED);
    D.getDiracMatrix(rbD, Base::col_BLUE,  Base::col_RED);
    D.getDiracMatrix(grD, Base::col_RED,   Base::col_GREEN);
    D.getDiracMatrix(gbD, Base::col_BLUE,  Base::col_GREEN);
    D.getDiracMatrix(ggD, Base::col_GREEN, Base::col_GREEN);
    D.getDiracMatrix(brD, Base::col_RED,   Base::col_BLUE);
    D.getDiracMatrix(bbD, Base::col_BLUE,  Base::col_BLUE);
    D.getDiracMatrix(bgD, Base::col_GREEN, Base::col_BLUE);

    std::complex< double > tmp[16];

    for (size_t idx=0; idx<alpha_f_perms.size(); idx++)
    {
      Base::DiracIndex const alpha_f = alpha_f_perms[idx];
      Base::DiracIndex const alpha_i = alpha_f_std[idx%4];


      std::complex< double > rrC[16];
      std::fill_n(rrC, 16, COMPLEX_0);
      ggU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      gbU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bgU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bbU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      ggU.outer_product(bbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      gbU.outer_product(bgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bgU.outer_product(gbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bbU.outer_product(ggD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      ggU.outer_product(bbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      gbU.outer_product(bgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bgU.outer_product(gbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bbU.outer_product(ggD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      ggU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());
      gbU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bgU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::minus< std::complex< double> >());
      bbU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rrC, rrC, std::plus< std::complex< double> >());

      std::complex< double > rgC[16];
      std::fill_n(rgC, 16, COMPLEX_0);
      gbU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      grU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      bbU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      brU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      gbU.outer_product(brD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      grU.outer_product(bbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      bbU.outer_product(grD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      brU.outer_product(gbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      gbU.outer_product(brD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      grU.outer_product(bbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      bbU.outer_product(grD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      brU.outer_product(gbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      gbU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());
      grU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      bbU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::minus< std::complex< double> >());
      brU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rgC, rgC, std::plus< std::complex< double> >());

      std::complex< double > rbC[16];
      std::fill_n(rbC, 16, COMPLEX_0);
      grU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      ggU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      brU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      bgU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      grU.outer_product(bgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      ggU.outer_product(brD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      brU.outer_product(ggD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      bgU.outer_product(grD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      grU.outer_product(bgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      ggU.outer_product(brD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      brU.outer_product(ggD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      bgU.outer_product(grD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      grU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());
      ggU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      brU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::minus< std::complex< double> >());
      bgU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, rbC, rbC, std::plus< std::complex< double> >());

      /* ----------------------- */

      std::complex< double > grC[16];
      std::fill_n(grC, 16, COMPLEX_0);
      rbU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rgU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bbU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bgU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rbU.outer_product(bgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rgU.outer_product(bbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bbU.outer_product(rgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bgU.outer_product(rbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rbU.outer_product(bgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rgU.outer_product(bbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bbU.outer_product(rgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bgU.outer_product(rbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rbU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());
      rgU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bbU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::minus< std::complex< double> >());
      bgU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, grC, grC, std::plus< std::complex< double> >());

      std::complex< double > ggC[16];
      std::fill_n(ggC, 16, COMPLEX_0);
      rrU.outer_product(bbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rbU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      brU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      bbU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rrU.outer_product(bbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rbU.outer_product(brD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      brU.outer_product(rbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      bbU.outer_product(rrD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rrU.outer_product(bbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rbU.outer_product(brD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      brU.outer_product(rbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      bbU.outer_product(rrD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rrU.outer_product(bbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());
      rbU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      brU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::minus< std::complex< double> >());
      bbU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, ggC, ggC, std::plus< std::complex< double> >());

      std::complex< double > gbC[16];
      std::fill_n(gbC, 16, COMPLEX_0);
      rgU.outer_product(brD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rrU.outer_product(bgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      bgU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      brU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rgU.outer_product(brD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rrU.outer_product(bgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      bgU.outer_product(rrD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      brU.outer_product(rgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rgU.outer_product(brD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rrU.outer_product(bgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      bgU.outer_product(rrD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      brU.outer_product(rgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rgU.outer_product(brD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());
      rrU.outer_product(bgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      bgU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::minus< std::complex< double> >());
      brU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, gbC, gbC, std::plus< std::complex< double> >());

      /* ----------------------- */

      std::complex< double > brC[16];
      std::fill_n(brC, 16, COMPLEX_0);
      rgU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rbU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      ggU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      gbU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rgU.outer_product(gbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rbU.outer_product(ggD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      ggU.outer_product(rbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      gbU.outer_product(rgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rgU.outer_product(gbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rbU.outer_product(ggD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      ggU.outer_product(rbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      gbU.outer_product(rgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rgU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());
      rbU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      ggU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::minus< std::complex< double> >());
      gbU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, brC, brC, std::plus< std::complex< double> >());

      std::complex< double > bgC[16];
      std::fill_n(bgC, 16, COMPLEX_0);
      rbU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rrU.outer_product(gbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      gbU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      grU.outer_product(rbD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rbU.outer_product(grD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rrU.outer_product(gbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      gbU.outer_product(rrD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      grU.outer_product(rbD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rbU.outer_product(grD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rrU.outer_product(gbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      gbU.outer_product(rrD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      grU.outer_product(rbD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rbU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());
      rrU.outer_product(gbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      gbU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::minus< std::complex< double> >());
      grU.outer_product(rbD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bgC, bgC, std::plus< std::complex< double> >());

      std::complex< double > bbC[16];
      std::fill_n(bbC, 16, COMPLEX_0);
      rrU.outer_product(ggD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rgU.outer_product(grD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      grU.outer_product(rgD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      ggU.outer_product(rrD, tmp, Dirac::order_FIRST_FIXED, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rrU.outer_product(ggD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rgU.outer_product(grD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      grU.outer_product(rgD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      ggU.outer_product(rrD, tmp, Dirac::order_BOTH_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rrU.outer_product(ggD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rgU.outer_product(grD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      grU.outer_product(rgD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      ggU.outer_product(rrD, tmp, Dirac::order_FIRST_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rrU.outer_product(ggD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());
      rgU.outer_product(grD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      grU.outer_product(rgD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::minus< std::complex< double> >());
      ggU.outer_product(rrD, tmp, Dirac::order_SECOND_OUTER_DELTA, alpha_i, alpha_f);
      std::transform(tmp, tmp+16, bbC, bbC, std::plus< std::complex< double> >());

      /* --------------------------------------------- */

      // the following part takes 9 Dirac::Matrixes and translates them to a Tensor (tmp_result)

      QCD::Tensor tmp_result;

      std::complex< double > *data_ptr(NULL);

      size_t index(-1);

      index = 12*Base::col_RED  + Base::col_RED;
      data_ptr = tmp_result.d_data + index;
      data_ptr[  0] =  rrC[ 0];
      data_ptr[  3] =  rrC[ 1];
      data_ptr[  6] =  rrC[ 2];
      data_ptr[  9] =  rrC[ 3];
      data_ptr[ 36] =  rrC[ 4];
      data_ptr[ 39] =  rrC[ 5];
      data_ptr[ 42] =  rrC[ 6];
      data_ptr[ 45] =  rrC[ 7];
      data_ptr[ 72] =  rrC[ 8];
      data_ptr[ 75] =  rrC[ 9];
      data_ptr[ 78] =  rrC[10];
      data_ptr[ 81] =  rrC[11];
      data_ptr[108] =  rrC[12];
      data_ptr[111] =  rrC[13];
      data_ptr[114] =  rrC[14];
      data_ptr[117] =  rrC[15];

      index = 12*Base::col_RED  + Base::col_GREEN;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  rgC[ 0];
      data_ptr[  3] =  rgC[ 1];
      data_ptr[  6] =  rgC[ 2];
      data_ptr[  9] =  rgC[ 3];
      data_ptr[ 36] =  rgC[ 4];
      data_ptr[ 39] =  rgC[ 5];
      data_ptr[ 42] =  rgC[ 6];
      data_ptr[ 45] =  rgC[ 7];
      data_ptr[ 72] =  rgC[ 8];
      data_ptr[ 75] =  rgC[ 9];
      data_ptr[ 78] =  rgC[10];
      data_ptr[ 81] =  rgC[11];
      data_ptr[108] =  rgC[12];
      data_ptr[111] =  rgC[13];
      data_ptr[114] =  rgC[14];
      data_ptr[117] =  rgC[15];

      index = 12*Base::col_RED  + Base::col_BLUE;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  rbC[ 0];
      data_ptr[  3] =  rbC[ 1];
      data_ptr[  6] =  rbC[ 2];
      data_ptr[  9] =  rbC[ 3];
      data_ptr[ 36] =  rbC[ 4];
      data_ptr[ 39] =  rbC[ 5];
      data_ptr[ 42] =  rbC[ 6];
      data_ptr[ 45] =  rbC[ 7];
      data_ptr[ 72] =  rbC[ 8];
      data_ptr[ 75] =  rbC[ 9];
      data_ptr[ 78] =  rbC[10];
      data_ptr[ 81] =  rbC[11];
      data_ptr[108] =  rbC[12];
      data_ptr[111] =  rbC[13];
      data_ptr[114] =  rbC[14];
      data_ptr[117] =  rbC[15];

      /* --------------------------------------------- */

      index = 12*Base::col_GREEN  + Base::col_RED;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  grC[ 0];
      data_ptr[  3] =  grC[ 1];
      data_ptr[  6] =  grC[ 2];
      data_ptr[  9] =  grC[ 3];
      data_ptr[ 36] =  grC[ 4];
      data_ptr[ 39] =  grC[ 5];
      data_ptr[ 42] =  grC[ 6];
      data_ptr[ 45] =  grC[ 7];
      data_ptr[ 72] =  grC[ 8];
      data_ptr[ 75] =  grC[ 9];
      data_ptr[ 78] =  grC[10];
      data_ptr[ 81] =  grC[11];
      data_ptr[108] =  grC[12];
      data_ptr[111] =  grC[13];
      data_ptr[114] =  grC[14];
      data_ptr[117] =  grC[15];

      index = 12*Base::col_GREEN  + Base::col_GREEN;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  ggC[ 0];
      data_ptr[  3] =  ggC[ 1];
      data_ptr[  6] =  ggC[ 2];
      data_ptr[  9] =  ggC[ 3];
      data_ptr[ 36] =  ggC[ 4];
      data_ptr[ 39] =  ggC[ 5];
      data_ptr[ 42] =  ggC[ 6];
      data_ptr[ 45] =  ggC[ 7];
      data_ptr[ 72] =  ggC[ 8];
      data_ptr[ 75] =  ggC[ 9];
      data_ptr[ 78] =  ggC[10];
      data_ptr[ 81] =  ggC[11];
      data_ptr[108] =  ggC[12];
      data_ptr[111] =  ggC[13];
      data_ptr[114] =  ggC[14];
      data_ptr[117] =  ggC[15];

      index = 12*Base::col_GREEN  + Base::col_BLUE;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  gbC[ 0];
      data_ptr[  3] =  gbC[ 1];
      data_ptr[  6] =  gbC[ 2];
      data_ptr[  9] =  gbC[ 3];
      data_ptr[ 36] =  gbC[ 4];
      data_ptr[ 39] =  gbC[ 5];
      data_ptr[ 42] =  gbC[ 6];
      data_ptr[ 45] =  gbC[ 7];
      data_ptr[ 72] =  gbC[ 8];
      data_ptr[ 75] =  gbC[ 9];
      data_ptr[ 78] =  gbC[10];
      data_ptr[ 81] =  gbC[11];
      data_ptr[108] =  gbC[12];
      data_ptr[111] =  gbC[13];
      data_ptr[114] =  gbC[14];
      data_ptr[117] =  gbC[15];

      /* --------------------------------------------- */

      index = 12*Base::col_BLUE  + Base::col_RED;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  brC[ 0];
      data_ptr[  3] =  brC[ 1];
      data_ptr[  6] =  brC[ 2];
      data_ptr[  9] =  brC[ 3];
      data_ptr[ 36] =  brC[ 4];
      data_ptr[ 39] =  brC[ 5];
      data_ptr[ 42] =  brC[ 6];
      data_ptr[ 45] =  brC[ 7];
      data_ptr[ 72] =  brC[ 8];
      data_ptr[ 75] =  brC[ 9];
      data_ptr[ 78] =  brC[10];
      data_ptr[ 81] =  brC[11];
      data_ptr[108] =  brC[12];
      data_ptr[111] =  brC[13];
      data_ptr[114] =  brC[14];
      data_ptr[117] =  brC[15];

      index = 12*Base::col_BLUE  + Base::col_GREEN;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  bgC[ 0];
      data_ptr[  3] =  bgC[ 1];
      data_ptr[  6] =  bgC[ 2];
      data_ptr[  9] =  bgC[ 3];
      data_ptr[ 36] =  bgC[ 4];
      data_ptr[ 39] =  bgC[ 5];
      data_ptr[ 42] =  bgC[ 6];
      data_ptr[ 45] =  bgC[ 7];
      data_ptr[ 72] =  bgC[ 8];
      data_ptr[ 75] =  bgC[ 9];
      data_ptr[ 78] =  bgC[10];
      data_ptr[ 81] =  bgC[11];
      data_ptr[108] =  bgC[12];
      data_ptr[111] =  bgC[13];
      data_ptr[114] =  bgC[14];
      data_ptr[117] =  bgC[15];

      index = 12*Base::col_BLUE  + Base::col_BLUE;
      data_ptr = tmp_result.d_data + index;

      data_ptr[  0] =  bbC[ 0];
      data_ptr[  3] =  bbC[ 1];
      data_ptr[  6] =  bbC[ 2];
      data_ptr[  9] =  bbC[ 3];
      data_ptr[ 36] =  bbC[ 4];
      data_ptr[ 39] =  bbC[ 5];
      data_ptr[ 42] =  bbC[ 6];
      data_ptr[ 45] =  bbC[ 7];
      data_ptr[ 72] =  bbC[ 8];
      data_ptr[ 75] =  bbC[ 9];
      data_ptr[ 78] =  bbC[10];
      data_ptr[ 81] =  bbC[11];
      data_ptr[108] =  bbC[12];
      data_ptr[111] =  bbC[13];
      data_ptr[114] =  bbC[14];
      data_ptr[117] =  bbC[15];

      data_ptr = NULL;

      tmp_result *= alpha_f_factors[idx];
      result += tmp_result;
    }
  }

}

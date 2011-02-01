#include "HYP.ih"

namespace Smear
{
  void HYP::smear(Core::Field< QCD::Gauge > &field) const
  {

    //std::complex< double > const COMPLEX_ZERO(0.0, 0.0);
    double const alpha3_half(d_alpha3/2.0);
    double const alpha2_fourth(d_alpha2/4.0);
    double const alpha1_sixth(d_alpha1/6.0);

    // ********************************************************
    // 1st level of smearing
    // ********************************************************

    Core::Field< QCD::Gauge > tmpField3(field);
    tmpField3 *= 1.0-d_alpha3;

    // in fatLinksMUNU, no links in directions MU nor NU enter
    // since fatLinksMUNU=fatLinksNUMU, we only use MU<NU
    Core::Field< QCD::Gauge > fatLinksXY(tmpField3);
    Core::Field< QCD::Gauge > fatLinksXZ(tmpField3);
    Core::Field< QCD::Gauge > fatLinksXT(tmpField3);
    Core::Field< QCD::Gauge > fatLinksYZ(tmpField3);
    Core::Field< QCD::Gauge > fatLinksYT(tmpField3);
    Core::Field< QCD::Gauge > fatLinksZT(tmpField3);

    {

    // naming convention: munuHYP involves the 2 staples in the mu-nu plane
    // that would be accomplished to a plaquette by the link in mu-direction

    // x-direction

    Core::Field< SU3::Matrix > xyHYP
           = Path::staple(field, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xyHYP += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xyHYP *= alpha3_half;

    Core::Field< SU3::Matrix > xzHYP
           = Path::staple(field, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xzHYP += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xzHYP *= alpha3_half;

    Core::Field< SU3::Matrix > xtHYP
           = Path::staple(field, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xtHYP += Path::staple(field, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xtHYP *= alpha3_half;

    // y-direction

    Core::Field< SU3::Matrix > yxHYP
           = Path::staple(field, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yxHYP += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yxHYP *= alpha3_half;

    Core::Field< SU3::Matrix > yzHYP
           = Path::staple(field, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yzHYP += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yzHYP *= alpha3_half;

    Core::Field< SU3::Matrix > ytHYP
           = Path::staple(field, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    ytHYP += Path::staple(field, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    ytHYP *= alpha3_half;

    // z-direction

    Core::Field< SU3::Matrix > zxHYP
           = Path::staple(field, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zxHYP += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zxHYP *= alpha3_half;

    Core::Field< SU3::Matrix > zyHYP
           = Path::staple(field, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zyHYP += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zyHYP *= alpha3_half;

    Core::Field< SU3::Matrix > ztHYP
           = Path::staple(field, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    ztHYP += Path::staple(field, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    ztHYP *= alpha3_half;

    // t-direction

    Core::Field< SU3::Matrix > txHYP
           = Path::staple(field, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    txHYP += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    txHYP *= alpha3_half;

    Core::Field< SU3::Matrix > tyHYP
           = Path::staple(field, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tyHYP += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tyHYP *= alpha3_half;

    Core::Field< SU3::Matrix > tzHYP
           = Path::staple(field, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tzHYP += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tzHYP *= alpha3_half;

    // now we assign the pieces according to the HYP blocking procedure

    fatLinksXY.component< SU3::Matrix >(Base::idx_Z) += ztHYP;
    fatLinksXZ.component< SU3::Matrix >(Base::idx_Y) += ytHYP;
    fatLinksXT.component< SU3::Matrix >(Base::idx_Y) += yzHYP;
    fatLinksYZ.component< SU3::Matrix >(Base::idx_X) += xtHYP;
    fatLinksYT.component< SU3::Matrix >(Base::idx_X) += xzHYP;
    fatLinksZT.component< SU3::Matrix >(Base::idx_X) += xyHYP;
    fatLinksXY.component< SU3::Matrix >(Base::idx_T) += tzHYP;
    fatLinksXZ.component< SU3::Matrix >(Base::idx_T) += tyHYP;
    fatLinksXT.component< SU3::Matrix >(Base::idx_Z) += zyHYP;
    fatLinksYZ.component< SU3::Matrix >(Base::idx_T) += txHYP;
    fatLinksYT.component< SU3::Matrix >(Base::idx_Z) += zxHYP;
    fatLinksZT.component< SU3::Matrix >(Base::idx_Y) += yxHYP;

    }


    Tool::reunitarize(&fatLinksXY);
    Tool::reunitarize(&fatLinksXZ);
    Tool::reunitarize(&fatLinksXT);
    Tool::reunitarize(&fatLinksXY);
    Tool::reunitarize(&fatLinksYZ);
    Tool::reunitarize(&fatLinksYT);
    Tool::reunitarize(&fatLinksXZ);
    Tool::reunitarize(&fatLinksYZ);
    Tool::reunitarize(&fatLinksZT);
    Tool::reunitarize(&fatLinksXT);
    Tool::reunitarize(&fatLinksYT);
    Tool::reunitarize(&fatLinksZT);


    // ********************************************************
    // 2nd level of smearing: create another set of fat links using the fat links from level 1
    // ********************************************************

    Core::Field< QCD::Gauge > tmpField2(field);
    tmpField2 *= 1.0-d_alpha2;

    // in fatLinksMU, no links in direction MU enter
    Core::Field< QCD::Gauge > fatLinksX(tmpField2);
    Core::Field< QCD::Gauge > fatLinksY(tmpField2);
    Core::Field< QCD::Gauge > fatLinksZ(tmpField2);
    Core::Field< QCD::Gauge > fatLinksT(tmpField2);

    {

    Core::Field< SU3::Matrix > xyHYP
           = Path::HYPstaple(fatLinksXY, fatLinksXZ, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    xyHYP += Path::HYPstaple(fatLinksXY, fatLinksXZ, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    xyHYP += Path::HYPstaple(fatLinksXY, fatLinksXT, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    xyHYP += Path::HYPstaple(fatLinksXY, fatLinksXT, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    xyHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > xzHYP
           = Path::HYPstaple(fatLinksXZ, fatLinksXY, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    xzHYP += Path::HYPstaple(fatLinksXZ, fatLinksXY, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    xzHYP += Path::HYPstaple(fatLinksXZ, fatLinksXT, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    xzHYP += Path::HYPstaple(fatLinksXZ, fatLinksXT, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    xzHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > xtHYP
           = Path::HYPstaple(fatLinksXT, fatLinksXY, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    xtHYP += Path::HYPstaple(fatLinksXT, fatLinksXY, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    xtHYP += Path::HYPstaple(fatLinksXT, fatLinksXZ, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    xtHYP += Path::HYPstaple(fatLinksXT, fatLinksXZ, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    xtHYP *= alpha2_fourth;


    // X component of fatLinksX will not be used
    //fatLinksX.component< SU3::Matrix >(Base::idx_X) *= COMPLEX_ZERO;
    fatLinksX.component< SU3::Matrix >(Base::idx_Y) += xyHYP;
    fatLinksX.component< SU3::Matrix >(Base::idx_Z) += xzHYP;
    fatLinksX.component< SU3::Matrix >(Base::idx_T) += xtHYP;

    Tool::reunitarize(&fatLinksX);

    }

    {

    Core::Field< SU3::Matrix > yxHYP
           = Path::HYPstaple(fatLinksXY, fatLinksYZ, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    yxHYP += Path::HYPstaple(fatLinksXY, fatLinksYZ, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    yxHYP += Path::HYPstaple(fatLinksXY, fatLinksYT, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    yxHYP += Path::HYPstaple(fatLinksXY, fatLinksYT, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    yxHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > yzHYP
           = Path::HYPstaple(fatLinksYZ, fatLinksXY, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    yzHYP += Path::HYPstaple(fatLinksYZ, fatLinksXY, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    yzHYP += Path::HYPstaple(fatLinksYZ, fatLinksYT, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    yzHYP += Path::HYPstaple(fatLinksYZ, fatLinksYT, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    yzHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > ytHYP
           = Path::HYPstaple(fatLinksYT, fatLinksXY, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    ytHYP += Path::HYPstaple(fatLinksYT, fatLinksXY, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    ytHYP += Path::HYPstaple(fatLinksYT, fatLinksYZ, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    ytHYP += Path::HYPstaple(fatLinksYT, fatLinksYZ, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    ytHYP *= alpha2_fourth;


    // Y component of fatLinksY will not be used
    //fatLinksY.component< SU3::Matrix >(Base::idx_Y) *= COMPLEX_ZERO;
    fatLinksY.component< SU3::Matrix >(Base::idx_X) += yxHYP;
    fatLinksY.component< SU3::Matrix >(Base::idx_Z) += yzHYP;
    fatLinksY.component< SU3::Matrix >(Base::idx_T) += ytHYP;

    Tool::reunitarize(&fatLinksY);

    }

    {

    Core::Field< SU3::Matrix > zxHYP
           = Path::HYPstaple(fatLinksXZ, fatLinksYZ, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    zxHYP += Path::HYPstaple(fatLinksXZ, fatLinksYZ, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    zxHYP += Path::HYPstaple(fatLinksXZ, fatLinksZT, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    zxHYP += Path::HYPstaple(fatLinksXZ, fatLinksZT, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    zxHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > zyHYP
           = Path::HYPstaple(fatLinksYZ, fatLinksXZ, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    zyHYP += Path::HYPstaple(fatLinksYZ, fatLinksXZ, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    zyHYP += Path::HYPstaple(fatLinksYZ, fatLinksZT, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    zyHYP += Path::HYPstaple(fatLinksYZ, fatLinksZT, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    zyHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > ztHYP
           = Path::HYPstaple(fatLinksZT, fatLinksXZ, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    ztHYP += Path::HYPstaple(fatLinksZT, fatLinksXZ, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    ztHYP += Path::HYPstaple(fatLinksZT, fatLinksYZ, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    ztHYP += Path::HYPstaple(fatLinksZT, fatLinksYZ, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    ztHYP *= alpha2_fourth;


    // Z component of fatLinksZ will not be used
    //fatLinksZ.component< SU3::Matrix >(Base::idx_Z) *= COMPLEX_ZERO;
    fatLinksZ.component< SU3::Matrix >(Base::idx_X) += zxHYP;
    fatLinksZ.component< SU3::Matrix >(Base::idx_Y) += zyHYP;
    fatLinksZ.component< SU3::Matrix >(Base::idx_T) += ztHYP;

    Tool::reunitarize(&fatLinksZ);

    }

    {

    Core::Field< SU3::Matrix > txHYP
           = Path::HYPstaple(fatLinksXT, fatLinksYT, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    txHYP += Path::HYPstaple(fatLinksXT, fatLinksYT, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    txHYP += Path::HYPstaple(fatLinksXT, fatLinksZT, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    txHYP += Path::HYPstaple(fatLinksXT, fatLinksZT, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    txHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > tyHYP
           = Path::HYPstaple(fatLinksYT, fatLinksXT, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    tyHYP += Path::HYPstaple(fatLinksYT, fatLinksXT, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    tyHYP += Path::HYPstaple(fatLinksYT, fatLinksZT, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    tyHYP += Path::HYPstaple(fatLinksYT, fatLinksZT, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    tyHYP *= alpha2_fourth;

    Core::Field< SU3::Matrix > tzHYP
           = Path::HYPstaple(fatLinksZT, fatLinksXT, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    tzHYP += Path::HYPstaple(fatLinksZT, fatLinksXT, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    tzHYP += Path::HYPstaple(fatLinksZT, fatLinksYT, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    tzHYP += Path::HYPstaple(fatLinksZT, fatLinksYT, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    tzHYP *= alpha2_fourth;


    // T component of fatLinksT will not be used
    //fatLinksT.component< SU3::Matrix >(Base::idx_T) *= COMPLEX_ZERO;
    fatLinksT.component< SU3::Matrix >(Base::idx_X) += txHYP;
    fatLinksT.component< SU3::Matrix >(Base::idx_Y) += tyHYP;
    fatLinksT.component< SU3::Matrix >(Base::idx_Z) += tzHYP;

    Tool::reunitarize(&fatLinksT);

    }

    // ********************************************************
    // 3rd and last level of smearing: simple APE-Smearing using the fat links from level 2
    // ********************************************************

    Core::Field< SU3::Matrix > xHYP
          = Path::HYPstaple(fatLinksX, fatLinksY, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xHYP += Path::HYPstaple(fatLinksX, fatLinksY, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xHYP += Path::HYPstaple(fatLinksX, fatLinksZ, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xHYP += Path::HYPstaple(fatLinksX, fatLinksZ, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xHYP += Path::HYPstaple(fatLinksX, fatLinksT, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xHYP += Path::HYPstaple(fatLinksX, fatLinksT, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xHYP *= alpha1_sixth;

    Core::Field< SU3::Matrix > yHYP
          = Path::HYPstaple(fatLinksY, fatLinksX, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksX, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksZ, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksZ, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksT, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksT, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yHYP *= alpha1_sixth;

    Core::Field< SU3::Matrix > zHYP
          = Path::HYPstaple(fatLinksZ, fatLinksX, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksX, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksY, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksY, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksT, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksT, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zHYP *= alpha1_sixth;

    Core::Field< SU3::Matrix > tHYP
          = Path::HYPstaple(fatLinksT, fatLinksX, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksX, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksY, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksY, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksZ, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksZ, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tHYP *= alpha1_sixth;

    field *= 1.0-d_alpha1;
    field.component< SU3::Matrix >(Base::idx_X) += xHYP;
    field.component< SU3::Matrix >(Base::idx_Y) += yHYP;
    field.component< SU3::Matrix >(Base::idx_Z) += zHYP;
    field.component< SU3::Matrix >(Base::idx_T) += tHYP;

    Tool::reunitarize(&field);
  }
}

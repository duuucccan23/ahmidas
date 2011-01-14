#include "HYP.ih"

namespace Smear
{
  void HYP::smear(Core::Field< QCD::Gauge > &field) const
  {

    // ********************************************************
    // 1st level of smearing
    // ********************************************************

    Core::Field< QCD::Gauge > fatLinksXY(field);
    Core::Field< QCD::Gauge > fatLinksXZ(field);
    Core::Field< QCD::Gauge > fatLinksXT(field);
    Core::Field< QCD::Gauge > fatLinksYX(field);
    Core::Field< QCD::Gauge > fatLinksYZ(field);
    Core::Field< QCD::Gauge > fatLinksYT(field);
    Core::Field< QCD::Gauge > fatLinksZX(field);
    Core::Field< QCD::Gauge > fatLinksZY(field);
    Core::Field< QCD::Gauge > fatLinksZT(field);
    Core::Field< QCD::Gauge > fatLinksTX(field);
    Core::Field< QCD::Gauge > fatLinksTY(field);
    Core::Field< QCD::Gauge > fatLinksTZ(field);

    {

    // x-direction

    Core::Field< SU3::Matrix > xyHYP
           = Path::staple(field, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xyHYP += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xyHYP *= d_alpha3;

    Core::Field< SU3::Matrix > xzHYP
           = Path::staple(field, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xzHYP += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xzHYP *= d_alpha3;

    Core::Field< SU3::Matrix > xtHYP
           = Path::staple(field, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xtHYP += Path::staple(field, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xtHYP *= d_alpha3;

    // y-direction

    Core::Field< SU3::Matrix > yxHYP
           = Path::staple(field, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yxHYP += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yxHYP *= d_alpha3;

    Core::Field< SU3::Matrix > yzHYP
           = Path::staple(field, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yzHYP += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yzHYP *= d_alpha3;

    Core::Field< SU3::Matrix > ytHYP
           = Path::staple(field, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    ytHYP += Path::staple(field, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    ytHYP *= d_alpha3;

    // z-direction

    Core::Field< SU3::Matrix > zxHYP
           = Path::staple(field, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zxHYP += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zxHYP *= d_alpha3;

    Core::Field< SU3::Matrix > zyHYP
           = Path::staple(field, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zyHYP += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zyHYP *= d_alpha3;

    Core::Field< SU3::Matrix > ztHYP
           = Path::staple(field, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    ztHYP += Path::staple(field, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    ztHYP *= d_alpha3;

    // t-direction

    Core::Field< SU3::Matrix > txHYP
           = Path::staple(field, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    txHYP += Path::staple(field, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    txHYP *= d_alpha3;

    Core::Field< SU3::Matrix > tyHYP
           = Path::staple(field, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tyHYP += Path::staple(field, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tyHYP *= d_alpha3;

    Core::Field< SU3::Matrix > tzHYP
           = Path::staple(field, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tzHYP += Path::staple(field, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tzHYP *= d_alpha3;

    // now we assign the pieces according to the HYP blocking procedure

    fatLinksXY.component< SU3::Matrix >(Base::idx_Z) += tzHYP;
    fatLinksXZ.component< SU3::Matrix >(Base::idx_Y) += tyHYP;
    fatLinksXT.component< SU3::Matrix >(Base::idx_Y) += zyHYP;
    fatLinksYX.component< SU3::Matrix >(Base::idx_Z) += tzHYP;
    fatLinksYZ.component< SU3::Matrix >(Base::idx_X) += txHYP;
    fatLinksYT.component< SU3::Matrix >(Base::idx_X) += zxHYP;
    fatLinksZX.component< SU3::Matrix >(Base::idx_Y) += tyHYP;
    fatLinksZY.component< SU3::Matrix >(Base::idx_X) += txHYP;
    fatLinksZT.component< SU3::Matrix >(Base::idx_X) += yxHYP;
    fatLinksTX.component< SU3::Matrix >(Base::idx_Y) += zyHYP;
    fatLinksTY.component< SU3::Matrix >(Base::idx_X) += zxHYP;
    fatLinksTZ.component< SU3::Matrix >(Base::idx_X) += yxHYP;
    fatLinksXY.component< SU3::Matrix >(Base::idx_T) += ztHYP;
    fatLinksXZ.component< SU3::Matrix >(Base::idx_T) += ytHYP;
    fatLinksXT.component< SU3::Matrix >(Base::idx_Z) += yzHYP;
    fatLinksYX.component< SU3::Matrix >(Base::idx_T) += ztHYP;
    fatLinksYZ.component< SU3::Matrix >(Base::idx_T) += xtHYP;
    fatLinksYT.component< SU3::Matrix >(Base::idx_Z) += xzHYP;
    fatLinksZX.component< SU3::Matrix >(Base::idx_T) += ytHYP;
    fatLinksZY.component< SU3::Matrix >(Base::idx_T) += xtHYP;
    fatLinksZT.component< SU3::Matrix >(Base::idx_Y) += xyHYP;
    fatLinksTX.component< SU3::Matrix >(Base::idx_Z) += yzHYP;
    fatLinksTY.component< SU3::Matrix >(Base::idx_Z) += xzHYP;
    fatLinksTZ.component< SU3::Matrix >(Base::idx_Y) += xyHYP;

    }


    Tool::reunitarize(&fatLinksXY);
    Tool::reunitarize(&fatLinksXZ);
    Tool::reunitarize(&fatLinksXT);
    Tool::reunitarize(&fatLinksYX);
    Tool::reunitarize(&fatLinksYZ);
    Tool::reunitarize(&fatLinksYT);
    Tool::reunitarize(&fatLinksZX);
    Tool::reunitarize(&fatLinksZY);
    Tool::reunitarize(&fatLinksZT);
    Tool::reunitarize(&fatLinksTX);
    Tool::reunitarize(&fatLinksTY);
    Tool::reunitarize(&fatLinksTZ);


    // ********************************************************
    // 2nd level of smearing: create another set of fat links using the fat links from level 1
    // ********************************************************

    Core::Field< QCD::Gauge > fatLinksX(field);
    Core::Field< QCD::Gauge > fatLinksY(field);
    Core::Field< QCD::Gauge > fatLinksZ(field);
    Core::Field< QCD::Gauge > fatLinksT(field);

    {

    Core::Field< SU3::Matrix > xyHYP
           = Path::HYPstaple(fatLinksYX, fatLinksZX, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xyHYP += Path::HYPstaple(fatLinksYX, fatLinksZX, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xyHYP += Path::HYPstaple(fatLinksYX, fatLinksTX, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xyHYP += Path::HYPstaple(fatLinksYX, fatLinksTX, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xyHYP *= d_alpha2;

    Core::Field< SU3::Matrix > xzHYP
           = Path::HYPstaple(fatLinksZX, fatLinksYX, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xzHYP += Path::HYPstaple(fatLinksZX, fatLinksYX, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xzHYP += Path::HYPstaple(fatLinksZX, fatLinksTX, Base::idx_T, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xzHYP += Path::HYPstaple(fatLinksZX, fatLinksTX, Base::idx_T, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xzHYP *= d_alpha2;

    Core::Field< SU3::Matrix > xtHYP
           = Path::HYPstaple(fatLinksTX, fatLinksYX, Base::idx_Y, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xtHYP += Path::HYPstaple(fatLinksTX, fatLinksYX, Base::idx_Y, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xtHYP += Path::HYPstaple(fatLinksTX, fatLinksZX, Base::idx_Z, Base::dir_UP,   Base::idx_X, Base::dir_UP);
    xtHYP += Path::HYPstaple(fatLinksTX, fatLinksZX, Base::idx_Z, Base::dir_DOWN, Base::idx_X, Base::dir_UP);
    xtHYP *= d_alpha2;


    // X component of fatLinksX will not be used
    fatLinksX.component< SU3::Matrix >(Base::idx_Y) += xyHYP;
    fatLinksX.component< SU3::Matrix >(Base::idx_Z) += xzHYP;
    fatLinksX.component< SU3::Matrix >(Base::idx_T) += xtHYP;

    Tool::reunitarize(&fatLinksX);

    }

    {

    Core::Field< SU3::Matrix > yxHYP
           = Path::HYPstaple(fatLinksXY, fatLinksZY, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yxHYP += Path::HYPstaple(fatLinksXY, fatLinksZY, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yxHYP += Path::HYPstaple(fatLinksXY, fatLinksTY, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yxHYP += Path::HYPstaple(fatLinksXY, fatLinksTY, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yxHYP *= d_alpha2;

    Core::Field< SU3::Matrix > yzHYP
           = Path::HYPstaple(fatLinksZY, fatLinksXY, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yzHYP += Path::HYPstaple(fatLinksZY, fatLinksXY, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yzHYP += Path::HYPstaple(fatLinksZY, fatLinksTY, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yzHYP += Path::HYPstaple(fatLinksZY, fatLinksTY, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yzHYP *= d_alpha2;

    Core::Field< SU3::Matrix > ytHYP
           = Path::HYPstaple(fatLinksTY, fatLinksXY, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    ytHYP += Path::HYPstaple(fatLinksTY, fatLinksXY, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    ytHYP += Path::HYPstaple(fatLinksTY, fatLinksZY, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    ytHYP += Path::HYPstaple(fatLinksTY, fatLinksZY, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    ytHYP *= d_alpha2;


    // Y component of fatLinksY will not be used
    fatLinksY.component< SU3::Matrix >(Base::idx_X) += yxHYP;
    fatLinksY.component< SU3::Matrix >(Base::idx_Z) += yzHYP;
    fatLinksY.component< SU3::Matrix >(Base::idx_T) += ytHYP;

    Tool::reunitarize(&fatLinksY);

    }

    {

    Core::Field< SU3::Matrix > zxHYP
           = Path::HYPstaple(fatLinksXZ, fatLinksYZ, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zxHYP += Path::HYPstaple(fatLinksXZ, fatLinksYZ, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zxHYP += Path::HYPstaple(fatLinksXZ, fatLinksTZ, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zxHYP += Path::HYPstaple(fatLinksXZ, fatLinksTZ, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zxHYP *= d_alpha2;

    Core::Field< SU3::Matrix > zyHYP
           = Path::HYPstaple(fatLinksYZ, fatLinksXZ, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zyHYP += Path::HYPstaple(fatLinksYZ, fatLinksXZ, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zyHYP += Path::HYPstaple(fatLinksYZ, fatLinksTZ, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zyHYP += Path::HYPstaple(fatLinksYZ, fatLinksTZ, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zyHYP *= d_alpha2;

    Core::Field< SU3::Matrix > ztHYP
           = Path::HYPstaple(fatLinksTZ, fatLinksXZ, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    ztHYP += Path::HYPstaple(fatLinksTZ, fatLinksXZ, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    ztHYP += Path::HYPstaple(fatLinksTZ, fatLinksYZ, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    ztHYP += Path::HYPstaple(fatLinksTZ, fatLinksYZ, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    ztHYP *= d_alpha2;


    // Z component of fatLinksZ will not be used
    fatLinksZ.component< SU3::Matrix >(Base::idx_X) += zxHYP;
    fatLinksZ.component< SU3::Matrix >(Base::idx_Y) += zyHYP;
    fatLinksZ.component< SU3::Matrix >(Base::idx_T) += ztHYP;

    Tool::reunitarize(&fatLinksZ);

    }

    {

    Core::Field< SU3::Matrix > txHYP
           = Path::HYPstaple(fatLinksXT, fatLinksYT, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    txHYP += Path::HYPstaple(fatLinksXT, fatLinksYT, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    txHYP += Path::HYPstaple(fatLinksXT, fatLinksZT, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    txHYP += Path::HYPstaple(fatLinksXT, fatLinksZT, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    txHYP *= d_alpha2;

    Core::Field< SU3::Matrix > tyHYP
           = Path::HYPstaple(fatLinksYT, fatLinksXT, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tyHYP += Path::HYPstaple(fatLinksYT, fatLinksXT, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tyHYP += Path::HYPstaple(fatLinksYT, fatLinksZT, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tyHYP += Path::HYPstaple(fatLinksYT, fatLinksZT, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tyHYP *= d_alpha2;

    Core::Field< SU3::Matrix > tzHYP
           = Path::HYPstaple(fatLinksZT, fatLinksXT, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tzHYP += Path::HYPstaple(fatLinksZT, fatLinksXT, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tzHYP += Path::HYPstaple(fatLinksZT, fatLinksYT, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tzHYP += Path::HYPstaple(fatLinksZT, fatLinksYT, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tzHYP *= d_alpha2;


    // T component of fatLinksT will not be used
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
    xHYP *= d_alpha1;

    Core::Field< SU3::Matrix > yHYP
          = Path::HYPstaple(fatLinksY, fatLinksX, Base::idx_X, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksX, Base::idx_X, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksZ, Base::idx_Z, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksZ, Base::idx_Z, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksT, Base::idx_T, Base::dir_UP,   Base::idx_Y, Base::dir_UP);
    yHYP += Path::HYPstaple(fatLinksY, fatLinksT, Base::idx_T, Base::dir_DOWN, Base::idx_Y, Base::dir_UP);
    yHYP *= d_alpha1;

    Core::Field< SU3::Matrix > zHYP
          = Path::HYPstaple(fatLinksZ, fatLinksX, Base::idx_X, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksX, Base::idx_X, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksY, Base::idx_Y, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksY, Base::idx_Y, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksT, Base::idx_T, Base::dir_UP,   Base::idx_Z, Base::dir_UP);
    zHYP += Path::HYPstaple(fatLinksZ, fatLinksT, Base::idx_T, Base::dir_DOWN, Base::idx_Z, Base::dir_UP);
    zHYP *= d_alpha1;

    Core::Field< SU3::Matrix > tHYP
          = Path::HYPstaple(fatLinksT, fatLinksX, Base::idx_X, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksX, Base::idx_X, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksY, Base::idx_Y, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksY, Base::idx_Y, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksZ, Base::idx_Z, Base::dir_UP,   Base::idx_T, Base::dir_UP);
    tHYP += Path::HYPstaple(fatLinksT, fatLinksZ, Base::idx_Z, Base::dir_DOWN, Base::idx_T, Base::dir_UP);
    tHYP *= d_alpha1;

    field.component< SU3::Matrix >(Base::idx_X) += xHYP;
    field.component< SU3::Matrix >(Base::idx_Y) += yHYP;
    field.component< SU3::Matrix >(Base::idx_Z) += zHYP;
    field.component< SU3::Matrix >(Base::idx_T) += tHYP;

    Tool::reunitarize(&field);
  }
}

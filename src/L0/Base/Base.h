#pragma once

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdint.h>

namespace Base
{
  // Definitions for clarity in refering to indices
  enum SpaceTimeIndex
  {
    idx_X = 0,
    idx_Y = 1,
    idx_Z = 2,
    idx_T = 3
  };

  // Definitions to be used in refering to directions for shifts
  enum Direction
  {
    dir_DOWN = -1,
    dir_UP = 1
  };

  enum ColourIndex
  {
    col_RED   = 0,
    col_GREEN = 1,
    col_BLUE  = 2
  };

  enum DiracIndex
  {
    gam_1 = 0,
    gam_2 = 1,
    gam_3 = 2,
    gam_4 = 3
  };

  // Definitions to be used in the characterization of sources
  enum SourcePolarization
  {
    sou_UNPOLARIZED      = 0,
    sou_PARTLY_POLARIZED = 1,
    sou_FULLY_POLARIZED  = 2
  };

  enum SourceColorState
  {
	  sou_WHITE   = 0,
	  sou_GENERIC = 1,
	  sou_PURE    = 2
  };

  enum SourceSpatialLayout
  {
	  sou_POINT,
	  sou_FULL,
	  sou_WALL
  };

  enum SourceStochasticTypeFlag
  {
	  sou_Z4 =  4, // this is Z(2) x Z(2) = 1/sqrt(2) *std::complex< double >(+/-1 , +/-1)
	  sou_Z2 =  2, // Z(2) = +/-1
	  sou_P1 =  1, // for test cases: all entries +1
	  sou_M1 = -1  // for test cases: all entries -1
  };


  enum weaveOperator
  {
	  wea_SUM,
	  wea_XOR
  };

  //   enum boundaryConditions
  //   {
  //     bc_ANTIPERIODIC_FIXED;
  //     bc_ANTIPERIODIC_UNIFORM;
  //     bc_PERIODIC_FIXED;
  //     bc_PERIODIC_UNIFORM;
  //   }

  // Definitions to be used for baryon contractions
  enum BaryonInterpolatingField
  {
	  bar_PROTON
  };

  enum BaryonPropagatorProjector
  {
	  proj_PARITY_PLUS_STD,
	  proj_PARITY_MINUS_STD,
	  proj_PARITY_PLUS_TM,
	  proj_PARITY_PLUS_TM_STAR,
	  proj_PARITY_MINUS_TM,
	  proj_1_MINUS_TM, // 1/2(1+gamma0)*i*gamma5*gamma1 in twisted basis (for d quark): 1/4(i*gamma5*gamma1-gamma0*gamma1)
	  proj_1_PLUS_TM,  // 1/2(1+gamma0)*i*gamma5*gamma1 in twisted basis (for u quark): 1/4(i*gamma5*gamma1+gamma0*gamma1)
	  proj_2_MINUS_TM, // same as above but with gamma2
	  proj_2_PLUS_TM,
	  proj_3_MINUS_TM, // same as above but with gamma3
	  proj_3_PLUS_TM,
	  proj_NO_PROJECTOR // this one does not do anything
  };

  enum Operator
  {
    op_UNITY = -1,
    op_GAMMA_1 = 1,
    op_GAMMA_2 = 2,
    op_GAMMA_3 = 3,
    op_GAMMA_4 = 4,
    op_GAMMA_5 = 5,
    op_GAMMA_15 = 6,
    op_GAMMA_25 = 7,
    op_GAMMA_35 = 8,
    op_GAMMA_45 = 9,
    op_GAMMA_12 = 10,
    op_GAMMA_13 = 11,
    op_GAMMA_14 = 12,
    op_GAMMA_23 = 13,
    op_GAMMA_24 = 14,
    op_GAMMA_34 = 15,
    op_CONSERVED_GAMMA_4 = 16,
    op_CONSERVED_GAMMA_1 = 17,
    op_CONSERVED_GAMMA_2 = 18,
    op_CONSERVED_GAMMA_3 = 19,
    op_O44 = 32,
    op_O11 = 33,
    op_O22 = 34,
    op_O33 = 35,
    op_O44_with_substraction = 36, // O44 - (O11 + O22 +O33)/3
    op_O41 = 141,
    op_O42 = 142,
    op_O43 = 143,
    op_O14 = 114,
    op_O12 = 112,
    op_O13 = 113,
    op_O24 = 124,
    op_O21 = 121,
    op_O23 = 123,
    op_O34 = 134,
    op_O31 = 131,
    op_O32 = 132
  };

  // definitions are given in the physical basis
  enum HermitianBilinearOperator
  { //                physical basis            |  twisted basis (changes only)
    // -----------------------------------------|-----------------------------------
    op_G_0 = 0,   //    gamma5                  |  i identity tau_3
    op_G_1 = 1,   //    gamma1                  |
    op_G_2 = 2,   //    gamma2                  |
    op_G_3 = 3,   //    gamma3                  |
    op_G_4 = 4,   // -i gamma0 gamma5           |
    op_G_5 = 5,   // -i gamma0 gamma1           |    gamma5 gamma0 gamma1 tau_3
    op_G_6 = 6,   // -i gamma0 gamma2           |    gamma5 gamma0 gamma2 tau_3
    op_G_7 = 7,   // -i gamma0 gamma3           |    gamma5 gamma0 gamma3 tau_3
    op_G_8 = 8,   //    identity                |  i gamma5 tau_3
    op_G_9 = 9,   // -i gamma_5 gamma1          |
    op_G_10 = 10, // -i gamma_5 gamma2          |
    op_G_11 = 11, // -i gamma_5 gamma3          |
    op_G_12 = 12, //    gamma_0                 |
    op_G_13 = 13, // -i gamma_5 gamma_0 gamma_1 |    gamma_0 gamma_1 tau_3
    op_G_14 = 14, // -i gamma_5 gamma_0 gamma_2 |    gamma_0 gamma_2 tau_3
    op_G_15 = 15  // -i gamma_5 gamma_0 gamma_3 |    gamma_0 gamma_3 tau_3
  };

  enum DiracOperator
  {
    Full = -1,
    A = 0,
    H = 1,
    B = 2,
    Bdagger = 3
  };


  extern bool const bigEndian;
}

#include "Base/Base.inlines"

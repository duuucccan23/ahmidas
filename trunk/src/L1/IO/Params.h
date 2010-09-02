#pragma once

namespace IO
{
  struct PropagatorFormat
  {
    int    flavours;
    int    prec;
    int    nx;
    int    ny;
    int    nz;
    int    nt;
    
    PropagatorFormat();
  };

  struct SourceFormat
  {
    size_t colours;
    size_t flavours;
    size_t prec;
    size_t nx;
    size_t ny;
    size_t nz;
    size_t nt;
    size_t spins;
    
    SourceFormat();
  };

  struct XlfInfo
  {
    char   date[64];
    char   package_version[32];

    double beta;
    double c2_rec;
    double epsilonbar;
    double kappa;
    double mu;
    double mubar;
    double plaq;
  };

  struct IldgFormat
  {
    size_t nx;
    size_t ny;
    size_t nz;
    size_t nt;
    size_t prec;
  };

  struct GaugeInfo
  {
    double plaquetteEnergy;
    bool gaugeRead;
    DML_Checksum checksum; // NOTE Change to Scidac!
    char * xlfInfo;
    char * ildg_data_lfn;
  };

  struct PropInfo
  {
    int splitted;
    int format;
    int precision;
    char * basename;
  };

  struct SourceInfo
  {
    int type;
    int splitted;
    int format;
    int precision;
    int t, x, y, z;
    int sample, nstore, ix;
    int no_flavours;
    char * basename;
  };
}

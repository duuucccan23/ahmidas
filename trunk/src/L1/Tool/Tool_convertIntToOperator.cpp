#include "Tool.ih"
#include <L0/Print.h>
#include <iostream>
#include <sstream>

Base::Operator Tool::convertIntToOperator(int const input)
{
  switch(input)
  {
    case Base::op_UNITY:
      return Base::op_UNITY;
      break;
    case Base::op_GAMMA_1:
      return Base::op_GAMMA_1;
      break;
    case Base::op_GAMMA_2:
      return Base::op_GAMMA_2;
      break;
    case Base::op_GAMMA_3:
      return Base::op_GAMMA_3;
      break;
    case Base::op_GAMMA_4:
      return Base::op_GAMMA_4;
      break;
    case Base::op_GAMMA_5:
      return Base::op_GAMMA_5;
      break;
    case Base::op_GAMMA_15:
      return Base::op_GAMMA_15;
      break;
    case Base::op_GAMMA_25:
      return Base::op_GAMMA_25;
      break;
    case Base::op_GAMMA_35:
      return Base::op_GAMMA_35;
      break;
    case Base::op_GAMMA_45:
      return Base::op_GAMMA_45;
      break;
    case Base::op_GAMMA_12:
      return Base::op_GAMMA_12;
      break;
    case Base::op_GAMMA_13:
      return Base::op_GAMMA_13;
      break;
    case Base::op_GAMMA_14:
      return Base::op_GAMMA_14;
      break;
    case Base::op_GAMMA_23:
      return Base::op_GAMMA_23;
      break;
    case Base::op_GAMMA_24:
      return Base::op_GAMMA_24;
      break;
    case Base::op_GAMMA_34:
      return Base::op_GAMMA_34;
      break;
    case Base::op_CONSERVED_GAMMA_4:
      return Base::op_CONSERVED_GAMMA_4;
      break;
    case Base::op_O44:
      return Base::op_O44;
      break;
    case Base::op_O11:
      return Base::op_O11;
      break;
    case Base::op_O22:
      return Base::op_O22;
      break;
    case Base::op_O33:
      return Base::op_O33;
      break;
    default:
    {
      std::ostringstream out;
      out << "Operator corresponding to integer " << input << " does not exist" << std::endl;
      Print(out.str(), std::cerr);
      exit(1);
    }
  }
}

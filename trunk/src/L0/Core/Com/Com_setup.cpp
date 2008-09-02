/* Publicly accessible static singleton creation function */

#include "../Com.h"

namespace Core
{
  Com &Com::setup()
  {
    if (!s_Com)
      s_Com = new Com();
    return *s_Com;
  }
}

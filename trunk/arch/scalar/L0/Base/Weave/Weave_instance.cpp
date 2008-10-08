/* Publicly accessible static singleton creation function */

#include "../Weave.h"

namespace Base
{
  Weave &Weave::instance()
  {
    if (!s_Weave)
      s_Weave = new Weave();
    return *s_Weave;
  }
}

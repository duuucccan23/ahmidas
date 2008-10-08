/* Initialization of the singleton before it is initialized */

#include "../Weave.h"

namespace Base
{
  Weave* Weave::s_Weave = 0;
}

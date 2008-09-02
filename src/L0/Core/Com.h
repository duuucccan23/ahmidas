#ifndef COM_H
#define COM_H

/* MPI Communication handling */

namespace Core
{
  class Com
  {
    static Com * s_Com;
    Com(); // Private constructor: singleton
    public:
      static Com &setup();


  };

}
#endif

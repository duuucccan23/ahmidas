# Introduction #

This page should be regarded as a very early attempt at introducing a set of naming conventions. Be bold in editing this page. Let's try and consistently implement a set of naming conventions that enables easy access to our source code for both ourselves and other users of our code.

# Details #

Several types of source files may be identified. The structuring of the files is quite simple, every namespace should correspond one to one to a folder/directory, and every class should also correspond one to one to a folder/directory. E.g.: expect to find the class `SU3::Matrix` in the folder `SU3/Matrix`.

  * Header files, internal header files are no-brainers: `Class.h` and `Class.ih`

  * Constructors are named after the class, with letters appended for different constructors. E.g.: `Class.cpp` or `Class_a.cpp`, `Class_b.cpp` etc. Constructors should in the header file be ordered from simplest to complex, with simpler being equivalent with requiring less arguments and more primitive arguments. The ordering in the header file dictates the alfabetical appendices. This is source for concern if new constructors turn out to be needed after a class was considered completed. Renaming of files seems inevitable then. The ordering may not need to be too restrictive.

  * Destructors should be as above.

  * Member functions should be named after the function name. In case of overloading, apply the usual ordering and alphabetical appendices. `Class_member.cpp` or `Class_member_a.cpp`, `Class_member_b.cpp`

  * Inline functions can go in the `Class.inlines`. One might think about a `Class.accessors` file as well, since almost all accessors are so trivial they can be safely grouped, while not all inline functions are as trivial.

  * Overloaded operators can go in the `Class.operators`. But this file can be quite big. Pity `Class_>>.cpp` is not allowed.
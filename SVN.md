# Introduction #

We use SVN, because we had trouble with CVS, and for now we are happy with centralized control. Once the code begins to be used in production, that may change, we might want to move to git or something similar. Until then, some best practices and rules will be gathered here.

# Keywords #

We want the executables to tell us their version, always. That is why we use keywords with them. A good example is main.cpp. The keyword of choice is $Id.

To set this keyword on other source files, use:
`svn propset svn:keywords "Id" myfile`
This is the README file for the Unix distribution
of the MERLIN class library source code
-------------------------------------------------
last mod: 22.12.04

1. Files included in the package
--------------------------------
The downloaded file (merlin_unix.tar) contains the following files: 

merlin-unix/README         - this file
merlin_buld/Makefile        - master gmake makefile
merlin-unix/make_depend    - script to build dependency files
merlin-unix/build_makefile - script to constuct makefiles


2. Installation Requirements
----------------------------
The MERLIN source code makes every attempt to adhere to  the ANSI/ISO C++ 
standard. However, MERLIN uses some advanced features of C++ that not all 
compilers support. This is particularly true for the C++ standard library, 
which Merlin makes heavy use of. In principle, any C++ compiler claiming 
to be a standard compiler should compile the library. To date, the 
following platforms/compilers have been tested: 

linux with GCC 3.2

It is known that gcc versions earlier than 3.0 failed to compile the 
library due to stdlib problems. 

GCC version 3.1 has not been tested.

The Makefiles are known to run correctly with GNU's gmake, but other  
variants of make have not been tested. It is therefore recommended  that 
you make sure that gmake is installed on your system. 

The current version of make_depend uses GCC to construct dependency files 
which are compatable with gmake (make). Again, version earlier that 3.0 
will fail. 

3. Configuring and compiling the library
----------------------------------------

1. extract the files in merlin_unix.tar
2. change directory to merlin-unix
3. edit Makefile to suite your preferences
   
   CPP                 change to your favorite compiler
   MERLIN_ROOT_DIR     should point to the fully qualified
                       Merlin source directory
                       e.g. /my_full_path/Merlin
   MERLIN_BUILD_DIR    should point to the fully qualified
                       merlin-unix directory
                       e.g. /my_full_path/merlin-unix
   MERLIN_INSTALL_DIR  the installation directory for the libraries

Note. If you choose not to use merlin-unix as the build directory, you 
should copy all the installation scripts and Makefile to the directory of 
choice, and modify MERLIN_BUILD_DIR accordingly.

4. type gmake install

gmake will now construct the required makefiles in all the Merlin source 
directories and sub-directories, build dependency files and finally compile
and construct the release version of the library. The library is installed
in MERLIN_INSTALL_DIR as libmerlin.a

Makeing the debug library
-------------------------

To build a separate debug library with no optimisation and debug symbols, 
type

  gmake install_debug

This will install libmerlin_debug.a in MERLIN_INSTALL_DIR

Good Luck!
Nick Walker
DESY
nicholas.walker@desy.de


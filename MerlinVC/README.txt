To build the Merlin library just open the Merlin solution file (Merlin.sln) 
in Visual C++ 2008 Express Edition (VC++ 9.0) and do a Build. 
Default is to build the Release version of the static library.

MerlinAll.sln is a solution to build the release and debug Merlin libraries 
and all examples. Open MerlinAll.sln and build the solution. The default is
static linking and all examples are buildt in Release mode. By default the
ROOT example will ne skipped. For this example the ROOT libraries and include files 
are needed and therefore the variable $(ROOTSYS) must be defined according to
your installation before building. You should get as final message:
========== Build: 15 succeeded, 0 failed, 1 up-to-date, 1 skipped ==========

MerlinAll.sln builds some_example.exe files in MerlinExamples/"some_example_subdir" . 
They can be started within their directory. 
The  MerlinAll solution is not meant for code debugging.


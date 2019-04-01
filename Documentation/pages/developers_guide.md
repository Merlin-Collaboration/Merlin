# Developers guide {#dev_guide}

[TOC]

# Code Style {#code_style}

Merlin++ has a large existing codebase, so not all code meets the current guidelines. New code should and old code should gradually be improved.

## General {#general}

* Class names should be fully descriptive. e.g.  SlicedMacroParticleTracker not SMPTracker.

* Comments should not get in the way of understanding the code - which should be mostly self-explanatory, given.

* Use standard C++11, newer features should be avoided as Scientific Linux with GCC 4.8 is still widely used.

* Prefer modern language features where appropriate, e.g. vectors over C style arrays.

* Prefer references to pointers where possible. For example arguments to methods should take references unless it is an optional parameter which will handle a `nullptr`.

* Do not add hard requirements beyond a modern C++ compiler and [CMake](https://cmake.org/) build system. [ROOT](https://root.cern.ch/), [Python](https://www.python.org/), etc must be optional for users.

## Memory {#memory}

* Use C++ smart pointers where possible. See [unique_ptr](http://en.cppreference.com/w/cpp/memory/unique_ptr), [shared_ptr](http://en.cppreference.com/w/cpp/memory/shared_ptr) and [weak_ptr](http://en.cppreference.com/w/cpp/memory/weak_ptr).

* Every use of `new` must have a matching delete, ideally in the same class.

* When using `new` consider if it could be avoided, for example by using a regular stack value or using a c++ container class instead of bare arrays.

## Formatting {#formatting}

* Merlin++ uses the Allman/BSD brace style (braces on separate lines) and tabs for indentation.

* Good formatting is enforced with a git hook that uses the script [uncrustify](http://uncrustify.sourceforge.net/). If you don't have it already, you can install uncrustify from its [website](http://uncrustify.sourceforge.net/) or through your linux package manager. The hook must be manually installed by copying or linking it from the `tools` to the `.git/hooks` directory:

`cp DeveloperTools/Tools/pre-commit .git/hooks`
or
`ln -sr DeveloperTools/Tools/pre-commit .git/hooks`

* Style can be checked, displayed or fixed with the the `style_wrapper.py` script:
`tools/style_wrapper.py check Merlin/AcceleratorComponent.cpp`
`tools/style_wrapper.py diff Merlin/AcceleratorComponent.cpp`
`tools/style_wrapper.py fix Merlin/AcceleratorComponent.cpp`

* If you are using the Eclipse IDE you can install the formatting style from `tools/eclipse_formatting_style.xml`. Note that Eclipse does not format identically to uncrustify, so do not run the auto formatter on existing code. Expect the commit hook to occasionally request changes.

# Version control {#version_control}

* Merlin++ source code is managed with [GIT](https://git-scm.com/) or [GitHub](https://github.com/MERLIN-Collaboration/merlin-cmake).

There is a large amount of documentation at the [Git website](https://www.git-scm.com/doc) and [GitHub](https://help.github.com/).

* The master branch should always be kept in a good state, buildable and passing all tests.

* New code should be developed in a branch in the developers repo, and a pull request when ready.

* Pull requests should be reviewed by another developer.

## Historic Versions {#historic}

* The main Merlin++ repository starts in 2009. To see older commits you can graft in the older history. From in side your Merlin++ git directory run:

`git remote add Merlin-DESY https://github.com/Sam-Tygier/Merlin-DESY.git`
`git fetch Merlin-DESY`
`echo "ba5f0d0a5eec186a2475c11c11833ea2241d8595 a80fe00a7205b87972f34771f3782fbe31d6d94a" >> .git/info/grafts`

## Making Changes {#changes}

If you want to make put some changes you have made to your Merlin++ program into the repository so that other users can benefit from your additions of bug fixes, the procedure is as follows

1. You will need to get a github account (if you're an expert you probably have one already, if not then it's useful for lots of stuff)
2. You need to get this made a member of the Merlin++ organization. Contact the team who will be happy to help.
3. On the github Merlin website, fork the Merlin project. This establishes a version of the code under `https://github.com/Your user name /Merlin`
4. Be absolutely certain that your code is correct and obeys the rules
5. Add the relevant routines to git's list of files to be staged. You can do this explicitly file by file, or wholesale for the entire directory, typically by by `git add`. If doing it wholesale make sure to do it in the source subdirectory, or all the intermediate binary files will get added to the list and mess up the repository. This need only be done once.
6. `git status` is always helpful
7. `git commit -m "useful message"` will commit your changed files to be uploaded. The useful message should briefly explain what you've done and why.
8. You may need to remind git of where the files were downloaded from with `git remote add origin https://github.com/ your username /Merlin`
9. Then `git push origin master` will make the changes in your github repository
10. Go to the github website, log in, go to the Merlin repository, click on `pull requests` and then `new pull request`. This will start the process of moving the changes in your fork into the mainstream. This requires approval by at least one other developer, who will now be alerted by emails. So it will not be instant. 

# Testing {#testing}

* Merlin++ use the ctest framework provided by cmake.

* Make sure the variable `BUILD_TESTING` is set to `ON`, either with the ccmake tool or by running: `cmake -DBUILD_TESTING=ON ..`
* After building run `ctest` 

* Nightly tests can be seen at the [ABP CDash server](http://abp-cdash.web.cern.ch/abp-cdash/index.php?project=MERLIN).

* New code should have tests.

# Documentation {#doc}

* Merlin++ uses [Doxygen](http://doxygen.nl/) for documentation. See the [Doxygen docs](http://www.doxygen.nl/manual/index.html) for examples of how to write them. 

* These comments will be inserted into the API documentation. Use these to explain what your classes do and how and why they do it - good doxygen is more than just machine generated lists of arguments.

* Merlin++ Doxygen comments are made available using the [web](http://www.accelerators.manchester.ac.uk/merlin/doxygen/)

* Specially marked comments will be inserted into the API documentation.

See the [Doxygen docs](http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html) for examples.

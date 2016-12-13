# force an out of source build to stop clutter
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(STATUS "Merlin requires an out-of-source build.")
	message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
	message(STATUS "CMakeCache.txt")
	message(STATUS "CMakeFiles")
	message(STATUS "Once these files are removed, create a separate directory")
	message(STATUS "and run CMake from there")
	message(FATAL_ERROR "in-source build detected")
endif()

#Set the cmake version required
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

# Set the project name and options
PROJECT(MERLIN)

OPTION(ENABLE_EXAMPLES "Build the example gslprograms. Default OFF" OFF)
OPTION(ENABLE_USER_RUNS "Build any user defined programs in the UserSim folder" OFF)
OPTION(BUILD_TESTING "Build the library test programs. Default ON" ON)
OPTION(ENABLE_OPENMP "Use OpenMP where possible. Default OFF" OFF)
OPTION(BUILD_DYNAMIC "Build Merlin as a dynamic library. Default ON" ON)
OPTION(BUILD_STATIC "Build Merlin as a static library. Default OFF" OFF)
OPTION(BUILD_DOCUMENTATION "Build doxygen documentation. Default OFF" OFF)
OPTION(ORBIT_DEBUG "Enable closed orbit debugging. Default OFF" OFF)
OPTION(INSTALL_HEADERS "Install the Merlin headers. Default OFF" OFF)
OPTION(LIBNUMA "Link to libnuma. See utility/CPUFeatures.h Default OFF" OFF)
OPTION(ENABLE_ROOT "Build the Root output example. Default OFF" OFF)

#Set the default build type
if(NOT CMAKE_BUILD_TYPE)
	MESSAGE(STATUS "No build type set. You probably want to run ccmake and pick one.")
endif()

#MPI default build options
SET( CMAKE_CXX_FLAGS_MPI "-DNDEBUG -Wall -O3 -pedantic -DENABLE_MPI -march=native" CACHE STRING
	"Flags used by the C++ compiler during MPI builds."
	FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MPI
	"-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
	"Flags used for linking binaries during MPI builds."
	FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MPI
	"-Wl,--warn-unresolved-symbols,--warn-once -fPIC" CACHE STRING
	"Flags used by the shared libraries linker during MPI builds."
	FORCE )
MARK_AS_ADVANCED(
	CMAKE_CXX_FLAGS_MPI
	#CMAKE_C_FLAGS_MPI
	CMAKE_EXE_LINKER_FLAGS_MPI
	CMAKE_SHARED_LINKER_FLAGS_MPI)

#MPI debug default build options
SET( CMAKE_CXX_FLAGS_MPIDEBUG "-Wall -g3 -ggdb -pedantic -DENABLE_MPI -march=native" CACHE STRING
	"Flags used by the C++ compiler during MPI Debug builds."
	FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MPIDEBUG
	"-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
	"Flags used for linking binaries during MPI Debug builds."
	FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MPIDEBUG
	"-Wl,--warn-unresolved-symbols,--warn-once -fPIC" CACHE STRING
	"Flags used by the shared libraries linker during MPI Debug builds."
	FORCE )
MARK_AS_ADVANCED(
	CMAKE_CXX_FLAGS_MPIDEBUG
	#CMAKE_C_FLAGS_MPIDEBUG
	CMAKE_EXE_LINKER_FLAGS_MPIDEBUG
	CMAKE_SHARED_LINKER_FLAGS_MPIDEBUG )

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs to include the MPI options
SET(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel MPI MPIDebug."
	FORCE)

#Check for MPI
if(CMAKE_BUILD_TYPE STREQUAL MPI OR CMAKE_BUILD_TYPE STREQUAL MPIDEBUG)
	FIND_PACKAGE(MPI REQUIRED)
	if(NOT MPI_CXX_FOUND)
		MESSAGE(FATAL_ERROR "MPI build requested but no MPI libraries found!")
	endif()
	if(MPI_CXX_FOUND)
		message(STATUS "Using MPI C++ compiler: ${MPI_CXX_COMPILER}")
		message(STATUS "Using MPI runtime: ${MPIEXEC}")
#		message(STATUS "Using MPI process number switch: ${MPIEXEC_NUMPROC_FLAG}")
#		message(STATUS "Using MPI pre-flags: ${MPIEXEC_PREFLAGS}")
#		message(STATUS "Using MPI post-flags: ${MPIEXEC_POSTFLAGS}")
		message(STATUS "Using MPI includes: ${MPI_CXX_INCLUDE_PATH}")
		message(STATUS "Using MPI libraries: ${MPI_CXX_LIBRARIES}")
		message(STATUS "Using MPI compile flags: ${MPI_CXX_COMPILE_FLAGS}")
		message(STATUS "Using MPI link flags: ${MPI_CXX_LINK_FLAGS}")

		#We need the MPI headers to be included!
		include_directories(${MPI_CXX_INCLUDE_PATH})
	endif()
endif()

#Set Merlin version info
set(${PROJECT_NAME}_VERSION "4.9.0")
set(${PROJECT_NAME}_VERSION_MAJOR "4")
set(${PROJECT_NAME}_VERSION_MINOR "9")
set(${PROJECT_NAME}_VERSION_PATCH "0")

# Add the source and includes
include_directories(${MERLIN_SOURCE_DIR}/Merlin/)

AUX_SOURCE_DIRECTORY(Merlin sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/ActiveMonitors sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Apertures sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Construction sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/ControlElements sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Frames sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Implementation sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Implementation/BPMChannels sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Miscellaneous sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/StdComponent sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/StdField sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/StdGeometry sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/Supports sources)
AUX_SOURCE_DIRECTORY(Merlin/AcceleratorModel/TrackingInterface sources)
AUX_SOURCE_DIRECTORY(Merlin/BasicTransport sources)

AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics sources)
AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics/CommonUtilities sources)
AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics/ParticleTracking sources)
AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics/ParticleTracking/Integrators sources)
AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics/ParticleTracking/Output sources)
AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics/SMPTracking sources)
AUX_SOURCE_DIRECTORY(Merlin/BeamDynamics/Utilities sources)

AUX_SOURCE_DIRECTORY(Merlin/BeamModel sources)
AUX_SOURCE_DIRECTORY(Merlin/Channels sources)
AUX_SOURCE_DIRECTORY(Merlin/Collimators sources)
AUX_SOURCE_DIRECTORY(Merlin/Collimators/Output sources)
AUX_SOURCE_DIRECTORY(Merlin/Corrections sources)
AUX_SOURCE_DIRECTORY(Merlin/EuclideanGeometry sources)
AUX_SOURCE_DIRECTORY(Merlin/Exception sources)
AUX_SOURCE_DIRECTORY(Merlin/GroundMotionModels sources)
AUX_SOURCE_DIRECTORY(Merlin/IO sources)
AUX_SOURCE_DIRECTORY(Merlin/MADInterface sources)
AUX_SOURCE_DIRECTORY(Merlin/NumericalUtils sources)
AUX_SOURCE_DIRECTORY(Merlin/Random sources)
AUX_SOURCE_DIRECTORY(Merlin/RingDynamics sources)
AUX_SOURCE_DIRECTORY(Merlin/stdext sources)
AUX_SOURCE_DIRECTORY(Merlin/TLAS sources)
AUX_SOURCE_DIRECTORY(Merlin/utility sources)

#Additional Compiler flags
#ADD_DEFINITIONS("-Wall -Wextra -Wno-unused-parameter -Wno-unused-but-set-parameter -Wno-unused-variable -ansi -pedantic")
ADD_DEFINITIONS("-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -pedantic -Wabi -std=c++11")

#Closed orbit debugging
if(ORBIT_DEBUG)
	ADD_DEFINITIONS("-DDEBUG_CLOSED_ORBIT")
endif(ORBIT_DEBUG)

#Check for Root
if(ENABLE_ROOT)
	#/usr/share/root/cmake/FindROOT.cmake
	SET(CMAKE_MODULE_PATH /usr/share/root/cmake/ ${CMAKE_MODULE_PATH})
	find_package(ROOT REQUIRED)
	include_directories(${ROOT_INCLUDE_DIR})
endif(ENABLE_ROOT)

#Check for OpenMP
if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif(ENABLE_OPENMP)

#Enable to build the MerlinExamples folder
if(ENABLE_EXAMPLES)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	add_subdirectory(MerlinExamples)
endif()

#Enable building of user created code in UserSim folder
if(ENABLE_USER_RUNS)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	add_subdirectory(UserSim)
endif()

#Enable to build the MerlinTests folder
if(BUILD_TESTING)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	enable_testing()
	add_subdirectory(MerlinTests)
endif()

#make the libmerlin.so shared library from the sources
if(BUILD_DYNAMIC)
	add_library(merlin SHARED ${sources})
endif()

#make libmerlin.a
if(BUILD_STATIC)
	add_library(merlin STATIC ${sources})
endif()

#We want to link to the mpi libraries if this is an MPI build.
if(CMAKE_BUILD_TYPE STREQUAL MPI OR CMAKE_BUILD_TYPE STREQUAL MPIDEBUG)
	#Add in linker and compile flags
	SET_TARGET_PROPERTIES(merlin PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
	SET_TARGET_PROPERTIES(merlin PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
	#and the lib link
	target_link_libraries(merlin ${MPI_CXX_LIBRARIES})
endif()

#If we have OpenMP enabled, set the compiler flags needed
if(OPENMP_FOUND)
	SET_TARGET_PROPERTIES(merlin PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
	SET_TARGET_PROPERTIES(merlin PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

if(LIBNUMA)
	SET_TARGET_PROPERTIES(merlin PROPERTIES COMPILE_FLAGS "-DLIBNUMA")
	target_link_libraries(merlin numa)
endif(LIBNUMA)

install(TARGETS merlin LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
install(DIRECTORY Merlin/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")

# Build doxygen documentation
# make doxygen
if(BUILD_DOCUMENTATION)
find_package(Doxygen)
if(DOXYGEN_FOUND)
	find_file(DOXYFILE_IN "Merlin.Doxyfile.in"
			PATHS "${MERLIN_SOURCE_DIR}/doxygen"
			NO_DEFAULT_PATH
			DOC "Path to the doxygen configuration template file")
	set(DOXYFILE "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(DOXYFILE_IN DEFAULT_MSG "DOXYFILE_IN")

	set(DOXYFILE "Merlin.Doxyfile")
	set(DOXYFILE_OUTPUT_DIR "doc")
	add_custom_target(doxygen
			COMMAND "${DOXYGEN_EXECUTABLE}"
				"${DOXYFILE}" 
			COMMENT "Writing documentation to ${DOXYFILE_OUTPUT_DIR}...")

	configure_file("${DOXYFILE_IN}" "${DOXYFILE}" @ONLY)

	set_property(DIRECTORY 
			APPEND PROPERTY
			ADDITIONAL_MAKE_CLEAN_FILES
			"${DOXYFILE_OUTPUT_DIR}/${DOXYFILE_HTML_DIR}")

endif()
endif()



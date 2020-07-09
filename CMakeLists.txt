# force an out of source build to stop clutter
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
	message(STATUS "Merlin++ requires an out-of-source build.")
	message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
	message(STATUS "CMakeCache.txt")
	message(STATUS "CMakeFiles")
	message(STATUS "Once these files are removed, create a separate directory")
	message(STATUS "and run CMake from there")
	message(FATAL_ERROR "in-source build detected")
endif()

#Set the cmake version required
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# Set the project name and options
PROJECT(MERLIN)

OPTION(ENABLE_EXAMPLES "Build the example programs. Default ON" ON)
OPTION(ENABLE_DEVTOOLS "Build the developer tools. Default OFF" OFF)
OPTION(ENABLE_USER_RUNS "Build any user defined programs in the UserSim folder" ON)
OPTION(BUILD_TESTING "Build the library test programs. Default ON" ON)
OPTION(ENABLE_OPENMP "Use OpenMP where possible. Default OFF" OFF)
OPTION(ENABLE_MPI "Use MPI where possible. Default OFF" OFF)
if(NOT WIN32)
	OPTION(BUILD_DYNAMIC "Build Merlin++ as a dynamic library. Default ON" ON)
	OPTION(BUILD_STATIC "Build Merlin++ as a static library. Default OFF" OFF)
else()
	OPTION(BUILD_DYNAMIC "Build Merlin++ as a dynamic library. Default OFF" OFF)
	OPTION(BUILD_STATIC "Build Merlin++ as a static library. Default ON" ON)
endif()
OPTION(BUILD_DOCUMENTATION "Build doxygen documentation. Default ON" ON)
OPTION(ORBIT_DEBUG "Enable closed orbit debugging. Default OFF" OFF)
OPTION(INSTALL_HEADERS "Install the Merlin++ headers. Default OFF" OFF)
OPTION(LIBNUMA "Link to libnuma. See utility/CPUFeatures.h Default OFF" OFF)
OPTION(ENABLE_ROOT "Build the Root output example. Default OFF" OFF)
OPTION(COVERAGE "Enable build flags for testing code coverage with gcov (only works with GNU compilers)" OFF)
SET(TEST_TIMEOUT "7200" CACHE STRING "Time allowed per test (seconds)")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wno-unused-parameter -Wno-unused-variable")
	# enables optimizations that do not interfere with debugging introduced in GCC 4.8
	SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wno-unused-parameter -Wno-unused-variable")
	# Hide an unnecessary clang warning
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
	SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
else()
	MESSAGE(WARNING "Unknown compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

#Set the default build type
if(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if(CMAKE_BUILD_TYPE MATCHES MPI)
	MESSAGE(FATAL_ERROR "To enable MPI use the Release or Debug build type and the ENABLE_MPI option")
endif()

if(ENABLE_MPI)
	FIND_PACKAGE(MPI REQUIRED)
	if(NOT MPI_CXX_FOUND)
		MESSAGE(FATAL_ERROR "MPI build requested but no MPI libraries found!")
	endif()

	message(STATUS "Using MPI C++ compiler: ${MPI_CXX_COMPILER}")
	message(STATUS "Using MPI runtime: ${MPIEXEC}")
	message(STATUS "Using MPI includes: ${MPI_CXX_INCLUDE_PATH}")
	message(STATUS "Using MPI libraries: ${MPI_CXX_LIBRARIES}")
	message(STATUS "Using MPI compile flags: ${MPI_CXX_COMPILE_FLAGS}")
	message(STATUS "Using MPI link flags: ${MPI_CXX_LINK_FLAGS}")

	include_directories(${MPI_CXX_INCLUDE_PATH})

	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MPI ")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--warn-unresolved-symbols,--warn-once ")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--warn-unresolved-symbols,--warn-once")
	SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

#Set Merlin version info
set(${PROJECT_NAME}_VERSION "5.03")
set(${PROJECT_NAME}_VERSION_MAJOR "5")
set(${PROJECT_NAME}_VERSION_MINOR "03")
set(${PROJECT_NAME}_VERSION_PATCH "0")

configure_file(Merlin++/MerlinVersion.cpp.in ${CMAKE_BINARY_DIR}/Merlin++/MerlinVersion.cpp)
AUX_SOURCE_DIRECTORY(${CMAKE_BINARY_DIR}/Merlin++ sources)

set(GITREV_FILE ${CMAKE_BINARY_DIR}/Merlin++/git_rev.h)
add_custom_command(
	OUTPUT gitrev_file
	COMMAND ${CMAKE_COMMAND} -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR} -DMERLIN_SOURCE_DIR=${MERLIN_SOURCE_DIR} -P ${MERLIN_SOURCE_DIR}/Merlin++/MerlinVersion.cmake
	VERBATIM
)
add_custom_target(gitrev ALL DEPENDS gitrev_file)


# Add the source and includes
include_directories(${MERLIN_SOURCE_DIR}/Merlin++/)
AUX_SOURCE_DIRECTORY(Merlin++ sources)

#Closed orbit debugging
if(ORBIT_DEBUG)
	ADD_DEFINITIONS("-DDEBUG_CLOSED_ORBIT")
endif(ORBIT_DEBUG)

#Check for Root
if(ENABLE_ROOT)
	#/usr/share/root/cmake/FindROOT.cmake
	SET(CMAKE_MODULE_PATH /usr/share/root/cmake/ ${CMAKE_MODULE_PATH})
	find_package(ROOT REQUIRED)
	include_directories(${ROOT_INCLUDE_DIRS})
endif(ENABLE_ROOT)

#Check for and set up OpenMP
if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
	if(NOT OPENMP_FOUND)
		MESSAGE(FATAL_ERROR "OpenMP build requested but no OpenMP libraries found!")
	endif()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(ENABLE_OPENMP)

#Enable to build the ExamplesTutorials folder
if(ENABLE_EXAMPLES)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	add_subdirectory(ExamplesTutorials)
endif()

#Enable to build the DeveloperTools/Tools folder
if(ENABLE_DEVTOOLS)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	add_subdirectory(DeveloperTools/Tools/multitrack)
endif()

#Enable building of user created code in UserSim folder
if(ENABLE_USER_RUNS)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	add_subdirectory(UserSim)
endif()

#Enable to build the DeveloperTools/CodeTesting folder
if(BUILD_TESTING)
	set(MERLIN_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
	enable_testing()
	add_subdirectory(DeveloperTools/CodeTesting)
endif()

#make the libmerlin.so shared library from the sources
if(BUILD_DYNAMIC)
	add_library(merlin++ SHARED ${sources})
endif()

#make libmerlin.a
if(BUILD_STATIC)
	add_library(merlin++ STATIC ${sources})
endif()

add_dependencies(merlin++ gitrev)

if(ENABLE_MPI)
	target_link_libraries(merlin++ ${MPI_CXX_LIBRARIES})
endif()


IF(COVERAGE)
	set(COVERAGE_FLAGS "-fprofile-arcs")
	target_compile_options(merlin++ PRIVATE -fprofile-arcs -ftest-coverage)
	SET_TARGET_PROPERTIES(merlin++ PROPERTIES LINK_FLAGS " ${COVERAGE_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
ENDIF(COVERAGE)

if(LIBNUMA)
	SET_TARGET_PROPERTIES(merlin++ PROPERTIES COMPILE_FLAGS "-DLIBNUMA")
	target_link_libraries(merlin++ numa)
endif(LIBNUMA)

install(TARGETS merlin++ LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
install(DIRECTORY Merlin++/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include FILES_MATCHING PATTERN "*.h")

# Build doxygen documentation
# make doxygen

find_package(Doxygen)

if(NOT DOXYGEN_FOUND AND BUILD_DOCUMENTATION)
	MESSAGE(WARNING "Doxygen not found found! Disabling")
	set(BUILD_DOCUMENTATION OFF)
	add_custom_target(doxygen COMMAND ${CMAKE_COMMAND} -E echo Please install doxygen and rerun cmake)
endif()

if(BUILD_DOCUMENTATION)
	find_file(DOXYFILE_IN "Merlin.Doxyfile.in"
			PATHS "${MERLIN_SOURCE_DIR}/Documentation"
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



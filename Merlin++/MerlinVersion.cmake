set(GITREV_FILE ${CMAKE_BINARY_DIR}/Merlin++/git_rev.h)

find_package(Git QUIET)
if(GIT_FOUND)
	execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always --abbrev=40 --dirty=+
	OUTPUT_VARIABLE MERLIN_GIT_VERSION WORKING_DIRECTORY ${MERLIN_SOURCE_DIR})
	string(REPLACE "\n" "" MERLIN_GIT_VERSION ${MERLIN_GIT_VERSION})
else()
	set(MERLIN_GIT_VERSION "")
endif()

file(WRITE ${GITREV_FILE}.tmp "const std::string GIT_REV=\"${MERLIN_GIT_VERSION}\";\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GITREV_FILE}.tmp ${GITREV_FILE} )
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_FILE}.tmp )

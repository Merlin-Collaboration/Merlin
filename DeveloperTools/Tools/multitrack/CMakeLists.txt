macro (merlin_py script_name)
	add_custom_target(${script_name} ALL)
	add_custom_command(TARGET ${script_name}
	                   DEPENDS ${script_name}
	                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${script_name} ${CMAKE_CURRENT_BINARY_DIR}/${script_name})
endmacro (merlin_py)


add_executable(merlin_track merlin_track.cpp)
target_link_libraries(merlin_track merlin++)

file(GLOB pyfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")

#message(FATAL_ERROR ${pyfiles})

foreach(script ${pyfiles})
	merlin_py(${script})
endforeach()

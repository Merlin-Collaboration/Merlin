# macro to install a test executable (need add_test() to call it)
macro (merlin_test dir exe_name)
	set(source_files "")
	foreach(arg IN ITEMS ${ARGN})
		list(APPEND source_files ${dir}/${arg})
	endforeach()
	add_executable(${exe_name} ${source_files})

	# To keep the directory structure, we must set the output_name, add a rule to make the directory
	set_target_properties( ${exe_name} PROPERTIES OUTPUT_NAME ${dir}/${exe_name} )
	set_target_properties( ${exe_name} PROPERTIES COMPILE_FLAGS -g )
	add_custom_command(TARGET ${exe_name} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${dir})
	target_link_libraries(${exe_name} merlin)

	IF(COVERAGE)
		target_compile_options(${exe_name} PRIVATE -fprofile-arcs -ftest-coverage)
		SET_TARGET_PROPERTIES(${exe_name} PROPERTIES LINK_FLAGS "-fprofile-arcs ${CMAKE_EXE_LINKER_FLAGS}")
	ENDIF(COVERAGE)

endmacro (merlin_test)

# macro to install a test python script (need add_test() to call it)
macro (merlin_test_py dir script_name)
	add_custom_target(${script_name} ALL)
	add_custom_command(TARGET ${script_name}
	                   DEPENDS ${dir}/${script_name}
	                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${script_name} ${CMAKE_CURRENT_BINARY_DIR}/${dir})
endmacro (merlin_test_py)


add_custom_target(TestDataFiles ALL)
add_custom_command(TARGET TestDataFiles
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data)
add_custom_command(TARGET TestDataFiles
                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/outputs)

merlin_test_py(BasicTests have_python.py)
add_test(have_python BasicTests/have_python.py)

merlin_test(BasicTests bunch_test bunch_test.cpp)
add_test(bunch_test BasicTests/bunch_test)

merlin_test(BasicTests landau_test landau_test.cpp)
add_test(landau_test BasicTests/landau_test)

merlin_test(BasicTests aperture_test aperture_test.cpp)
add_test(aperture_test BasicTests/aperture_test)

merlin_test(BasicTests collimate_particle_process_test collimate_particle_process_test.cpp)
add_test(collimate_particle_process_test BasicTests/collimate_particle_process_test)

merlin_test(OpticsTests lhc_optics_test lhc_optics_test.cpp)
add_test(lhc_optics_test OpticsTests/lhc_optics_test)

merlin_test(OpticsTests lhc_fft_tune_test lhc_fft_tune_test.cpp)
add_test(lhc_fft_tune_test OpticsTests/lhc_fft_tune_test)

merlin_test(ScatteringTests cu50_test cu50_test.cpp)
add_test(cu50_test ScatteringTests/cu50_test)

merlin_test(ScatteringTests lhc_collimation_test lhc_collimation_test.cpp)
add_test(lhc_collimation_test ScatteringTests/lhc_collimation_test)
merlin_test(ScatteringTests rand_test rand_test.cpp)

merlin_test(ScatteringTests cu50_sixtrack_test cu50_sixtrack_test.cpp)
add_test(cu50_sixtrack_test ScatteringTests/cu50_sixtrack_test)

find_program(MEMORYCHECK_COMMAND NAMES valgrind)
INCLUDE(Dart)



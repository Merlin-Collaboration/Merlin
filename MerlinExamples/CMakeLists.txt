macro (merlin_example dir exe_name)
	set(source_files "")
	foreach(arg IN ITEMS ${ARGN})
		list(APPEND source_files ${dir}/${arg})
	endforeach()
	add_executable(${exe_name} ${source_files})

	# To keep the directory structure, we must set the output_name, add a rule to make the directory
	set_target_properties( ${exe_name} PROPERTIES OUTPUT_NAME ${dir}/${exe_name} )
	add_custom_command(TARGET ${exe_name} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${dir})
	target_link_libraries(${exe_name} merlin)
endmacro (merlin_example)

# macro to install a test python script (need add_test() to call it)
macro (merlin_example_py dir script_name)
	add_custom_target(${script_name} ALL)
	add_custom_command(TARGET ${script_name}
	                   DEPENDS ${dir}/${script_name}
	                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${script_name} ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${script_name})
endmacro (merlin_example_py)

add_custom_target(ExampleDataFiles ALL)
add_custom_command(TARGET ExampleDataFiles
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/lattices ${CMAKE_CURRENT_BINARY_DIR}/lattices)

merlin_example(AcceleratorErrors AccErr main.cpp)
merlin_example(ClosedOrbit ClosedOrbit ClosedOrbit_Main.cpp)
merlin_example(HaloCollimation HaloCollimation HaloCollimation_main.cpp ConstructModel.cpp HaloTracker.cpp QuadIntegrator.cpp)
merlin_example(ILCDFS ILCDFS src/Accelerator.cpp src/AcceleratorWithErrors.cpp src/ConstantGradientAdjustment.cpp src/ConstructSrot.cpp src/DFSApp.cpp src/DFSCorrection.cpp src/DFSOutput.cpp src/EnergyAdjustmentPolicy.cpp src/ILCDFS_IO.cpp src/ILCDFS_main.cpp src/KSEnergyAdjustment.cpp src/ModelConstruction.cpp src/OneToOneCorrection.cpp src/OptimumDispersionMatch.cpp src/ParticleTrackingModel.cpp src/ResponseMatrixGenerator.cpp src/SMPTrackingModel.cpp src/TeslaWakePotential.cpp src/XTFFInterface_1.cpp)
merlin_example(ILCML ILCML ILCML_main.cpp model_construction.cpp TrackingOutput.cpp TeslaWakePotential.cpp QuadReferenceOutput.cpp)
merlin_example(LatticeFunctions LatticeFunctions LatticeFunctions_Main.cpp)
#merlin_example(LHC lhc lhc.cpp) # does not build
#merlin_example(ModelConstruction ConstructModel ConstructModel.cpp) # does not build
#merlin_example(ROOT_trackingOutput trackingOutput MerlinROOT_main.cpp TeslaWakePotential.cpp TrackingOutputROOT.cpp) # Needs work to build with ROOT
merlin_example(SpinTracking SpinTracking SpinTracking_Main.cpp)
merlin_example(StrayFieldClosedOrbit StrayFieldClosedOrbit StrayFieldClosedOrbit_Main.cpp ConstantStrayFieldProcess.cpp)
merlin_example(Tracking Tracking Tracking_Main.cpp)
#merlin_example_py(Tracking Tracking_Main.py)
merlin_example(Trajectory Trajectory Trajectory_Main.cpp)
#merlin_example_py(Trajectory Trajectory_Main.py)
merlin_example(Wakefields wakefile_ex1 ex1.cpp TeslaWakePotential.cpp)
merlin_example(Wakefields wakefile_ex2 ex2.cpp TeslaWakePotential.cpp)
merlin_example(Wakefields wakefile_ex3 ex3.cpp TeslaWakePotential.cpp)



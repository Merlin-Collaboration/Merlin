macro (merlin_example dir exe_name)
	set(source_files "")
	foreach(arg IN ITEMS ${ARGN})
		list(APPEND source_files ${dir}/${arg})
	endforeach()
	add_executable(${exe_name} ${source_files})

	# To keep the directory structure, we must set the output_name, add a rule to make the directory
	set_target_properties( ${exe_name} PROPERTIES OUTPUT_NAME ${dir}/${exe_name} )
	add_custom_command(TARGET ${exe_name} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${dir})
	target_link_libraries(${exe_name} merlin++)
endmacro (merlin_example)

macro (merlin_tutorial dir exe_name)
	set(source_files "")
	foreach(arg IN ITEMS ${ARGN})
		list(APPEND source_files ${dir}/${arg})
	endforeach()
	add_executable(${exe_name} ${source_files})

	# To keep the directory structure, we must set the output_name, add a rule to make the directory
	set_target_properties( ${exe_name} PROPERTIES OUTPUT_NAME ${dir}/${exe_name} )
	add_custom_command(TARGET ${exe_name} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${dir})
	target_link_libraries(${exe_name} merlin++)
endmacro (merlin_tutorial)

# macro to install a test python script (need add_test() to call it)
macro (merlin_examples_py dir script_name)
	add_custom_target(${script_name} ALL)
	add_custom_command(TARGET ${script_name}
	                   DEPENDS ${dir}/Examples/${script_name}
	                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${script_name} ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${script_name})
endmacro (merlin_examples_py)

# macro to install a test python script (need add_test() to call it)
macro (merlin_tutorial_py dir script_name)
	add_custom_target(${script_name} ALL)
	add_custom_command(TARGET ${script_name}
	                   DEPENDS ${dir}/Tutorials/${script_name}
	                   COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${script_name} ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${script_name})
endmacro (merlin_tutorial_py)


add_custom_target(ExampleDataFiles ALL)
add_custom_command(TARGET ExampleDataFiles
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/Examples/lattices ${CMAKE_CURRENT_BINARY_DIR}/Examples/lattices)

add_custom_target(TutorialDataFiles ALL)
add_custom_command(TARGET TutorialDataFiles
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/Tutorials/input ${CMAKE_CURRENT_BINARY_DIR}/Tutorials/input)

merlin_example(Examples/Wakefields ex1 ex1.cpp)
if(ENABLE_ROOT)
merlin_example(Examples/CollimatorMaterials CollimatorMaterials CollimatorMaterials.cpp)
target_link_libraries(CollimatorMaterials ${ROOT_LIBRARIES})
endif()
merlin_example(Examples/AcceleratorErrors AcceleratorErrors main.cpp)
merlin_example(Examples/ClosedOrbit ClosedOrbit ClosedOrbit_Main.cpp)
merlin_example(Examples/HaloCollimation HaloCollimation HaloCollimation_main.cpp ConstructModel.cpp HaloTracker.cpp QuadIntegrator.cpp)
merlin_example(Examples/ILCDFS ILCDFS src/Accelerator.cpp src/AcceleratorWithErrors.cpp src/ConstantGradientAdjustment.cpp src/ConstructSrot.cpp src/DFSApp.cpp src/DFSCorrection.cpp src/DFSOutput.cpp src/EnergyAdjustmentPolicy.cpp src/ILCDFS_IO.cpp src/ILCDFS_main.cpp src/KSEnergyAdjustment.cpp src/ModelConstruction.cpp src/OneToOneCorrection.cpp src/OptimumDispersionMatch.cpp src/ParticleTrackingModel.cpp src/ResponseMatrixGenerator.cpp src/SMPTrackingModel.cpp src/TeslaWakePotential.cpp src/XTFFInterface_1.cpp)
merlin_example(Examples/ILCML ILCML ILCML_main.cpp model_construction.cpp TrackingOutput.cpp TeslaWakePotential.cpp QuadReferenceOutput.cpp)
merlin_example(Examples/LatticeFunctions LatticeFunctions LatticeFunctions_Main.cpp)
#merlin_example(Examples/LHC lhc lhc.cpp) # depreciated
#merlin_example(Examples/ModelConstruction ConstructModel ConstructModel.cpp) # depreciated
#merlin_example(Examples/ROOT_trackingOutput trackingOutput MerlinROOT_main.cpp TeslaWakePotential.cpp TrackingOutputROOT.cpp) # Requires ROOT
merlin_example(Examples/SpinTracking SpinTracking SpinTracking_Main.cpp)
merlin_example(Examples/StrayFieldClosedOrbit StrayFieldClosedOrbit StrayFieldClosedOrbit_Main.cpp ConstantStrayFieldProcess.cpp)
merlin_example(Examples/Tracking Tracking Tracking_Main.cpp)
merlin_example(Examples/Trajectory Trajectory Trajectory_Main.cpp)
merlin_example(Examples/Wakefields wakefile_ex1 ex1.cpp TeslaWakePotential.cpp)
merlin_example(Examples/Wakefields wakefile_ex2 ex2.cpp TeslaWakePotential.cpp)
merlin_example(Examples/Wakefields wakefile_ex3 ex3.cpp TeslaWakePotential.cpp)

merlin_tutorial(Tutorials/1-LatticeConfiguration LatticeConfiguration LatticeConfiguration.cpp)
merlin_tutorial_py(Tutorials/1-LatticeConfiguration LatticeConfiguration.py)

merlin_tutorial(Tutorials/2-LatticeConfigurationMAD LatticeConfigurationMAD LatticeConfigurationMAD.cpp)
merlin_tutorial_py(Tutorials/2-LatticeConfigurationMAD LatticeConfigurationMAD.py)

merlin_tutorial(Tutorials/3-LatticeManipulation LatticeManipulation LatticeManipulation.cpp)
merlin_tutorial_py(Tutorials/3-LatticeManipulation LatticeManipulation.py)

merlin_tutorial(Tutorials/4-ParticleTracking ParticleTracking ParticleTracking.cpp)
merlin_tutorial_py(Tutorials/4-ParticleTracking ParticleTracking.py)

merlin_tutorial(Tutorials/5-ParticleBunchTracking ParticleBunchTracking ParticleBunchTracking.cpp)
merlin_tutorial_py(Tutorials/5-ParticleBunchTracking ParticleBunchTracking.py)

merlin_tutorial(Tutorials/6-LHCLattice LHCLattice LHCLattice.cpp)
merlin_tutorial_py(Tutorials/6-LHCLattice LHCLattice.py)

merlin_tutorial(Tutorials/7-CollimationAndScattering CollimationAndScattering CollimationAndScattering.cpp)
merlin_tutorial_py(Tutorials/7-CollimationAndScattering CollimationAndScattering.py)

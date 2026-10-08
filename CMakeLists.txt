cmake_minimum_required(VERSION 3.19)
project(SiliconeLab VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Versions come from the same record as the byte-for-byte source imports.
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/upstream.json" upstream_record)
foreach(system oil network coating)
    string(JSON ${system}_commit GET "${upstream_record}" sources ${system} integrated_commit)
endforeach()
configure_file(src/source_versions.hpp.in generated/source_versions.hpp @ONLY)

add_executable(siliconelab_generator
    src/generator.cpp src/backends/oil.cpp src/backends/elastomer.cpp
    src/backends/coating_v22.cpp src/backends/coating_v35.cpp)
target_include_directories(siliconelab_generator PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
install(TARGETS siliconelab_generator RUNTIME DESTINATION bin)

add_executable(siliconelab_static src/analysis/static.cpp)
add_executable(siliconelab_topology src/analysis/topology.cpp)
install(TARGETS siliconelab_static siliconelab_topology RUNTIME DESTINATION bin)
configure_file(src/siliconelab_analyzer.py siliconelab_analyzer COPYONLY)
install(PROGRAMS src/siliconelab_analyzer.py DESTINATION bin RENAME siliconelab_analyzer)

add_executable(siliconelab_spatial_dynamics src/analysis/spatial_dynamics.cpp)
add_executable(siliconelab_profiles src/analysis/profiles.cpp)
add_executable(siliconelab_phase vendor/Silicone_Coating/Analysis/phase_separation_analyzer.cpp)
install(TARGETS siliconelab_spatial_dynamics siliconelab_profiles siliconelab_phase RUNTIME DESTINATION bin)

foreach(name part2 part3)
    configure_file(src/siliconelab_${name}.py siliconelab_${name} COPYONLY)
    install(PROGRAMS src/siliconelab_${name}.py DESTINATION bin RENAME siliconelab_${name})
endforeach()
configure_file(src/siliconelab_part2.py siliconelab_part2.py COPYONLY)
install(FILES src/siliconelab_part2.py DESTINATION bin)
configure_file(src/siliconelab_analysis_common.py siliconelab_analysis_common.py COPYONLY)
install(FILES src/siliconelab_analysis_common.py DESTINATION bin)
install(FILES vendor/Silicone_Oil/Analysis/atsc4i.py DESTINATION bin RENAME siliconelab_atsc4i.py)

include(CTest)
if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(NAME film_workflow_contract
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_film_workflow.py")
    add_test(NAME extended_analyzer_contract
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_extended_analyzers.py"
        --generator $<TARGET_FILE:siliconelab_generator>
        --backend $<TARGET_FILE:siliconelab_spatial_dynamics>)
    # Independently built native programs provide a scientific parity reference.
    add_executable(reference_oil vendor/Silicone_Oil/Generator/oil_generator.cpp)
    add_executable(reference_elastomer vendor/PDMS_Elastomer/Generator/pdms_elastomer_generator.cpp)
    add_executable(reference_coating_v22 vendor/Silicone_Coating/V22/v22_generator.cpp)
    add_executable(reference_coating_v35 vendor/Silicone_Coating/V35/v35_generator.cpp)
    add_executable(reference_basic_network vendor/PDMS_Elastomer/Analysis/basic_network_analyzer.cpp)
    add_test(NAME analyzer_contract
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_analyzer.py"
        --generator $<TARGET_FILE:siliconelab_generator>
        --driver "${CMAKE_CURRENT_SOURCE_DIR}/src/siliconelab_analyzer.py"
        --backend-dir $<TARGET_FILE_DIR:siliconelab_static>
        --reference-basic $<TARGET_FILE:reference_basic_network>)
    add_test(NAME generator_contract
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_generator.py"
        --generator $<TARGET_FILE:siliconelab_generator>
        --reference-dir $<TARGET_FILE_DIR:reference_oil>
        --source-dir "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

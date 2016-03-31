cmake_minimum_required(VERSION 3.1)
add_executable(testUpdateFlow testUpdateFlow.cpp)
target_link_libraries(testUpdateFlow PUBLIC multiStateCommon)

if(WIN32)
	set_target_properties(testUpdateFlow PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

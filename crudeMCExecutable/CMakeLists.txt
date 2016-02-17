cmake_minimum_required(VERSION 3.1)
add_executable(crudeMCExecutable crudeMCExecutable.cpp)
target_link_libraries(crudeMCExecutable PUBLIC multiStateCommon)

if(WIN32)
	set_target_properties(crudeMCExecutable PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

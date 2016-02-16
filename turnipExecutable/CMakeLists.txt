cmake_minimum_required(VERSION 3.1)
add_executable(turnipExecutable turnipExecutable.cpp)
target_link_libraries(turnipExecutable PUBLIC multiStateCommon)

if(WIN32)
	set_target_properties(turnipExecutable PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

cmake_minimum_required(VERSION 3.1)
add_executable(turnip turnip.cpp)
target_link_libraries(turnip PUBLIC multiStateCommon)

if(WIN32)
	set_target_properties(turnip PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} PDB_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

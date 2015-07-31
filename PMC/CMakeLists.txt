cmake_minimum_required(VERSION 3.1)
add_executable(PMC PMC.cpp)
target_link_libraries(PMC PUBLIC multiStateCommon)

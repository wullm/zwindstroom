cmake_minimum_required(VERSION 3.11)
project(class)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(GSL REQUIRED)
mark_as_advanced(pkgcfg_lib_GSL_gsl pkgcfg_lib_GSL_gslcblas pkgcfg_lib_GSL_m)

file(GLOB 3FA_SOURCES src/*.c)

add_library(lib3fa ${3FA_SOURCES})
add_library(3fa::lib3fa ALIAS lib3fa)
set_target_properties(lib3fa PROPERTIES OUTPUT_NAME 3fa)
set_target_properties(lib3fa PROPERTIES C_STANDARD 11)
set_property(TARGET lib3fa PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(lib3fa PRIVATE m)
target_link_libraries(lib3fa PRIVATE GSL::gsl)
target_include_directories(lib3fa PUBLIC include)

add_executable(test_3fa tests/test.c)
target_link_libraries(test_3fa PRIVATE 3fa::lib3fa)

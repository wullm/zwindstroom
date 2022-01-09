cmake_minimum_required(VERSION 3.11)
project(3fa VERSION 1.0.0 DESCRIPTION "Library for cosmic neutrino calculations")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(GSL REQUIRED)
mark_as_advanced(pkgcfg_lib_GSL_gsl pkgcfg_lib_GSL_gslcblas pkgcfg_lib_GSL_m)

file(GLOB 3FA_SOURCES src/*.c)

# Static library
add_library(3fa_static STATIC ${3FA_SOURCES})
add_library(3fa::3fa_static ALIAS 3fa_static)
set_target_properties(3fa_static PROPERTIES OUTPUT_NAME 3fa)
set_target_properties(3fa_static PROPERTIES C_STANDARD 11)
set_property(TARGET 3fa_static PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(3fa_static PRIVATE m)
target_link_libraries(3fa_static PRIVATE GSL::gsl)
target_include_directories(3fa_static PUBLIC include)

# Shared library
add_library(3fa_shared SHARED ${3FA_SOURCES})
add_library(3fa::3fa_shared ALIAS 3fa_shared)
set_target_properties(3fa_shared PROPERTIES OUTPUT_NAME 3fa)
set_target_properties(3fa_shared PROPERTIES C_STANDARD 11)
set_property(TARGET 3fa_shared PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(3fa_shared PRIVATE m)
target_link_libraries(3fa_shared PRIVATE GSL::gsl)
target_include_directories(3fa_shared PUBLIC include)

add_executable(test_3fa tests/test.c)
target_link_libraries(test_3fa PRIVATE 3fa::3fa_static)

cmake_minimum_required(VERSION 3.15)
project(SecChem LANGUAGES CXX)

add_library(SecChem)
find_package(Eigen3 5 REQUIRED)
target_link_libraries(SecChem PUBLIC Eigen3::Eigen)
target_sources(SecChem PRIVATE SecChem/Molecule.cpp)

add_executable(test main.cpp)
target_link_libraries(test PRIVATE SecChem)

add_subdirectory(tests)

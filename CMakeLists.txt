cmake_minimum_required(VERSION 3.15)
project(SecChem LANGUAGES CXX)

if (NOT DEFINED CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 17)
    set(CMAKE_CXX_STANDARD 17)
endif ()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(SecChem INTERFACE)
add_library(SecChem::SecChem ALIAS SecChem)
include(GNUInstallDirs)
target_include_directories(SecChem INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

#-----------------------------------------------------------------------------------

if (TARGET Eigen3::Eigen)
    message(STATUS "SecChem will use the Eigen3 target defined by ancestor project")
else ()
    if (NOT EIGEN3_SOURCE AND NOT EIGEN3_REPO AND NOT Eigen3_FOUND)
        find_package(Eigen3 5 NO_MODULE)
        if (Eigen3_FOUND)
            message(STATUS "Found Eigen3 in a default directory from `${EIGEN3_INCLUDE_DIR}`")
            list(APPEND CMAKE_INCLUDE_PATH "${EIGEN3_INCLUDE_DIR}")
        else ()
            if (NOT EIGEN3_SOURCE AND NOT EIGEN3_REPO)
                set(EIGEN3_SOURCE https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.bz2)
            endif ()
        endif ()
    endif ()

    if (NOT Eigen3_FOUND)
        # Newer versions of Eigen will only generate EigenConfig.cmake when it was top level project on default,
        # see commit 9700fc84. Without EigenConfig.cmake, the following error will occur
        #     export called with target "SecChem" which requires target "eigen" that is not in any export set.
        # The following line should fix it
        set(EIGEN_BUILD_CMAKE_PACKAGE TRUE CACHE BOOL "Enables the creation of EigenConfig.cmake and related files" FORCE)
        set(BUILD_TESTING FALSE CACHE BOOL "Enable creation of Eigen tests" FORCE)

        message(STATUS "Eigen3 (version 5.0) not found and will be fetched from `${EIGEN3_SOURCE}`")

        include(FetchContent)
        FetchContent_Declare(
                eigen
                URL ${EIGEN3_SOURCE}
                GIT_REPOSITORY ${EIGEN3_REPO}
                GIT_TAG ${EIGEN3_TAG}
                GIT_PROGRESS
                PREFIX ${CMAKE_CURRENT_BINARY_DIR}/eigen)
        FetchContent_MakeAvailable(eigen)
    endif ()
endif ()
target_link_libraries(SecChem INTERFACE Eigen3::Eigen)

if (TARGET range-v3::range-v3)
    message(STATUS "SecChem will use the Eigen3 target defined by ancestor project")
else ()
    if (NOT RANGE_V3_SOURCE AND NOT RANGE_V3_REPO AND NOT RANGE_V3_TAG)
        find_package(range-v3)
    endif ()

    if (NOT range-v3_FOUND)
        if (NOT RANGE_V3_SOURCE AND NOT RANGE_V3_REPO AND NOT RANGE_V3_TAG)
             set(RANGE_V3_SOURCE https://github.com/ericniebler/range-v3/archive/refs/tags/0.12.0.tar.gz)
            # set(RANGE_V3_REPO https://github.com/ericniebler/range-v3.git)
            # set(RANGE_V3_TAG 0.12.0)
        endif ()

        include(FetchContent)
        FetchContent_Declare(
                range-v3
                URL ${RANGE_V3_SOURCE}
                GIT_REPOSITORY ${RANGE_V3_REPO}
                GIT_TAG ${RANGE_V3_TAG}
                UPDATE_DISCONNECTED ON
        )
        FetchContent_MakeAvailable(range-v3)
    endif ()
endif ()
target_link_libraries(SecChem INTERFACE range-v3::range-v3)

enable_testing()
add_subdirectory(tests)

#-----------------------------------------------------------------------------------

set_target_properties(SecChem PROPERTIES EXPORT_NAME SecChem)
install(TARGETS SecChem EXPORT SecChemTargets)
export(TARGETS SecChem NAMESPACE SecChem:: FILE SecChemTargets.cmake)
export(PACKAGE SecChem)
install(EXPORT SecChemTargets NAMESPACE SecChem:: DESTINATION "${CMAKE_INSTALL_DATADIR}/SecChem/cmake")
install(DIRECTORY "SecChem" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

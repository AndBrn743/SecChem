if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif ()

option(USE_SYSTEM_CATCH2 NO)

if (USE_SYSTEM_CATCH2)
    find_package(Catch2 3.7.1 CONFIG)
endif ()

if (NOT USE_SYSTEM_CATCH2 OR NOT Catch2_FOUND)
    if (NOT USE_SYSTEM_CATCH2)
        message(STATUS "Building Catch2 Locally")
    else ()
        message(STATUS "Could NOT Find Catch2 (Building Locally)")
    endif ()
    include(FetchContent)
    if (DEFINED CATCH2_SOURCE)
        message(STATUS "Building Catch2 with user specified source")
        FetchContent_Declare(
                catch2
                URL ${CATCH2_SOURCE})
    else ()
        message(STATUS "Building Catch2 with source from internet")
        FetchContent_Declare(
                catch2
                GIT_REPOSITORY https://github.com/catchorg/Catch2.git
                GIT_PROGRESS TRUE
                GIT_TAG v3.7.1)
    endif ()
    FetchContent_MakeAvailable(catch2)
else ()
    message(STATUS "Found Catch2 version ${Catch2_VERSION} from ${Catch2_DIR}")
endif ()

if (USE_SYSTEM_NLOHMANN_JSON)
    find_package(nlohmann_json 3.12.0)
endif ()

if (NOT USE_SYSTEM_NLOHMANN_JSON OR NOT nlohmann_json_FOUND)
    include(FetchContent)

    if (NOT DEFINED NLOHMANN_JSON_SOURCE OR NLOHMANN_JSON_SOURCE STREQUAL "")
        set(NLOHMANN_JSON_SOURCE "https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz")
    endif ()

    FetchContent_Declare(NlohmannJson URL ${NLOHMANN_JSON_SOURCE})
    FetchContent_MakeAvailable(NlohmannJson)
    #    target_link_libraries(foo PRIVATE nlohmann_json::nlohmann_json)
else ()
    message(STATUS "Found nlohmann/json version ${nlohmann_json_VERSION}")
endif ()


add_library(BseJsonParser BasisSetExchangeJsonParserExample.cpp)
target_link_libraries(BseJsonParser PUBLIC SecChem nlohmann_json::nlohmann_json)

if (NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-unused-variable")

    if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    endif ()
endif ()

include(CTest)
include(Catch)

# Helper function to create test executable and register with CTest
function(add_secchem_test test_name source_file)
    add_executable(${test_name} ${source_file})
    target_link_libraries(${test_name} PRIVATE Catch2::Catch2WithMain SecChem)

    # For tests that need BSE JSON parser
    if("${test_name}" STREQUAL "BasisParserTest" OR "${test_name}" STREQUAL "MolecularBasisSetTest")
        target_link_libraries(${test_name} PRIVATE BseJsonParser)
    endif()

    catch_discover_tests(${test_name})
endfunction()

# Utility tests
add_secchem_test(PolyfillFlatMapTest PolyfillFlatMapTest.cpp)
add_secchem_test(ProcessTest ProcessTest.cpp)

# Quantum chemistry type tests
add_secchem_test(AzimuthalQuantumNumberTest AzimuthalQuantumNumberTest.cpp)
add_secchem_test(ElectronicSubShellTest ElectronicSubShellTest.cpp)
add_secchem_test(AtomicElectronConfigurationTest AtomicElectronConfigurationTest.cpp)

# Core chemistry tests
add_secchem_test(ElementTest ElementTest.cpp)
add_secchem_test(AtomTest AtomTest.cpp)
add_secchem_test(MoleculeTest MoleculeTest.cpp)

# Geometry/input parsing tests
add_secchem_test(InternalCoordinateLineParseTest InternalCoordinateLineParseTest.cpp)
add_secchem_test(GeometryInputTest GeometryInputTest.cpp)

# Basis set tests
add_secchem_test(ContractedRadialOrbitalSetTest ContractedRadialOrbitalSetTest.cpp)
add_secchem_test(SemiLocalEcpTest SemiLocalEcpTest.cpp)
add_secchem_test(AngularMomentumBlockTest AngularMomentumBlockTest.cpp)
add_secchem_test(BasisSetTest BasisSetTest.cpp)
add_secchem_test(BasisParserTest BasisParserTest.cpp)
add_secchem_test(MolecularBasisSetTest MolecularBasisSetTest.cpp)

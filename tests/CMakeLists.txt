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
target_link_libraries(BseJsonParser PRIVATE SecChem nlohmann_json::nlohmann_json)

if (NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-unused-variable")
endif ()

include(CTest)
include(Catch)

add_executable(PolyfillFlatMapTest PolyfillFlatMapTest.cpp)
target_link_libraries(PolyfillFlatMapTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME PolyfillFlatMapTest COMMAND PolyfillFlatMapTest)

add_executable(AzimuthalQuantumNumberTest AzimuthalQuantumNumberTest.cpp)
target_link_libraries(AzimuthalQuantumNumberTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME AzimuthalQuantumNumberTest COMMAND AzimuthalQuantumNumberTest)

add_executable(ElectronicSubShellTest ElectronicSubShellTest.cpp)
target_link_libraries(ElectronicSubShellTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME ElectronicSubShellTest COMMAND ElectronicSubShellTest)

add_executable(AtomicElectronConfigurationTest AtomicElectronConfigurationTest.cpp)
target_link_libraries(AtomicElectronConfigurationTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME AtomicElectronConfigurationTest COMMAND AtomicElectronConfigurationTest)

add_executable(ElementTest ElementTest.cpp)
target_link_libraries(ElementTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME ElementTest COMMAND ElementTest)

add_executable(AtomTest AtomTest.cpp)
target_link_libraries(AtomTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME AtomTest COMMAND AtomTest)

add_executable(MoleculeTest MoleculeTest.cpp)
target_link_libraries(MoleculeTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME MoleculeTest COMMAND MoleculeTest)

add_executable(InternalCoordinateLineParseTest InternalCoordinateLineParseTest.cpp)
target_link_libraries(InternalCoordinateLineParseTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME InternalCoordinateLineParseTest COMMAND InternalCoordinateLineParseTest)

add_executable(ContractedRadialOrbitalSetTest ContractedRadialOrbitalSetTest.cpp)
target_link_libraries(ContractedRadialOrbitalSetTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME ContractedRadialOrbitalSetTest COMMAND ContractedRadialOrbitalSetTest)

add_executable(SemiLocalEcpTest SemiLocalEcpTest.cpp)
target_link_libraries(SemiLocalEcpTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME SemiLocalEcpTest COMMAND SemiLocalEcpTest)

add_executable(AngularMomentumBlockTest AngularMomentumBlockTest.cpp)
target_link_libraries(AngularMomentumBlockTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME AngularMomentumBlockTest COMMAND AngularMomentumBlockTest)

add_executable(BasisSetTest BasisSetTest.cpp)
target_link_libraries(BasisSetTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME BasisSetTest COMMAND BasisSetTest)

add_executable(BasisParserTest BasisParserTest.cpp)
target_link_libraries(BasisParserTest PRIVATE Catch2::Catch2WithMain SecChem BseJsonParser)
add_test(NAME BasisParserTest COMMAND BasisParserTest)

add_executable(MolecularBasisSetTest MolecularBasisSetTest.cpp)
target_link_libraries(MolecularBasisSetTest PRIVATE Catch2::Catch2WithMain SecChem BseJsonParser)
add_test(NAME MolecularBasisSetTest COMMAND MolecularBasisSetTest)

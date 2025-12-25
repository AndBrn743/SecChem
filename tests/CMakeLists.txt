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
    message(STATUS "Found Catch2 VERSION ${Catch2_VERSION} DIR ${Catch2_DIR}")
endif ()

include(CTest)
include(Catch)

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

add_executable(ContractedRadialOrbitalSetTest ContractedRadialOrbitalSetTest.cpp)
target_link_libraries(ContractedRadialOrbitalSetTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME ContractedRadialOrbitalSetTest COMMAND ContractedRadialOrbitalSetTest)

add_executable(SemiLocalEcpTest SemiLocalEcpTest.cpp)
target_link_libraries(SemiLocalEcpTest PRIVATE Catch2::Catch2WithMain SecChem)
add_test(NAME SemiLocalEcpTest COMMAND SemiLocalEcpTest)

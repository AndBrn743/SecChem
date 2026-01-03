//
// Created by Andy on 11/23/2025.
//

#include <sstream>
#include <catch2/catch_test_macros.hpp>
#include <SecChem/AtomicElectronConfiguration.hpp>

using namespace SecChem;

TEST_CASE("AtomicElectronConfiguration filling", "[AtomicElectronConfiguration]") {
    SECTION("Hydrogen (1 electron)") {
        AtomicElectronConfiguration config(1);
        REQUIRE(config[ElectronicSubShell(1,0)] == 1); // 1s
        REQUIRE(config[ElectronicSubShell(2,0)] == 0); // 2s
    }

    SECTION("Helium (2 electrons)") {
        AtomicElectronConfiguration config(2);
        REQUIRE(config[ElectronicSubShell(1,0)] == 2); // 1s
    }

    SECTION("Lithium (3 electrons)") {
        AtomicElectronConfiguration config(3);
        REQUIRE(config[ElectronicSubShell(1,0)] == 2); // 1s
        REQUIRE(config[ElectronicSubShell(2,0)] == 1); // 2s
    }

    SECTION("Carbon (6 electrons)") {
        AtomicElectronConfiguration config(6);
        REQUIRE(config[ElectronicSubShell(1,0)] == 2); // 1s
        REQUIRE(config[ElectronicSubShell(2,0)] == 2); // 2s
        REQUIRE(config[ElectronicSubShell(2,1)] == 2); // 2p
    }

    SECTION("Oxygen (8 electrons)") {
        AtomicElectronConfiguration config(8);
        REQUIRE(config[ElectronicSubShell(2,1)] == 4); // 2p partially filled
    }

    SECTION("Moving electrons") {
        AtomicElectronConfiguration config(6);
        auto src = ElectronicSubShell(2,0); // 2s
        auto dst = ElectronicSubShell(2,1); // 2p

        int beforeSrc = config[src];
        int beforeDst = config[dst];

        config.MoveElectron(src, dst);
        REQUIRE(config[src] == beforeSrc - 1);
        REQUIRE(config[dst] == beforeDst + 1);

        // Move multiple electrons
        config.MoveElectrons(src, dst, 5); // should only move up to available or capacity
        REQUIRE(config[src] == 0);
        REQUIRE(config[dst] == 4);

        // Move negative electrons is noop
        config.MoveElectrons(src, dst, -5); // should only move up to available or capacity
        REQUIRE(config[src] == 0);
        REQUIRE(config[dst] == 4);
    }
}

TEST_CASE("Equality operators", "[AtomicElectronConfiguration]") {
    AtomicElectronConfiguration c1(6);
    AtomicElectronConfiguration c2(6);
    AtomicElectronConfiguration c3(8);

    REQUIRE(c1 == c2);
    REQUIRE(c1 != c3);
}

TEST_CASE("operator<<", "[AtomicElectronConfiguration]")
{
	std::ostringstream oss{};
        AtomicElectronConfiguration config(53);
	oss << config;
	CHECK(oss.str() == "1s(2) 2s(2) 2p(6) 3s(2) 3p(6) 3d(10) 4s(2) 4p(6) 4d(10) 5s(2) 5p(5) ");
}

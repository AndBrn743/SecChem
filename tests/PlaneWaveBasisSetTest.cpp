#include <catch2/catch_test_macros.hpp>

#include <SecChem/BasisSet/PlaneWave.hpp>


TEST_CASE("PlaneWave::BasisSet", "[PlaneWave]")
{
	SecChem::BasisSet::PlaneWave::BasisSet bss{3.12, 42.69};
}

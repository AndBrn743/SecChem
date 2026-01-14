// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#include <SecChem/Utility/Process.hpp>
#include <catch2/catch_test_macros.hpp>

#include <string>

TEST_CASE("Process: ExecuteCommandStdout simple command", "[process]")
{
	const std::string output = SecUtility::ExecuteCommandStdout("echo hello world");
	REQUIRE(output.find("hello world") != std::string::npos);
}

TEST_CASE("Process: ExecuteCommandStdout with exit code", "[process]")
{
	int exitCode = -999;
	const std::string output = SecUtility::ExecuteCommandStdout("echo test", exitCode);
	REQUIRE(output.find("test") != std::string::npos);
	// Exit code should be valid (non -999)
	REQUIRE(exitCode != -999);
}

TEST_CASE("Process: ExecuteCommandStdout - list files", "[process]")
{
#if defined(_WIN32)
	const std::string output = SecUtility::ExecuteCommandStdout("dir");
#else
	const std::string output = SecUtility::ExecuteCommandStdout("ls");
#endif
	// Should have some output
	REQUIRE(!output.empty());
}

TEST_CASE("Process: ExecuteCommand - basic functionality", "[process]")
{
	const auto result = SecUtility::ExecuteCommand("echo stdout test");
	REQUIRE(result.StandardOutput.find("stdout test") != std::string::npos);
	REQUIRE(result.StandardError.empty());
}

#if defined(__linux__)
TEST_CASE("Process: ExecuteCommand - command that fails", "[process]")
{
	const auto result = SecUtility::ExecuteCommand("ls /nonexistent/path/that/does/not/exist");
	// Should have some output in stdout (stderr redirected to stdout)
	REQUIRE(!result.StandardOutput.empty());
}

TEST_CASE("Process: ExecuteCommandStdout - grep pattern", "[process]")
{
	// Use printf to avoid shell interpretation issues
	const std::string output =
	        SecUtility::ExecuteCommandStdout(R"(printf "apple\nbanana\ncherry\n" | grep banana)");
	REQUIRE(output.find("banana") != std::string::npos);
	REQUIRE(output.find("apple") == std::string::npos);
}

#endif

#if defined(_WIN32)
TEST_CASE("Process: Windows - exit code from failing command", "[process][windows]")
{
	const auto result = SecUtility::ExecuteCommand("cmd.exe /c exit 1");
	REQUIRE(result.ExitCode == 1);
}

TEST_CASE("Process: Windows - environment variable access", "[process][windows]")
{
	const std::string output = SecUtility::ExecuteCommandStdout("echo %USERNAME%");
	REQUIRE(!output.empty());
}
#endif

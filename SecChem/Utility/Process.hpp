// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <array>
#include <cstdio>  // for popen/_popen/pclose/_pclose
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#if defined(_WIN32)
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#endif


namespace SecUtility
{
	/// <summary>
	/// Represents result of executing a command line process.
	/// </summary>
	struct ProcessResult
	{
		std::string StandardOutput;  ///< Output written to stdout
		std::string StandardError;   ///< Output written to stderr
		int ExitCode;                ///< Exit code returned by the process
	};

	namespace Detail
	{
#if defined(_WIN32)
		/// <summary>
		/// RAII wrapper for Windows HANDLE types.
		/// Ensures handles are properly closed on scope exit.
		/// </summary>
		struct HandleDeleter
		{
			void operator()(HANDLE handle) const noexcept
			{
				if (handle != nullptr && handle != INVALID_HANDLE_VALUE)
				{
					CloseHandle(handle);
				}
			}
		};

		using Handle = std::unique_ptr<std::remove_pointer_t<HANDLE>, HandleDeleter>;

		/// <summary>
		/// RAII wrapper for FILE* handles.
		/// Ensures the file is properly closed on scope exit.
		/// </summary>
		struct FileDeleter
		{
			void operator()(FILE* file) const noexcept
			{
				if (file != nullptr)
				{
					std::fclose(file);
				}
			}
		};

		using FileHandle = std::unique_ptr<FILE, FileDeleter>;

		/// <summary>
		/// Creates a pipe and returns read and write handles.
		/// </summary>
		inline bool CreatePipeInternal(Handle& readHandle, Handle& writeHandle)
		{
			SECURITY_ATTRIBUTES sa{};
			sa.nLength = sizeof(SECURITY_ATTRIBUTES);
			sa.bInheritHandle = TRUE;
			sa.lpSecurityDescriptor = nullptr;

			HANDLE read = nullptr;
			HANDLE write = nullptr;

			if (!CreatePipe(&read, &write, &sa, 0))
			{
				return false;
			}

			// Ensure the read handle is not inherited
			SetHandleInformation(read, HANDLE_FLAG_INHERIT, 0);

			readHandle.reset(read);
			writeHandle.reset(write);
			return true;
		}

		/// <summary>
		/// Executes a command on Windows using CreateProcess and captures stdout/stderr.
		/// </summary>
		inline ProcessResult ExecuteWindows(const std::string& command)
		{
			ProcessResult result{};

			Handle stdoutRead, stdoutWrite;
			Handle stderrRead, stderrWrite;

			if (!CreatePipeInternal(stdoutRead, stdoutWrite))
			{
				result.ExitCode = -1;
				result.StandardError = "Failed to create stdout pipe";
				return result;
			}

			if (!CreatePipeInternal(stderrRead, stderrWrite))
			{
				result.ExitCode = -1;
				result.StandardError = "Failed to create stderr pipe";
				return result;
			}

			STARTUPINFOA si{};
			PROCESS_INFORMATION pi{};

			si.cb = sizeof(si);
			si.dwFlags = STARTF_USESTDHANDLES | STARTF_USESHOWWINDOW;
			si.hStdOutput = stdoutWrite.get();
			si.hStdError = stderrWrite.get();
			si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
			si.wShowWindow = SW_HIDE;

			// Create command line with cmd.exe
			std::string cmdLine = "cmd.exe /c " + command;

			if (!CreateProcessA(nullptr,
			                    const_cast<LPSTR>(cmdLine.c_str()),
			                    nullptr,
			                    nullptr,
			                    TRUE,
			                    CREATE_NO_WINDOW,
			                    nullptr,
			                    nullptr,
			                    &si,
			                    &pi))
			{
				result.ExitCode = -1;
				result.StandardError = "Failed to create process: " + command;
				return result;
			}

			// Close write ends of pipes so we can read from them
			stdoutWrite.reset();
			stderrWrite.reset();

			Handle processHandle(pi.hProcess);
			Handle threadHandle(pi.hThread);

			// Read stdout
			{
				FileHandle stdoutFile(_fdopen(_open_osfhandle(reinterpret_cast<intptr_t>(stdoutRead.release()), _O_RDONLY), "r"));
				if (stdoutFile)
				{
					std::ostringstream output;
					static constexpr int BufferSize = 128;
					std::array<char, BufferSize> buffer{};
					while (std::fgets(buffer.data(), BufferSize, stdoutFile.get()) != nullptr)
					{
						output << buffer.data();
					}
					result.StandardOutput = output.str();
				}
			}

			// Read stderr
			{
				FileHandle stderrFile(_fdopen(_open_osfhandle(reinterpret_cast<intptr_t>(stderrRead.release()), _O_RDONLY), "r"));
				if (stderrFile)
				{
					std::ostringstream output;
					static constexpr int BufferSize = 128;
					std::array<char, BufferSize> buffer{};
					while (std::fgets(buffer.data(), BufferSize, stderrFile.get()) != nullptr)
					{
						output << buffer.data();
					}
					result.StandardError = output.str();
				}
			}

			// Wait for process to finish
			WaitForSingleObject(processHandle.get(), INFINITE);

			// Get exit code
			DWORD exitCode = 0;
			if (GetExitCodeProcess(processHandle.get(), &exitCode))
			{
				result.ExitCode = static_cast<int>(exitCode);
			}
			else
			{
				result.ExitCode = -1;
			}

			return result;
		}

#else
		/// <summary>
		/// RAII wrapper for FILE* handles returned by popen().
		/// Ensures the pipe is properly closed on scope exit.
		/// </summary>
		struct PipeDeleter
		{
			void operator()(FILE* pipe) const noexcept
			{
				if (pipe != nullptr)
				{
					pclose(pipe);
				}
			}
		};

		using PipeHandle = std::unique_ptr<FILE, PipeDeleter>;

		/// <summary>
		/// Opens a pipe for reading command output on Unix systems.
		/// </summary>
		inline PipeHandle OpenPipe(const char* command, const char* mode)
		{
			return PipeHandle(popen(command, mode));
		}

		/// <summary>
		/// Reads all available data from a file stream into a string.
		/// </summary>
		inline std::string ReadStreamToString(FILE* stream)
		{
			std::ostringstream output;
			static constexpr int BufferSize = 128;
			std::array<char, BufferSize> buffer{};
			while (std::fgets(buffer.data(), BufferSize, stream) != nullptr)
			{
				output << buffer.data();
			}
			return output.str();
		}

		/// <summary>
		/// Executes a command on Unix systems and captures both stdout and stderr.
		/// </summary>
		inline ProcessResult ExecuteUnix(const std::string& command)
		{
			ProcessResult result{};

			const std::string stdoutCmd = command + " 2>&1";
			auto pipe = OpenPipe(stdoutCmd.c_str(), "r");
			if (pipe == nullptr)
			{
				result.ExitCode = -1;
				result.StandardError = "Failed to open pipe for command: " + command;
				return result;
			}

			result.StandardOutput = ReadStreamToString(pipe.get());
			result.ExitCode = pclose(pipe.release());
			return result;
		}
#endif
	}  // namespace Detail

	/// <summary>
	/// Executes a command line command and returns both stdout and stderr.
	/// </summary>
	/// <param name="command">The command to execute as a string</param>
	/// <returns>ProcessResult containing StandardOutput, StandardError, and ExitCode</returns>
	/// <remarks>
	/// This function is portable across Windows, Linux, and macOS.
	/// - On Windows: Uses CreateProcess/GetExitCodeProcess for accurate exit codes
	/// - On Unix: Uses pclose to get the exit status, with stderr redirected to stdout
	/// </remarks>
	inline ProcessResult ExecuteCommand(const std::string& command)
	{
#if defined(_WIN32)
		return Detail::ExecuteWindows(command);
#else
		return Detail::ExecuteUnix(command);
#endif
	}

	/// <summary>
	/// Executes a command line command and returns only stdout.
	/// </summary>
	/// <param name="command">The command to execute as a string</param>
	/// <param name="exitCode">Output parameter that receives the exit code</param>
	/// <returns>String containing the stdout output</returns>
	/// <remarks>
	/// This is a simplified interface when you only need stdout and want to check the exit code.
	/// </remarks>
	inline std::string ExecuteCommandStdout(const std::string& command, int& exitCode)
	{
		const auto result = ExecuteCommand(command);
		exitCode = result.ExitCode;
		return std::move(result).StandardOutput;
	}

	/// <summary>
	/// Executes a command line command and returns only stdout (no exit code).
	/// </summary>
	/// <param name="command">The command to execute as a string</param>
	/// <returns>String containing the stdout output</returns>
	/// <remarks>
	/// This is the simplest interface when you only care about the output and not the exit code.
	/// </remarks>
	inline std::string ExecuteCommandStdout(const std::string& command)
	{
		return ExecuteCommand(command).StandardOutput;
	}
}  // namespace SecUtility

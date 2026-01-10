# Local CI with act

This document explains how to use `act` (GitHub Actions emulator) to run CI workflows locally on your Arch Linux system.

## Prerequisites

Install required system packages:

```bash
sudo pacman -S act cmake ninja gcc clang libc++ lcov ccache
```

## Running Workflows

### Quick Start

Use the provided helper script:

```bash
# Run all local CI jobs
./scripts/act-local.sh

# Run specific job
./scripts/act-local.sh -j local-gcc
./scripts/act-local.sh -j local-clang-libstdcxx
./scripts/act-local.sh -j local-clang-libcxx
./scripts/act-local.sh -j local-coverage
```

### Using act Directly

If you prefer to use act directly:

```bash
# Run all jobs in local workflow
act -W .github/workflows/local.yml

# Run specific job
act -W .github/workflows/local.yml -j local-gcc

# Run with verbose output
act -W .github/workflows/local.yml -j local-gcc -v

# Run without dry-run (actually execute the workflow)
act -W .github/workflows/local.yml -j local-gcc --dry-run=false
```

## Available Jobs

### local-gcc
Runs tests with system GCC compiler.

### local-clang-libstdcxx
Runs tests with system Clang compiler using libstdc++.

### local-clang-libcxx
Runs tests with system Clang compiler using libc++.

### local-coverage
Builds with coverage flags and generates a coverage report using lcov.

## Architecture

The `local.yml` workflow is designed to:
- Use system compilers (no Docker containers required)
- Work offline (all dependencies are either system packages or fetched by CMake)
- Be fast and lightweight for rapid iteration

## Troubleshooting

### "act: command not found"
Install act: `yay -S act` or `sudo pacman -S act`

### "Missing required packages"
Install missing packages listed in the error message with pacman.

### Coverage report not generated
Ensure lcov is installed: `sudo pacman -S lcov`

### Workflow fails to find CMake
CMake should be installed system-wide: `sudo pacman -S cmake`

## Testing Specific Platforms

To test workflows for other platforms (Linux, macOS, Windows), you can run those workflows directly with act:

```bash
# Test Linux workflow
act -W .github/workflows/linux.yml -j gcc

# Test macOS workflow (requires macOS)
act -W .github/workflows/macos.yml

# Test Windows workflow (requires Windows)
act -W .github/workflows/windows.yml
```

Note: macOS and Windows workflows require those operating systems to work properly.

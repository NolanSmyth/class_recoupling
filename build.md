# Building/Runng the CLASS Library Using CMake

To build the CLASS library using cmake, first install cmake.
This can be done by downloading from their [website](https://cmake.org) or,
using a package manager:
```bash
# MacOS using HomeBrew
brew install cmake
# ArchLinux
sudo pacman -S cmake
# Ubuntu
sudo apt install cmake
```
Next, make a build directory called `cmake-build-debug` (for a release build,
you can make a directory called `cmake-build-release`.) Move into this directory
and run the following to configure the project:
```bash
# Move into the build directory
cd cmake-build-debug
# Initialize cmake (for release mode, change Debug to Release)
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_BUILD_TYPE=Debug ..
```
To build the project run:
```bash
# Run this inside your build directory
make -j
```
This will build the `main/main.c` file into an executable called `class_main`.
To run the main class program, run the following from the **ROOT** (this is
important!) directory:
```bash
# Run this in the root directory of the class code (for release mode, change debug to release)
./cmake-build-debug/class_main explanatory.ini
```
This will run CLASS with your parameters from the `.ini` file.


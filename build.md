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
This will build the `main/main.c` file into an executable called 
`class_main` (located in the project root directory.) To run the 
main class program, run the following from the **root**  directory:
```bash
# Run this in the root directory of the class code (for release mode, change debug to release)
./class_main explanatory.ini
```
This will run CLASS with your parameters from the `.ini` file.

# VSCode
Much of the above is automated by using VSCode. However, to get things to work in VSCode, you will need the following extensions installed:

- `C/C++`
- `CMake`
- `CMake Tools`
- `clangd` (optional, but great c/c++ linting if you have `clang-tidy` install on your machine)

When you open VSCode, the extensions will be triggered and VSCode will begin to configure the CLASS project.

## Debugging CLASS in VSCode
VSCode makes it simple to attach a debugger to the `class_main` executable. To attach a debugger, we will need to add a launch file. To do so, select `Add Configuration` under the `Run` tab and add the following to the `launch.json` file:
```json
{
    "version": "0.2.0",
    "configurations": [
        
        {
            "type": "lldb",
            "request": "launch",
            "name": "Debug",
            "program": "${workspaceRoot}/class_main",
            "args": [
                "explanatory.ini"
            ],
            "cwd": "${workspaceFolder}"
        }
    ]
}
```
Then, add a breakpoint to the `class.c` file by clicking on left side next to the line numbers on the line you wish the debugger to stop. Then, select `Start Debugging` under the `Run` tab.
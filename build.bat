@echo off

set BUILD_TYPE=%1

set cpp_files=*.cpp
::  -fopenmp=libomp
set common_flags=-fopenmp=libomp -mavx2 -Wall -Wextra -Wunused-function -Wuninitialized -Wsign-compare -Wno-writable-strings -fno-exceptions -fno-rtti -fno-asynchronous-unwind-tables -D"_CRT_SECURE_NO_WARNINGS"
set debug_build_flags=-O0 -g -fno-inline-functions
set release_build_flags=-O3

:: -fsanitize=address
:: default build type is debug
set build_flags=%common_flags% %debug_build_flags%

IF "%BUILD_TYPE%" == "release" echo RELEASE BUILD
IF "%BUILD_TYPE%" == "debug"   echo DEBUG BUILD

IF "%BUILD_TYPE%" == "release" set "build_flags=%common_flags% %release_build_flags%"
IF "%BUILD_TYPE%" == "debug"   set "build_flags=%common_flags% %debug_build_flags%"

:: echo %build_flags%

clang *.c %build_flags% -lkernel32.lib

echo DONE!
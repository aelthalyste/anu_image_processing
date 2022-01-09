@echo off

set BUILD_TYPE=%1

:: set common_flags=-fopenmp=libomp -mavx2 -Wall -Wextra -Wunused-function -Wuninitialized -Wsign-compare -Wno-writable-strings -fno-exceptions -fno-rtti -fno-asynchronous-unwind-tables -D"_CRT_SECURE_NO_WARNINGS"
:: set debug_build_flags=-O0 -g -fsanitize=address -fno-inline-functions
:: set release_build_flags=-O3

set common_flags=/arch:AVX2 /W4 /wd4244 -D_CRT_SECURE_NO_WARNINGS=1 /INCREMENTAL:NO /nologo /openmp /ISDL SDL2main.lib SDL2.lib shell32.lib
set debug_build_flags=/Od /DEBUG:FULL /Zi /fsanitize=address
set release_build_flags=/O2 /GL /GS-

:: -fsanitize=address
:: default build type is debug
set build_flags=%common_flags% %debug_build_flags%

IF "%BUILD_TYPE%" == "release" echo RELEASE BUILD
IF "%BUILD_TYPE%" == "debug"   echo DEBUG BUILD

IF "%BUILD_TYPE%" == "release" set "build_flags=%common_flags% %release_build_flags%"
IF "%BUILD_TYPE%" == "debug"   set "build_flags=%common_flags% %debug_build_flags%"

:: echo %build_flags%
call vcvars64.bat

cl *.c %build_flags%

echo DONE!
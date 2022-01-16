#pragma once

#include <Windows.h>
#define inline __inline

#define Assert(exp) do {if (!(exp)) __debugbreak();} while(0);
//#define Assert(exp)


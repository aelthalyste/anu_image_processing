#pragma once

#include <Windows.h>

#define Assert(exp) do {if (!(exp)) __debugbreak();} while(0);
//#define Assert(exp)


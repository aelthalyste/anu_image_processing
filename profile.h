#pragma once
#include <Windows.h>

// some debug stuff
static inline int64_t start_prof(){
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return li.QuadPart;
}
static inline int64_t get_perf_freq(){
    static int64_t cache = 0;
    if(cache == 0){
        LARGE_INTEGER i;
        QueryPerformanceFrequency(&i);
        cache = i.QuadPart;
    }
    return cache;
}

static inline double time_elapsed_ms(int64_t start, int64_t end) {
    return 1000.0*(end-start)/(double)get_perf_freq();
}

static inline double time_elapsed_us(int64_t start, int64_t end) {
    return 1000000.0*(end-start)/(double)get_perf_freq();
}


struct Profile_Scope {
    const char *name;
    double us_elapsed;
    int64_t begin;
    int64_t end;
};
typedef struct Profile_Scope Profile_Scope;

#define PROF_DEFINE(scope, sname) Profile_Scope scope; scope.name=sname; scope.us_elapsed=0;
#define PROF_BEGIN(scope)        scope.begin = start_prof();
#define PROF_END(scope)          scope.end   = start_prof();
#define PROF_ELAPSED_MS(scope)   time_elapsed_ms(scope.begin, scope.end)
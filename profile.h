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

// time elapsed in ms
static inline double time_elapsed_ms(int64_t start){
    return 1000.0*(start_prof() - start)/(double)get_perf_freq();
}


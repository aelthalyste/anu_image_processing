#pragma once

#include <stdint.h>
#include <string.h>
#include <windows.h>

//#define Assert(exp) do {if (!(exp)) __debugbreak();} while(0);
#define Assert(exp) 
#define yield_execution() SwitchToThread();

// from left upper to right bottom
struct RectangleU16 {
    uint16_t x, y, w, h;
};
typedef struct RectangleU16 RectangleU16;

struct Job {
    DWORD (*proc)(void *arg);
    void *argument;
};
typedef struct Job Job;



struct Draw_Rect_Async_Params {
    uint8_t *image;
    int64_t w, h;
    int64_t channels;
    RectangleU16 r;
};
typedef struct Draw_Rect_Async_Params Draw_Rect_Async_Params;


struct Job_List {
    Job *jobs;
    volatile int64_t len;
    volatile int64_t begin;
    volatile int64_t end;
    volatile int64_t rem;
    volatile int32_t terminate_all;
    CRITICAL_SECTION add_job_section;
    CRITICAL_SECTION consume_job_section;
};
typedef struct Job_List Job_List;





struct Process_Image_Result {

    struct {
        uint8_t *mask;
        uint64_t w, h;
    } downsampled_mask;

    struct {
        uint8_t *mask;
        uint64_t w, h;
    } pure_mask;

};
typedef struct Process_Image_Result Process_Image_Result;


struct Downsample_Task {
    // input boundaries
    uint64_t xs, ys;
    uint64_t ye, xe;

    uint64_t window_width;
    uint64_t window_height;

    uint64_t output_indice;
};
typedef struct Downsample_Task Downsample_Task;


struct Downsample_Task_List {
    Downsample_Task *tasks;
    uint64_t output_width;
    uint64_t output_height;
    uint64_t count;
};
typedef struct Downsample_Task_List Downsample_Task_List;


struct YCBCR_Means {
    uint64_t y, cb, cr;
};
typedef struct YCBCR_Means YCBCR_Means;


struct Linear_Allocator {
    void *memory;
    uint64_t size;
    uint64_t used;
    uint64_t aligment;
};
typedef struct Linear_Allocator Linear_Allocator;


struct Allocator_Mark {
    uint64_t internal_mark;
};
typedef struct Allocator_Mark Allocator_Mark;



static inline void init_linear_allocator(Linear_Allocator *allocator, void *memory, uint64_t memory_size, uint64_t aligment) {
    allocator->memory   = memory;
    allocator->size     = memory_size;
    allocator->aligment = aligment ? aligment : 16;
    allocator->used     = 0;
}

static inline void * linear_allocate_aligned(Linear_Allocator *allocator, uint64_t s, uint64_t al) {

    uint64_t al_bonus = (uint64_t)((uint8_t *)allocator->memory + allocator->used) % al;
    if (al_bonus)
        al_bonus = al - al_bonus;
    
    void *result = NULL;
    if (s + al_bonus + allocator->used < allocator->size) {
        result = (uint8_t *)allocator->memory + allocator->used + al_bonus;
        allocator->used += al_bonus + s;
    }
    Assert(result);
    return result;
}

static inline void * linear_allocate(Linear_Allocator *allocator, uint64_t s) {
    return linear_allocate_aligned(allocator, s, allocator->aligment);   
}

static inline void * linear_allocate_zero_aligned(Linear_Allocator *allocator, uint64_t s, uint64_t al) {
    void *result = linear_allocate_aligned(allocator, s, al);
    memset(result, 0, s);
    return result;
}

static inline void * linear_allocate_zero(Linear_Allocator *allocator, uint64_t s) {
    void *result = linear_allocate(allocator, s);
    memset(result, 0, s);
    return result;         
}


static inline void linear_allocator_mark(Linear_Allocator *allocator, Allocator_Mark *mark) {
    mark->internal_mark = allocator->used;
}

static inline void linear_allocator_restore(Linear_Allocator *allocator, Allocator_Mark mark) {
    allocator->used = mark.internal_mark;
}

#include <Windows.h>
static void print(const char* msg, ...) {
    va_list args;
    va_start(args, msg);

    char buffer[1024];
    DWORD length = wvsprintfA(buffer, msg, args);

    DWORD written;
    WriteFile(GetStdHandle(STD_OUTPUT_HANDLE), buffer, length, &written, NULL);

    va_end(args);
}


void rgb_to_ycbcr(uint8_t *__restrict rgb, uint8_t *__restrict ycbcr, YCBCR_Means *means, uint64_t w, uint64_t h, uint64_t channels);
void filter_rgb(uint8_t *__restrict input_rgb, uint8_t *__restrict output_mask, uint64_t w, uint64_t h, uint64_t input_channel_count);
void filter_ycbcr_means(uint8_t *input_ycbcr, uint8_t *output_mask, YCBCR_Means means, int64_t cbcrdiff_threshold, uint64_t w, uint64_t h, uint64_t input_channel_count);

Process_Image_Result process_image(Job_List *jl, uint8_t *image, uint64_t w, uint64_t h, uint64_t channel_count, Linear_Allocator *allocator);

void prepare_downsample_task(uint64_t w, uint64_t h, uint64_t window_w, uint64_t window_h, Downsample_Task_List *output, Linear_Allocator *allocator);
void downsample_mask(uint8_t *mask_input, uint8_t *output, uint64_t w, Downsample_Task_List *list);


void draw_rectangle(uint8_t *image, int64_t w, int64_t h, uint64_t channels, RectangleU16 rect);
void draw_line_hor(uint8_t *image, int64_t w, int64_t h, uint64_t channels, int64_t y, int64_t xs, int64_t xe);
void draw_line_ver(uint8_t *image, int64_t w, int64_t h, uint64_t channels, int64_t x, int64_t ys, int64_t ye);


void init_job_list(Job_List *jl, uint64_t queue_cap, Linear_Allocator *allocator);
void queue_job(Job_List *jl, DWORD (*proc)(void *arg), void *argument);
void synchronize_jobs(Job_List *jl);
DWORD job_poll_thread(void *arg);


void draw_rectangle_async(Job_List *jl, uint8_t *image, int64_t w, int64_t h, int64_t channels, RectangleU16 rect, Linear_Allocator *allocator);
DWORD draw_rectangle_proc_conv(void *arg);
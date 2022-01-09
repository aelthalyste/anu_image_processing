#pragma once

#include <stdint.h>
#include <string.h>
#include <windows.h>

#include "memory.h"
#include "job.h"



// from left upper to right bottom
struct RectangleU {
    uint64_t x, y, w, h;
};
typedef struct RectangleU RectangleU;



struct Draw_Rect_Async_Params {
    uint8_t *image;
    int64_t w, h;
    int64_t channels;
    RectangleU *r;
    uint64_t rc;
};
typedef struct Draw_Rect_Async_Params Draw_Rect_Async_Params;



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
    int64_t y, cb, cr;
};
typedef struct YCBCR_Means YCBCR_Means;


struct RGB_YCBCR_Conversion_Job {
    uint8_t *rgb;
    uint8_t *ycbcr;
    int64_t w, h;
    int64_t channels;
    YCBCR_Means *means;
};
typedef struct RGB_YCBCR_Conversion_Job RGB_YCBCR_Conversion_Job;

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
void rgb_to_ycbcr_avx2(uint8_t *__restrict rgb, uint8_t *__restrict ycbcr, YCBCR_Means *means, uint64_t w, uint64_t h, uint64_t channels);
void filter_rgb(uint8_t *__restrict input_rgb, uint8_t *__restrict output_mask, uint64_t w, uint64_t h, uint64_t input_channel_count);
void filter_ycbcr_means(uint8_t *input_ycbcr, uint8_t *output_mask, YCBCR_Means means, int64_t cbcrdiff_threshold, uint64_t w, uint64_t h, uint64_t input_channel_count);

Process_Image_Result process_image(uint8_t *image, uint64_t w, uint64_t h, uint64_t channel_count, Linear_Allocator *allocator);

void prepare_downsample_task(uint64_t w, uint64_t h, uint64_t window_w, uint64_t window_h, Downsample_Task_List *output, Linear_Allocator *allocator);
void downsample_mask(uint8_t *mask_input, uint8_t *output, uint64_t w, Downsample_Task_List *list);


void draw_rectangle(uint8_t *image, int64_t w, int64_t h, uint64_t channels, RectangleU rect);
void draw_line_hor(uint8_t *image, int64_t w, int64_t h, uint64_t channels, int64_t y, int64_t xs, int64_t xe);
void draw_line_ver(uint8_t *image, int64_t w, int64_t h, uint64_t channels, int64_t x, int64_t ys, int64_t ye);

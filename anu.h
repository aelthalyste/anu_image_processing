#pragma once

#include <stdint.h>
#include "bg.hpp"

// from left upper to right bottom
struct RectangleU16 {
    u16 x, y, w, h;
};


template<typename T>
struct Stack {
    T *values = NULL;
    u64 count = 0;
    void push(T val) {
        values[count++] = val;
    }
    T pop() {
        if (count > 0)
            return values[--count];
    }
};


struct Process_Image_Result {

    struct {
        u8 *mask;
        u64 w, h;
    }downsampled_mask;

    struct {
        u8 *mask;
        u64 w, h;
    }pure_mask;

};


struct Downsample_Task {
    // input boundaries
    u64 xs, ys;
    u64 ye, xe;

    u64 window_width;
    u64 window_height;

    u64 output_indice;
};

struct Downsample_Task_List {
    Downsample_Task *tasks;
    u64 output_width;
    u64 output_height;
    u64 count;
};


struct RGB_YCBCR_Means {
    u64 r, g, b;
    u64 y, cb, cr;
};

struct Image_Size {
    u64 w, h;
};


struct Linear_Allocator {

    void *memory;
    unsigned long long cap;
    unsigned long long used;

    void *allocate(u64 size);
    void reset();

};

Linear_Allocator
init_linear_allocator(void *memory, u64 size);

void
rgb_to_ycbcr(u8 *__restrict rgb, u8 *__restrict ycbcr, RGB_YCBCR_Means *means, u64 w, u64 h, u64 channels);

void
filter_rgb(u8 *__restrict input_rgb, u8 *__restrict output_mask, u64 w, u64 h, u64 input_channel_count);

void
filter_ycbcr_means(u8 *input_ycbcr, u8 *output_mask, RGB_YCBCR_Means means, s64 cbcrdiff_threshold, u64 w, u64 h, u64 input_channel_count);

Process_Image_Result
process_image(u8 *image, u64 w, u64 h, u64 channel_count, Linear_Allocator *allocator);

void
prepare_downsample_task(u64 w, u64 h, u64 window_w, u64 window_h, Downsample_Task_List *output, Linear_Allocator *allocator);


void
downsample_mask(u8 *mask_input, u8 *output, u64 w, Downsample_Task_List *list);


void
draw_rectangle(u8 *image, s64 w, s64 h, u64 channels, RectangleU16 rect);

void
draw_line_hor(u8 *image, s64 w, s64 h, u64 channels, s64 y, s64 xs, s64 xe);

void
draw_line_ver(u8 *image, s64 w, s64 h, u64 channels, s64 x, s64 ys, s64 ye);
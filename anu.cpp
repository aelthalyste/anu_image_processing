#include "anu.h"

#include <immintrin.h>


Linear_Allocator
init_linear_allocator(void *memory, u64 size) {
    Linear_Allocator result;
    result.memory = memory;
    result.used = 0;
    result.cap = size;
    return result;
}


void * Linear_Allocator::allocate(u64 size) {

    if (size % 64 != 0)
        size += (64 - (size % 64));

    BG_ASSERT(size % 64 == 0);
    BG_ASSERT(size + used < cap);
    if (size + used > cap)
        return NULL;

    void *result = (u8 *)memory + used;
    used += size;
    return result;

}

void Linear_Allocator::reset() {
    used = 0;
}

void 
rgb_to_ycbcr_avx2(u8 *__restrict rgb, u8 *__restrict ycbcr, RGB_YCBCR_Means *means, u64 w, u64 h, u64 channels) {
    s64 tr = 0, tg = 0, tb = 0;
    s64 ty = 0, tcb = 0, tcr = 0;

#define movdqu(a, b)  (a) = _mm_loadu_si128((__m128i *)(b))

#pragma omp parallel for num_threads(10)
    for (s64 y = 0; y < (s64)h; ++y) {


        u8 *__restrict input  = &rgb[(y)*w * channels];
        u8 *__restrict output = &ycbcr[(y)*w * channels];

        // 8 pixels - 32bytes, __m256i
        // rgba rgba rgba rgba rgba rgba rgba rgba

        // 4 pixels - 16bytes, __m128i
        // rgba rgba rgba rgba


        s16 y_constants[] ={
            66, 129, 25,0,
            66, 129, 25,0,
            66, 129, 25,0,
            66, 129, 25,0
        };

        s16 cb_constants[] ={
            -38, -75, 112, 0,
            -38, -75, 112, 0,
            -38, -75, 112, 0,
            -38, -75, 112, 0
        };

        s16 cr_constants[] ={
            112, -94, -18, 0,
            112, -94, -18, 0,
            112, -94, -18, 0,
            112, -94, -18, 0
        };

        s16 adder_table[] ={
            16, 16, 128, 128, 128, 128, 0, 0,
            16, 16, 128, 128, 128, 128, 0, 0
        };

        u8 shuffle8_layout[]  ={ 0, 4, 8, 12, 2, 6, 10, 14, 0, 4, 8, 12, 2, 6, 10, 14, 0, 4, 8, 12, 2, 6, 10, 14, 0, 4, 8, 12, 2, 6, 10, 14 };
        u8 shuffle16_layout[] ={ 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15 };

        // rgba-rgba
        __m256i rgb_mean_vector = _mm256_setzero_si256();
        // ybra-ybra
        __m256i mean_vector     = _mm256_setzero_si256();
        __m256i yc_vector;
        __m256i cbc_vector;
        __m256i crc_vector;
        __m256i add_vec;
        __m256i shuffle8_vec;
        __m256i shuffle16_vec;

        yc_vector   = _mm256_loadu_si256((__m256i *) & y_constants);
        cbc_vector  = _mm256_loadu_si256((__m256i *) & cb_constants);
        crc_vector  = _mm256_loadu_si256((__m256i *) & cr_constants);
        add_vec     = _mm256_loadu_si256((__m256i *) & adder_table);
        shuffle8_vec = _mm256_loadu_si256((__m256i *) & shuffle8_layout);
        shuffle16_vec = _mm256_loadu_si256((__m256i *) & shuffle16_layout);

        for (u64 x = 0; x < w; x+=4) {


            __m128i lane;

            lane = _mm_loadu_si128((__m128i *) & input[x * channels]);


            // there are 4 pixels
            __m256i pixels = _mm256_cvtepu8_epi16(lane);


            __m256i l0 = _mm256_mullo_epi16(pixels, yc_vector);
            __m256i l1 = _mm256_mullo_epi16(pixels, cbc_vector);
            __m256i l2 = _mm256_mullo_epi16(pixels, crc_vector);
            __m256i lt = _mm256_setzero_si256();

            /*
                yyyyyyyyyyyyyyyy
                bbbbbbbbbbbbbbbb

                yyyybbbbyyyybbbb
            */
            __m256i l01 = _mm256_hadd_epi16(l0, l1);


            /*
                rrrrrrrrrrrrrrrr
                0000000000000000

                rrrr0000rrrr0000
            */
            __m256i l2t = _mm256_hadd_epi16(l2, lt);

            /*
                yyyybbbbyyyybbbb
                rrrr0000rrrr0000

                yybbrr00yybbrr00
            */
            __m256i cpa = _mm256_hadd_epi16(l01, l2t);
            __m256i shifted_ = _mm256_srai_epi16(cpa, 8);

            // do this!!
            // _mm256_packs_epi16

            // layout
            /*
            * first two pixel's luminance(y), then first two pixel's blue chrominence(b), then red chrominance(r) + zero alfa,
            * same thing repeats for 3-4'th pixel
                yybbrr00yybbrr00
                0123456789012345
            */
            __m256i result_ycbcr = _mm256_add_epi16(shifted_, add_vec);


            /*
           FOR j := 0 to 15
                i := j*8
                IF b[i+7] == 1
                   dst[i+7:i] := 0
                ELSE
                   index[3:0] := b[i+3:i]
                   dst[i+7:i] := a[index*8+7:index*8]
                FI
                IF b[128+i+7] == 1
                   dst[128+i+7:128+i] := 0
                ELSE
                   index[3:0] := b[128+i+3:128+i]
                   dst[128+i+7:128+i] := a[128+index*8+7:128+index*8]
                FI
           ENDFOR
           dst[MAX:256] := 0
           000000001

           0123456789012345
           0a0a0b0b0c0c0d0d 0a0a0b0b0c0c0d0d  A
           xxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxx  B
           abcdabcdabcdabcd 0000000000000000  output


           0123456789012345
           aaaabbbbccccdddd aaaabbbbccccdddd  A
           xxxxxxxxxxxxxxxx xxxxxxxxxxxxxxxx  B
           output
           */

           // 01 45 89 23
           // 23 67 01 45

            __m256i r = _mm256_shuffle_epi8(result_ycbcr, shuffle8_vec);
            // __m256i r2 = _mm256_shuffle_epi8(result_ycbcr, shuffle16_vec);
            ((u64 *)output)[0] = ((u64 *)&r)[0];
            ((u64 *)output)[1] = ((u64 *)&r)[1];

            // memcpy(output, &r, 16);

#if 0
            output[(x + 0) * channels + 0] = result_ycbcr.m256i_i16[0];
            output[(x + 0) * channels + 1] = result_ycbcr.m256i_i16[2];
            output[(x + 0) * channels + 2] = result_ycbcr.m256i_i16[4];
            output[(x + 0) * channels + 3] = 0;

            output[(x + 1) * channels + 0] = result_ycbcr.m256i_i16[1];
            output[(x + 1) * channels + 1] = result_ycbcr.m256i_i16[3];
            output[(x + 1) * channels + 2] = result_ycbcr.m256i_i16[5];
            output[(x + 1) * channels + 3] = 0;


            output[(x + 2) * channels + 0] = result_ycbcr.m256i_i16[8];
            output[(x + 2) * channels + 1] = result_ycbcr.m256i_i16[10];
            output[(x + 2) * channels + 2] = result_ycbcr.m256i_i16[12];
            output[(x + 2) * channels + 3] = 0;

            output[(x + 3) * channels + 0] = result_ycbcr.m256i_i16[9];
            output[(x + 3) * channels + 1] = result_ycbcr.m256i_i16[11];
            output[(x + 3) * channels + 2] = result_ycbcr.m256i_i16[13];
            output[(x + 3) * channels + 3] = 0;
#endif

            __m128i low_r = _mm256_extracti128_si256(result_ycbcr, 0); // low lane,  2 pixel
            __m128i hig_r = _mm256_extracti128_si256(result_ycbcr, 1); // high line, 2 pixel
            __m256i t0 = _mm256_cvtepi16_epi32(low_r);
            __m256i t1 = _mm256_cvtepi16_epi32(hig_r);
            __m256i t2 = _mm256_add_epi32(t0, t1);
            mean_vector     = _mm256_add_epi32(mean_vector, t2);

            __m128i pixels_low  = _mm256_extracti128_si256(pixels, 0); // 2 pixel
            __m128i pixels_high = _mm256_extracti128_si256(pixels, 1); // 2 pixel
            rgb_mean_vector = _mm256_add_epi32(rgb_mean_vector, _mm256_add_epi32(_mm256_cvtepi16_epi32(pixels_low), _mm256_cvtepi16_epi32(pixels_high)));

        }

        u64 tlr = 0;
        u64 tlg = 0;
        u64 tlb = 0;
        u64 tly = 0;
        u64 tlcb = 0;
        u64 tlcr = 0;


        // tlr  = rgb_mean_vector.m256i_i32[0] + rgb_mean_vector.m256i_i32[4];
        // tlg  = rgb_mean_vector.m256i_i32[1] + rgb_mean_vector.m256i_i32[5];
        // tlb  = rgb_mean_vector.m256i_i32[2] + rgb_mean_vector.m256i_i32[6];
        // tly  = mean_vector.m256i_i32[0] + mean_vector.m256i_i32[1];
        // tlcb = mean_vector.m256i_i32[2] + mean_vector.m256i_i32[3];
        // tlcr = mean_vector.m256i_i32[4] + mean_vector.m256i_i32[5];

#if 1

#pragma omp atomic
        tr += tlr;

#pragma omp atomic
        tg += tlg;

#pragma omp atomic
        tb += tlb;

#pragma omp atomic
        ty += tly;

#pragma omp atomic
        tcb += tlcb;

#pragma omp atomic
        tcr+=tlcr;
#endif

    }

    u64 pixel_count = w * h;

    tr /= pixel_count;
    tg /= pixel_count;
    tb /= pixel_count;

    ty  /= pixel_count;
    tcb /= pixel_count;
    tcr /= pixel_count;

    means->r = tr;
    means->g = tg;
    means->b = tb;

    means->y = ty;
    means->cb = tcb;
    means->cr = tcr;

}


void
rgb_to_ycbcr(u8 *__restrict rgb, u8 *__restrict ycbcr, RGB_YCBCR_Means *means, u64 w, u64 h, u64 channels) {
    s64 tr = 0, tg = 0, tb = 0;
    s64 ty = 0, tcb = 0, tcr = 0;

    
#pragma omp parallel for
    for (s64 y = 0; y < (s64)h; ++y) {
        u8 *__restrict input  = &rgb[y * w * channels];
        u8 *__restrict output = &ycbcr[y * w * channels];
        // thread locals
        u64 tlr = 0, tlg = 0, tlb = 0;
        u64 tly = 0, tlcb = 0, tlcr = 0;


        for (u64 x = 0; x < w; ++x) {



            s32 r = input[x * channels + 0];
            s32 g = input[x * channels + 1];
            s32 b = input[x * channels + 2];
            s32 _y  = 16 + ((r * 66 + g * 129 + b * 25) >> 8);
            s32 cb  = 128 + ((-(r * 38) - (g * 75) + b * 112) >> 8);
            s32 cr  = 128 + ((r * 112 - (g * 94) - (b * 18)) >> 8);

            BG_ASSERT(cr < 256);
            BG_ASSERT(cb < 256);
            BG_ASSERT(_y < 256);
            BG_ASSERT(_y >= 16);

            output[x * channels + 0] = _y;
            output[x * channels + 1] = cb;
            output[x * channels + 2] = cr;

            tlr += r;
            tlg += g;
            tlb += b;
            tly += _y;
            tlcb += cb;
            tlcr += cr;
        }
        
#if 1

#pragma omp atomic
        tr += tlr;

#pragma omp atomic
        tg += tlg;

#pragma omp atomic
        tb += tlb;

#pragma omp atomic
        ty += tly;

#pragma omp atomic
        tcb += tlcb;

#pragma omp atomic
        tcr+=tlcr;
#endif

    }

    u64 pixel_count = w * h;

    tr /= pixel_count;
    tg /= pixel_count;
    tb /= pixel_count;

    ty  /= pixel_count;
    tcb /= pixel_count;
    tcr /= pixel_count;

    means->r = tr;
    means->g = tg;
    means->b = tb;

    means->y = ty;
    means->cb = tcb;
    means->cr = tcr;

}

void
filter_rgb(u8 *__restrict input_rgb, u8 *__restrict output_mask, u64 w, u64 h, u64 input_channel_count) {

#pragma omp parallel for num_threads(12)
    for (s64 _y = 0; _y < (s64)h; ++_y) {
        u8 *__restrict input   = &input_rgb[_y * w * input_channel_count];
        u8 *__restrict output  = &output_mask[_y * w];

        for (u64 x = 0; x < w; ++x) {
            u8 r = input[x * input_channel_count + 0];
            u8 g = input[x * input_channel_count + 1];
            u8 b = input[x * input_channel_count + 2];
            u8 mask = (r > 210) & (g < 200) & (b < 200);
            output[x] &= (mask * 0xff);
        }
    }

}

void
filter_ycbcr_means(u8 *input_ycbcr, u8 *output_mask, RGB_YCBCR_Means in_means, s64 cbcrdiff_threshold, u64 w, u64 h, u64 input_channel_count) {

#pragma omp parallel for num_threads(12)
    for (s64 _y = 0; _y < (s64)h; ++_y) {

        u8 *__restrict input   = &input_ycbcr[_y * w * input_channel_count];
        u8 *__restrict output  = &output_mask[_y * w];
        RGB_YCBCR_Means means = in_means;

        for (u64 x = 0; x < w; ++x) {
            u8 y  = input[x * input_channel_count + 0];
            u8 cb = input[x * input_channel_count + 1];
            u8 cr = input[x * input_channel_count + 2];

            u8 mask = 0;
            s64 cbcrdiff = (s64)cb - (s64)cr;
            if (cbcrdiff < 0)
                cbcrdiff = -cbcrdiff;

            mask = (y > means.y) & (cb < means.cb) & (cr > means.cr) & (cbcrdiff > cbcrdiff_threshold);
            output[x] = mask * 0xff;
        }

    }

}

void
apply_binary_mask_to_image(u8 *image, u8 *binary_mask, u64 w, u64 h, u64 image_channel_count) {
    
#pragma omp parallel for num_threads(12)
    for (s64 _y = 0; _y < (s64)h; ++_y) {

        u8 *__restrict input   = &binary_mask[_y * w];
        u8 *__restrict output  = &image[_y * w * image_channel_count];

        for (u64 x = 0; x < w; ++x) {
            if (input[x]) {
                output[x * image_channel_count + 0] = 0xff; //r
                output[x * image_channel_count + 1] = 0x00; //g
                output[x * image_channel_count + 2] = 0xff; //b
            }
        }

    }
}




void
prepare_downsample_task(u64 w, u64 h, u64 window_w, u64 window_h, Downsample_Task_List *output, Linear_Allocator *allocator) {
    
    u64 xcount = (w + window_w - 1)/ window_w;
    u64 ycount = (h + window_h -1 )/ window_h;
    u64 tile_count = ycount * xcount;

    output->output_width = xcount;
    output->output_height = ycount;

    output->tasks = (Downsample_Task *)allocator->allocate(tile_count * sizeof(Downsample_Task));
    output->count = tile_count;

    u64 ti = 0;

    for (u64 yi = 0; 
        yi < ycount; 
        ++yi) 
    {
        u64 ystart = yi * window_h;
        u64 yend   = ystart + window_h;
        if (yend > h)
            yend = h;

        for (u64 xi = 0;
            xi < xcount;
            ++xi) 
        {
            u64 xstart = xi * window_w;
            u64 xend   = xstart + window_w;
            if (xend > w)
                xend = w;

            output->tasks[ti].xs = xstart;
            output->tasks[ti].xe = xend;

            output->tasks[ti].ys = ystart;
            output->tasks[ti].ye = yend;

            output->tasks[ti].output_indice = ti;
            
            output->tasks[ti].window_height = window_h;
            output->tasks[ti].window_width  = window_w;
            ++ti;
        }

    }

}


void
downsample_mask(u8 *mask_input, u8 *output, u64 w, Downsample_Task_List *list) {
    
    s64 task_count = list->count;

#pragma omp parallel for num_threads(12)
    for (s64 _ti = 0; _ti < task_count; ++_ti) {

        Downsample_Task *task = &list->tasks[_ti];

        u8 found = 0;

        for (u64 ys = task->ys; ys < task->ye; ++ys) {
            u8 *__restrict input = &mask_input[ys * w];

            for (u64 xs = task->xs; xs < task->xe; ++xs) {
                found = input[xs];
                if (found)
                    goto early_exit;
            }
        }
        
    early_exit:
        output[task->output_indice] = found;
    }
    

}




RectangleU16 *
extract_rectangles(u8 *mask, u64 w, u64 h, u64 *rectangle_count_output, Linear_Allocator *allocator) {
    
    // find rectangles with floodfill algorithm
    u8 *temp_mask = (u8 *)allocator->allocate(w * h);
    memcpy(temp_mask, mask, w * h);
    
    u8 visited = 0x01;
    

    struct Coordinates{
        s16 x, y;
    };
    

    u64 rc = 0;
    u64 max_rc_count = w * h / 8;
    RectangleU16 *rectangles = (RectangleU16 *)allocator->allocate((max_rc_count) * sizeof(*rectangles));

    Stack<Coordinates> stack;
    stack.values = (Coordinates *)allocator->allocate(w * h * sizeof(stack.values[0]));
    stack.count = 0;

    for (u64 y = 0; y < h; ++y) {

        u8 *__restrict input = &temp_mask[y * w];

        for (u64 x = 0; x < w; ++x) {
            
            if (input[x] == visited)
                continue;
            
            if (input[x] == 0)
                continue;

#if 1
            rectangles[rc++] = RectangleU16{ (u16)x, (u16)y, 1, 1};
#else

            // safe to convert
            stack.push(Coordinates{ (s16)x, (s16)y });
            
            Coordinates left_upper   = { x,y };
            Coordinates right_bottom = { x,y };

            for (;stack.count;) {
                
                Coordinates c = stack.pop();
                u8 m = temp_mask[c.y * w + c.x];
                if (m == visited) {
                    continue;
                }
                if (m == 0) {
                    continue;
                }

                temp_mask[c.y * w + c.x] = visited;

                // higher y means lower in the picture
                if (c.x < left_upper.x)
                    left_upper.x = c.x;
                if (c.y < left_upper.y)
                    left_upper.y = c.y;


                // lower y means higher in the picture
                if (c.x > right_bottom.x)
                    right_bottom.x = c.x;
                if (c.y > right_bottom.y)
                    right_bottom.y = c.y;

                /*
                    right, upper, left, down
                */
                // right
                if (c.x+1 < w)
                    stack.push(Coordinates{ (s16)(c.x +  1), (s16)(c.y + 0) });
                // upper
                if (c.x+1 <w && c.y != 0)
                    stack.push(Coordinates{ (s16)(c.x +  1), (s16)(c.y - 1) });
                // left
                if (c.x != 0)
                    stack.push(Coordinates{ (s16)(c.x - 1), (s16)(c.y + 0) });
                // down
                if (c.y < h)
                    stack.push(Coordinates{ (s16)(c.x +  0), (s16)(c.y + 1) });

                rectangles[rc++] = RectangleU16{ (u16)c.x, (u16)c.y, 1, 1 };

            }

            BG_ASSERT(rc + 1 < max_rc_count);
            // RectangleU16 r;
          // r.x = left_upper.x;
          // r.y = left_upper.y;
          // r.w = (right_bottom.x - left_upper.x) + 1;
          // r.h = (right_bottom.y - left_upper.y) + 1;

          // BG_ASSERT(r.w);
          // BG_ASSERT(r.h);
          // 
          // rectangles[rc] = r;
          // rc+=1;

#endif

      
        }

    }

    *rectangle_count_output = rc;
    return rectangles;

}




/*
avx
4.59
4.63
4.58

normal
15.34
15.23
15.23


avx2 + openmp
2.78
2.56
2.63

normal + openmp
4.62
4.64
4.49
*/

Process_Image_Result
process_image(u8 *image, u64 w, u64 h, u64 channel_count, Linear_Allocator *allocator) {

    Process_Image_Result result;
    zero_memory(&result, sizeof(result));
    result.pure_mask.w = w;
    result.pure_mask.h = h;
    
    u8 *ycbcr      = (u8 *)allocator->allocate(w * h * channel_count);
    u8 *ycbcr_avx2 = (u8 *)allocator->allocate(w * h * channel_count);
    result.pure_mask.mask = (u8 *)allocator->allocate(w * h);
    
    zero_memory(ycbcr     , w * h * channel_count);
    zero_memory(ycbcr_avx2, w * h * channel_count);

    Downsample_Task_List downsample_tasks;
    u64 subsample_w_window = 16;
    u64 subsample_h_window = 16;
    
    BG_ASSERT(ycbcr);

    RGB_YCBCR_Means means;
    RGB_YCBCR_Means means_avx2;

    rgb_to_ycbcr_avx2(image, ycbcr_avx2, &means_avx2, w, h, channel_count);
    //rgb_to_ycbcr(image, ycbcr, &means, w, h, channel_count);

    filter_ycbcr_means(ycbcr, result.pure_mask.mask, means, 38, w, h, channel_count);
    filter_rgb(image, result.pure_mask.mask, w, h, channel_count);
    apply_binary_mask_to_image(image, result.pure_mask.mask, w, h, channel_count);
    
#if 0
    prepare_downsample_task(w, h, subsample_w_window, subsample_h_window, &downsample_tasks, allocator);
    result.downsampled_mask.w = downsample_tasks.output_width;
    result.downsampled_mask.h = downsample_tasks.output_height;

    result.downsampled_mask.mask = (u8 *)allocator->allocate(downsample_tasks.output_height * downsample_tasks.output_width);
    downsample_mask(result.pure_mask.mask, result.downsampled_mask.mask, w, &downsample_tasks);
#endif
    // u64 rect_count = 0;
    // RectangleU16 *rectangles = extract_rectangles(result.downsampled_mask.mask, result.downsampled_mask.w, result.downsampled_mask.h, &rect_count, allocator);

    
    // for (u64 i =0; i < rect_count; ++i) {
    //     
    //     RectangleU16 r = rectangles[i];
    //     
    //     r.x = (u64)r.x * subsample_w_window;
    //     r.y = (u64)r.y * subsample_h_window;
    //     
    //     r.w = (u64)r.w * subsample_w_window;
    //     r.h = (u64)r.h * subsample_h_window;
    // 
    //     draw_rectangle(image, w, h, channel_count, r);
    // }

    return result;
}

void
draw_rectangle(u8 *image, s64 w, s64 h, u64 channels, RectangleU16 rect) {

    draw_line_hor(image, w, h, channels, rect.y         , rect.x, rect.x + rect.w);
    draw_line_hor(image, w, h, channels, rect.y + rect.h, rect.x, rect.x + rect.w);

    draw_line_ver(image, w, h, channels, rect.x         , rect.y, rect.y + rect.h);
    draw_line_ver(image, w, h, channels, rect.x + rect.w, rect.y, rect.y + rect.h);

}

void
draw_line_hor(u8 *image, s64 w, s64 h, u64 channels, s64 y, s64 xs, s64 xe) {
    
    u64 line_width = 5;
    
    s64 ly = y - line_width / 2;
    s64 uy = y + line_width / 2;
    if (ly < 0)
        ly = 0;
    if (uy > h)
        uy = h;

    if (xe > w)
        xe = w;

    for (s64 _y = ly; _y < uy; ++_y) {
        u8 *input = &image[_y * w * channels];

        for (s64 x = xs; x < xe; ++x) {
            input[x * channels + 0] = 0x00;
            input[x * channels + 1] = 0xff;
            input[x * channels + 2] = 0xbb;
        }

    }


}

void 
draw_line_ver(u8 *image, s64 w, s64 h, u64 channels, s64 x, s64 ys, s64 ye) {
    
    u64 line_width = 5;
    
    s64 lx = x - line_width / 2;
    s64 ux = x + line_width / 2;
    if (lx < 0)
        lx = 0;
    if (ux > w)
        ux = w;

    if (ye > h)
        ye = h;

    for (s64 _y = ys; _y < ye; ++_y) {
        u8 *input = &image[_y * w * channels];
        for (s64 _x = lx; _x < ux; ++_x) {
            input[_x * channels + 0] = 0x00;
            input[_x * channels + 1] = 0xff;
            input[_x * channels + 2] = 0xbb;
        }
    }

}
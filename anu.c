#include "anu.h"

#include <string.h>
#include <immintrin.h>
#include <intrin.h>


struct RGB_YCBCR_Conversion_Job {
    uint8_t *rgb;
    uint8_t *ycbcr;
    int64_t w, h;
    int64_t channels;
};


void rgb_to_ycbcr_avx2(uint8_t *__restrict rgb, uint8_t *__restrict ycbcr, YCBCR_Means *means, uint64_t w, uint64_t h, uint64_t channels) {
    

    int64_t ty = 0, tcb = 0, tcr = 0;

#pragma omp parallel for
    for (int64_t y = 0; y < (int64_t)h; ++y) {

        uint8_t *__restrict input  = &rgb[(y)*w * channels];
        uint8_t *__restrict output = &ycbcr[(y)*w * channels];

        // 8 pixels - 32bytes, __m256i
        // rgba rgba rgba rgba rgba rgba rgba rgba

        // 4 pixels - 16bytes, __m128i
        // rgba rgba rgba rgba


        int16_t y_constants[] ={
            66, 129, 25,0,
            66, 129, 25,0,
            66, 129, 25,0,
            66, 129, 25,0
        };

        int16_t cb_constants[] ={
            -38, -75, 112, 0,
            -38, -75, 112, 0,
            -38, -75, 112, 0,
            -38, -75, 112, 0
        };

        int16_t cr_constants[] ={
            112, -94, -18, 0,
            112, -94, -18, 0,
            112, -94, -18, 0,
            112, -94, -18, 0
        };

        int16_t adder_table[] ={
            16, 16, 128, 128, 128, 128, 0, 0,
            16, 16, 128, 128, 128, 128, 0, 0
        };

        uint8_t shuffle8_layout[]  ={ 0, 4, 8, 12, 2, 6, 10, 14, 0, 4, 8, 12, 2, 6, 10, 14, 0, 4, 8, 12, 2, 6, 10, 14, 0, 4, 8, 12, 2, 6, 10, 14 };
        
        // ybra-ybra
        __m256i ycbcr_mean_vector = _mm256_setzero_si256();
        
        __m256i yc_vector;
        __m256i cbc_vector;
        __m256i crc_vector;
        __m256i add_vec;
        __m256i shuffle8_vec;

        yc_vector     = _mm256_loadu_si256((__m256i *) & y_constants);
        cbc_vector    = _mm256_loadu_si256((__m256i *) & cb_constants);
        crc_vector    = _mm256_loadu_si256((__m256i *) & cr_constants);
        add_vec       = _mm256_loadu_si256((__m256i *) & adder_table);
        shuffle8_vec  = _mm256_loadu_si256((__m256i *) & shuffle8_layout);

        for (uint64_t x = 0; x < w; x+=4) {

            __m128i lane;

            lane = _mm_loadu_si128((__m128i *) &input[x * channels]);


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
            // @TODO : must replace this dumb hadd with blend or permute, since this loop is port5 related and next hadd does nothing special about
            // hadd, its just fancy blend.
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
            *((uint64_t *)(&output[x * channels]) + 0) = ((uint64_t *)&r)[0];
            *((uint64_t *)(&output[x * channels]) + 1) = ((uint64_t *)&r)[1];
            
            //t0 = _mm256_add_epi32(_mm256_cvtepuint8_t_epi32(low_r), _mm256_cvtepuint8_t_epi32(hig_r));
            //ycbcr_mean_vector = _mm_add_epi32(ycbcr_mean_vector, _mm_add_epi32(_mm256_extracti128_si256(t0, 0), _mm256_extracti128_si256(t0, 1)));


            ycbcr_mean_vector = _mm256_add_epi32(ycbcr_mean_vector, _mm256_add_epi32(_mm256_cvtepu8_epi32(_mm256_extracti128_si256(r, 0)), _mm256_cvtepu8_epi32(_mm256_extracti128_si256(r, 1))));

        }

        uint64_t tly  = ((int32_t *)&ycbcr_mean_vector)[0] + ((int32_t *)&ycbcr_mean_vector)[4];
        uint64_t tlcb = ((int32_t *)&ycbcr_mean_vector)[1] + ((int32_t *)&ycbcr_mean_vector)[5];
        uint64_t tlcr = ((int32_t *)&ycbcr_mean_vector)[2] + ((int32_t *)&ycbcr_mean_vector)[6];


#if 0
#pragma omp atomic
        ty += tly;

#pragma omp atomic
        tcb += tlcb;

#pragma omp atomic
        tcr+=tlcr;
#else
        InterlockedAdd64(&ty , tly);
        InterlockedAdd64(&tcb, tlcb);
        InterlockedAdd64(&tcr, tlcr);
#endif

    }

    uint64_t pixel_count = w * h;

    ty  /= pixel_count;
    tcb /= pixel_count;
    tcr /= pixel_count;

    means->y = ty;
    means->cb = tcb;
    means->cr = tcr;

}


void
rgb_to_ycbcr(uint8_t *__restrict rgb, uint8_t *__restrict ycbcr, YCBCR_Means *means, uint64_t w, uint64_t h, uint64_t channels) {
    int64_t ty = 0, tcb = 0, tcr = 0;

#pragma omp parallel for
    for (int64_t y = 0; y < (int64_t)h; ++y) {
        uint8_t *__restrict input  = &rgb[y * w * channels];
        uint8_t *__restrict output = &ycbcr[y * w * channels];
        // thread locals
        uint64_t tly = 0, tlcb = 0, tlcr = 0;


        for (uint64_t x = 0; x < w; ++x) {

            int32_t r = input[x * channels + 0];
            int32_t g = input[x * channels + 1];
            int32_t b = input[x * channels + 2];
            int32_t _y  = 16 + ((r * 66 + g * 129 + b * 25) >> 8);
            int32_t cb  = 128 + ((-(r * 38) - (g * 75) + b * 112) >> 8);
            int32_t cr  = 128 + ((r * 112 - (g * 94) - (b * 18)) >> 8);

            Assert(cr < 256);
            Assert(cb < 256);
            Assert(_y < 256);
            Assert(_y >= 16);

            
            output[x * channels + 0] = _y;            
            output[x * channels + 1] = cb;
            output[x * channels + 2] = cr;
            

            tly += _y;
            tlcb += cb;
            tlcr += cr;



        }
        

#if 0
#pragma omp atomic
        ty += tly;

#pragma omp atomic
        tcb += tlcb;

#pragma omp atomic
        tcr+=tlcr;
#else
        InterlockedAdd64(&ty , tly);
        InterlockedAdd64(&tcb, tlcb);
        InterlockedAdd64(&tcr, tlcr);
#endif

    }

    uint64_t pixel_count = w * h;

    ty  /= pixel_count;
    tcb /= pixel_count;
    tcr /= pixel_count;

    means->y = ty;
    means->cb = tcb;
    means->cr = tcr;

}

void
filter_rgb_avx2(uint8_t *__restrict input_rgb, uint8_t *__restrict output_mask, uint64_t w, uint64_t h, uint64_t input_channel_count) {
    
#pragma omp parallel for
    for (int64_t _y = 0; _y < (int64_t)h; ++_y) {

        uint8_t *__restrict input   = &input_rgb[_y * w * input_channel_count];
        uint8_t *__restrict output  = &output_mask[_y * w];

        // we can operate 32 bytes which is 8 pixels at a time.
        for (uint64_t x = 0; x < w; ++x) {
            uint8_t r = input[x * input_channel_count + 0];
            uint8_t g = input[x * input_channel_count + 1];
            uint8_t b = input[x * input_channel_count + 2];
            uint8_t mask = (r > 220) & (g < 190) & (b < 190);
            output[x] &= (mask * 0xff);
        }
    }

}


void
filter_rgb(uint8_t *__restrict input_rgb, uint8_t *__restrict output_mask, uint64_t w, uint64_t h, uint64_t input_channel_count) {

#pragma omp parallel for
    for (int64_t _y = 0; _y < (int64_t)h; ++_y) {
        uint8_t *__restrict input   = &input_rgb[_y * w * input_channel_count];
        uint8_t *__restrict output  = &output_mask[_y * w];

        for (uint64_t x = 0; x < w; ++x) {
            uint8_t r = input[x * input_channel_count + 0];
            uint8_t g = input[x * input_channel_count + 1];
            uint8_t b = input[x * input_channel_count + 2];
            uint8_t mask = (r > 230) & (g < 200) & (b < 200);
            //uint8_t mask = (r > 128) & (g < 128) & (b < 60);
            output[x] &= (mask * 0xff);
        }
    }

}

void filter_ycbcr_means(uint8_t *input_ycbcr, uint8_t *output_mask, YCBCR_Means in_means, int64_t cbcrdiff_threshold, uint64_t w, uint64_t h, uint64_t input_channel_count) {
    
    

#pragma omp parallel for
    for (int64_t _y = 0; _y < (int64_t)h; ++_y) {

        uint8_t *__restrict input   = &input_ycbcr[_y * w * input_channel_count];
        uint8_t *__restrict output  = &output_mask[_y * w];
        YCBCR_Means means = in_means;

        for (uint64_t x = 0; x < w; ++x) {
            uint8_t y  = input[x * input_channel_count + 0];
            uint8_t cb = input[x * input_channel_count + 1];
            uint8_t cr = input[x * input_channel_count + 2];

            uint8_t mask = 0;
            int64_t cbcrdiff = (int64_t)cb - (int64_t)cr;
            if (cbcrdiff < 0)
                cbcrdiff = -cbcrdiff;

            mask = (y > means.y) & (cb < means.cb) & (cr > means.cr) & (cbcrdiff > cbcrdiff_threshold) & (cb < 120) & (cr > 150);
            mask &= (y > cb);
            output[x] = mask * 0xff;
        }

    }

}

void
apply_binary_mask_to_image(uint8_t *image, uint8_t *binary_mask, uint64_t w, uint64_t h, uint64_t image_channel_count) {
    
    

#pragma omp parallel for
    for (int64_t _y = 0; _y < (int64_t)h; ++_y) {

        uint8_t *__restrict input   = &binary_mask[_y * w];
        uint8_t *__restrict output  = &image[_y * w * image_channel_count];

        for (uint64_t x = 0; x < w; ++x) {
            if (input[x]) {
                output[x * image_channel_count + 0] = 0x00; //r
                output[x * image_channel_count + 1] = 0xff; //g
                output[x * image_channel_count + 2] = 0x00; //b
            }
        }

    }
}




void
prepare_downsample_task(uint64_t w, uint64_t h, uint64_t window_w, uint64_t window_h, Downsample_Task_List *output, Linear_Allocator *allocator) {
    
    uint64_t xcount = (w + window_w - 1)/ window_w;
    uint64_t ycount = (h + window_h -1 )/ window_h;
    uint64_t tile_count = ycount * xcount;

    output->output_width = xcount;
    output->output_height = ycount;

    output->tasks = (Downsample_Task *)linear_allocate(allocator, tile_count * sizeof(Downsample_Task));
    output->count = tile_count;

    uint64_t ti = 0;

    for (uint64_t yi = 0; 
        yi < ycount; 
        ++yi) 
    {
        uint64_t ystart = yi * window_h;
        uint64_t yend   = ystart + window_h;
        if (yend > h)
            yend = h;

        for (uint64_t xi = 0;
            xi < xcount;
            ++xi) 
        {
            uint64_t xstart = xi * window_w;
            uint64_t xend   = xstart + window_w;
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
downsample_mask(uint8_t *mask_input, uint8_t *output, uint64_t w, Downsample_Task_List *list) {
    
    int64_t task_count = list->count;

#pragma omp parallel for
    for (int64_t _ti = 0; _ti < task_count; ++_ti) {

        Downsample_Task *task = &list->tasks[_ti];
        Assert(task->xe > task->xs);
        Assert(task->ye > task->ys);

        uint64_t fc=0;

        for (uint64_t ys = task->ys; ys < task->ye; ++ys) {
            uint8_t *__restrict input = &mask_input[ys * w];            
            for (uint64_t xs = task->xs; xs < task->xe; ++xs)
                fc += (input[xs] > 0);                
            
        }
        
        uint64_t pc   = (task->xe - task->xs) * (task->ye - task->ys);
        output[task->output_indice] = (100 * fc) / pc;
    }
    

}




RectangleU16 *
extract_rectangles(uint8_t *mask, uint64_t w, uint64_t h, uint64_t *rectangle_count_output, Linear_Allocator *allocator) {
    
    // find rectangles with floodfill algorithm
    uint8_t *temp_mask = (uint8_t *)linear_allocate(allocator, w * h);
    memcpy(temp_mask, mask, w * h);
    

    struct Coordinates {
        int16_t x, y;
    };
    

    uint64_t rc = 0;
    uint64_t max_rc_count = w * h / 2;
    RectangleU16 *rectangles = (RectangleU16 *)linear_allocate(allocator, (max_rc_count) * sizeof(*rectangles));


    for (uint64_t y = 0; y < h; ++y) {

        uint8_t *__restrict input = &temp_mask[y * w];

        for (uint64_t x = 0; x < w; ++x) {
        
            RectangleU16 r;
            r.x = x;
            r.y = y;
            r.w = 1;
            r.h = 1;
#if 1
            if (input[x] > 5)
                rectangles[rc++] = r;
#else

            if (input[x] == visited)
                continue;
            
            if (input[x] == 0)
                continue;

            // safe to convert
            stack.push(Coordinates{ (int16_t)x, (int16_t)y });
            
            Coordinates left_upper   = { x,y };
            Coordinates right_bottom = { x,y };

            for (;stack.count;) {
                
                Coordinates c = stack.pop();
                uint8_t m = temp_mask[c.y * w + c.x];
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
                    stack.push(Coordinates{ (int16_t)(c.x +  1), (int16_t)(c.y + 0) });
                // upper
                if (c.x+1 <w && c.y != 0)
                    stack.push(Coordinates{ (int16_t)(c.x +  1), (int16_t)(c.y - 1) });
                // left
                if (c.x != 0)
                    stack.push(Coordinates{ (int16_t)(c.x - 1), (int16_t)(c.y + 0) });
                // down
                if (c.y < h)
                    stack.push(Coordinates{ (int16_t)(c.x +  0), (int16_t)(c.y + 1) });

                rectangles[rc++] = RectangleU16{ (uint16_t)c.x, (uint16_t)c.y, 1, 1 };

            }

            Assert(rc + 1 < max_rc_count);
            // RectangleU16 r;
          // r.x = left_upper.x;
          // r.y = left_upper.y;
          // r.w = (right_bottom.x - left_upper.x) + 1;
          // r.h = (right_bottom.y - left_upper.y) + 1;

          // Assert(r.w);
          // Assert(r.h);
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
process_image(Job_List *jl, uint8_t *image, uint64_t w, uint64_t h, uint64_t channel_count, Linear_Allocator *allocator) {
    
    (void)(jl);

    Process_Image_Result result;
    memset(&result, 0, sizeof(result));
    result.pure_mask.w = w;
    result.pure_mask.h = h;
    
    uint8_t *ycbcr      = (uint8_t *)linear_allocate(allocator, w * h * channel_count);

    result.pure_mask.mask = (uint8_t *)linear_allocate(allocator, w * h);    

    Downsample_Task_List downsample_tasks;
    uint64_t subsample_w_window = 32;
    uint64_t subsample_h_window = 32;
    
    Assert(ycbcr);

    // bool avx2_supported = true;
    YCBCR_Means means;
    int64_t avxclocks    = 0;
    int64_t normalclocks = 0;
    int64_t ycbcrfilterclocks = 0;
    int64_t downsampleclocks = 0;
    int64_t binarymaskclocks = 0;

    uint32_t tvforrdt = 0;
    
    avxclocks -=__rdtscp(&tvforrdt);
    rgb_to_ycbcr_avx2(image, ycbcr, &means, w, h, channel_count);
    avxclocks += __rdtscp(&tvforrdt);

    normalclocks -= __rdtscp(&tvforrdt);
    // rgb_to_ycbcr(image, ycbcr, &means, w, h, channel_count);
    normalclocks += __rdtscp(&tvforrdt);


#if 1

    ycbcrfilterclocks -= __rdtscp(&tvforrdt);
    filter_ycbcr_means(ycbcr, result.pure_mask.mask, means, 40, w, h, channel_count);
    ycbcrfilterclocks += __rdtscp(&tvforrdt);
    
    // filter_rgb(image, result.pure_mask.mask, w, h, channel_count);
    binarymaskclocks -= __rdtscp(&tvforrdt);
    apply_binary_mask_to_image(image, result.pure_mask.mask, w, h, channel_count);
    binarymaskclocks += __rdtscp(&tvforrdt);



    prepare_downsample_task(w, h, subsample_w_window, subsample_h_window, &downsample_tasks, allocator);
    result.downsampled_mask.w = downsample_tasks.output_width;
    result.downsampled_mask.h = downsample_tasks.output_height;

    result.downsampled_mask.mask = (uint8_t *)linear_allocate(allocator, downsample_tasks.output_height * downsample_tasks.output_width);

    downsampleclocks -= __rdtscp(&tvforrdt);
    downsample_mask(result.pure_mask.mask, result.downsampled_mask.mask, w, &downsample_tasks);
    downsampleclocks += __rdtscp(&tvforrdt);

    double million = 1000000.0;
    print("avx clocks %10.4f, normal %10.4f, gain : %.4fx, ycbcrfilter %10.4f, binarymask %10.4f, downsample mask %10.4f\n", avxclocks/million, normalclocks/million, (double)normalclocks/(double)avxclocks, ycbcrfilterclocks/million, binarymaskclocks/million, downsampleclocks/million);

    //print("avx clocks %10llu, normal %10llu, gain : %.4fx, ycbcrfilter %10llu, binarymask %10llu, downsample mask %10llu\n", avxclocks, normalclocks, (double)normalclocks/(double)avxclocks, ycbcrfilterclocks, binarymaskclocks, downsampleclocks);

    uint64_t rect_count = 0;
    RectangleU16 *rectangles = extract_rectangles(result.downsampled_mask.mask, result.downsampled_mask.w, result.downsampled_mask.h, &rect_count, allocator);

    uint64_t i=0;
    uint64_t rb = 200;

    for (; i+rb-1<rect_count; i+=rb) {
        for(uint64_t j=0;j<rb;++j) {
            rectangles[i + j].x = (uint64_t)rectangles[i + j].x * subsample_w_window;
            rectangles[i + j].y = (uint64_t)rectangles[i + j].y * subsample_h_window;
            
            rectangles[i + j].w = (uint64_t)rectangles[i + j].w * subsample_w_window;
            rectangles[i + j].h = (uint64_t)rectangles[i + j].h * subsample_h_window;            
        }
        draw_rectangle_async(jl, image, w, h, channel_count, &rectangles[i], rb, allocator);
    }


    for (; i<rect_count; ++i) {
        RectangleU16 r = rectangles[i];
        
        r.x = (uint64_t)r.x * subsample_w_window;
        r.y = (uint64_t)r.y * subsample_h_window;
        
        r.w = (uint64_t)r.w * subsample_w_window;
        r.h = (uint64_t)r.h * subsample_h_window;
    
        draw_rectangle(image, w, h, channel_count, r);
    }

    synchronize_jobs(jl);

#if 0
    for (uint64_t i =0; i < rect_count; ++i) {  
        RectangleU16 r = rectangles[i];
        
        r.x = (uint64_t)r.x * subsample_w_window;
        r.y = (uint64_t)r.y * subsample_h_window;
        
        r.w = (uint64_t)r.w * subsample_w_window;
        r.h = (uint64_t)r.h * subsample_h_window;
    
        draw_rectangle(image, w, h, channel_count, r);
    }
#endif

#endif
    

    return result;
}

void draw_rectangle_async(Job_List *jl, uint8_t *image, int64_t w, int64_t h, int64_t channels, RectangleU16 *rects, uint64_t rect_count, Linear_Allocator *allocator) {
    Draw_Rect_Async_Params *p = linear_allocate(allocator, sizeof(*p));
    p->image = image;
    p->w = w;
    p->h = h;
    p->channels = channels;
    p->r  = rects;
    p->rc = rect_count;
    Assert(p->rc < 256);
    queue_job(jl, draw_rectangle_proc_conv, p);
}

DWORD draw_rectangle_proc_conv(void *arg) {
    Draw_Rect_Async_Params *p = arg;
    for(uint64_t i=0; i<p->rc; ++i)
        draw_rectangle(p->image, p->w, p->h, p->channels, p->r[i]);
    return 0;
}

void draw_rectangle(uint8_t *image, int64_t w, int64_t h, uint64_t channels, RectangleU16 rect) {

    draw_line_hor(image, w, h, channels, rect.y         , rect.x, rect.x + rect.w);
    draw_line_hor(image, w, h, channels, rect.y + rect.h, rect.x, rect.x + rect.w);

    draw_line_ver(image, w, h, channels, rect.x         , rect.y, rect.y + rect.h);
    draw_line_ver(image, w, h, channels, rect.x + rect.w, rect.y, rect.y + rect.h);

}

void
draw_line_hor(uint8_t *image, int64_t w, int64_t h, uint64_t channels, int64_t y, int64_t xs, int64_t xe) {
    
    uint64_t line_width = 5;
    
    int64_t ly = y - line_width / 2;
    int64_t uy = y + line_width / 2;
    if (ly < 0)
        ly = 0;
    if (uy > h)
        uy = h;

    if (xe > w)
        xe = w;

    for (int64_t _y = ly; _y < uy; ++_y) {
        uint8_t *input = &image[_y * w * channels];

        for (int64_t x = xs; x < xe; ++x) {
            input[x * channels + 0] = 0x00;
            input[x * channels + 1] = 0xff;
            input[x * channels + 2] = 0xaa;
        }

    }


}

void 
draw_line_ver(uint8_t *image, int64_t w, int64_t h, uint64_t channels, int64_t x, int64_t ys, int64_t ye) {

    uint64_t line_width = 5;
    
    int64_t lx = x - line_width / 2;
    int64_t ux = x + line_width / 2;
    if (lx < 0)
        lx = 0;
    if (ux > w)
        ux = w;

    if (ye > h)
        ye = h;

    for (int64_t _y = ys; _y < ye; ++_y) {
        uint8_t *input = &image[_y * w * channels];
        for (int64_t _x = lx; _x < ux; ++_x) {
            input[_x * channels + 0] = 0x00;
            input[_x * channels + 1] = 0xff;
            input[_x * channels + 2] = 0xbb;
        }
    }

}
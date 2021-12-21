#include <stdio.h>
#include "anu.h"

#define BG_IMPLEMENTATION
#include "bg.hpp"

#include "anu_stb_include.h"

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
    return 1000.0*((double)start_prof() - (double)start)/(double)get_perf_freq();
}



int main() {
        
    //stbi_set_flip_vertically_on_load(true);
    //stbi_flip_vertically_on_write(true);
    int n =0;

    const char *input_path   = "C:\\W\\anu_tests\\data\\ttt";
    const char *output_path  = "C:\\W\\anu_tests\\output\\";
	bg_unused(output_path);
	bg_unused(input_path);

    auto files = get_file_paths_in_directory(input_path);
    
    u64 memory_size = 1024 * 1024 * 128;
    u8 *memory = (u8 *)_aligned_malloc(1024 * 1024 * 128, 64);
    Linear_Allocator allocator = init_linear_allocator(memory, memory_size);

    Array<u8 *> images;
    Array<Image_Size> image_sizes;

    arrreserve(&images, files.len);
    arrreserve(&image_sizes, files.len);

    int channels = 4;
    for_array (fi, files) {
        int w, h;
        arrput(&images, (u8*)stbi_load(files[fi], &w, &h, &n, channels));
        Image_Size size;
        size.w = w;
        size.h = h;
        arrput(&image_sizes, size);
    }
    
    double total_ms = 0;
    u64 process_image_count = 0;
    u64 total_bytes_processed = 0;
    
    LOG_INFO("Warmup started!");
    for_n (warmup, 0) {
        allocator.reset();
        int w = image_sizes[0].w;
        int h = image_sizes[0].h;
        Process_Image_Result output = process_image(images[0], w, h, channels, &allocator);                	
    	bg_unused(output);
    }
    LOG_INFO("Warmup done!");

    for_n (_overloop, 100) {
        for_array (_i, images) {
			s64 start = start_prof();
		
		    allocator.reset();
            int w = image_sizes[_i].w;
            int h = image_sizes[_i].h;
            Process_Image_Result output = process_image(images[_i], w, h, channels, &allocator);
            // draw_rectangle(images[_i], w, h, channels, { 0,0,200,200 });
            process_image_count += 1;
            total_bytes_processed += (u64)w * h * channels;

            // char bf[1024];
            // memset(bf, 0, sizeof(bf));
            // sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
            // stbi_write_bmp(bf, w, h, channels, images[_i]);
            // 
            // sprintf(bf, "%s\\mini_output_%s", output_path, path_file_name(files[_i]));
            // stbi_write_bmp(bf, output.downsampled_mask.w, output.downsampled_mask.h, 1, output.downsampled_mask.mask);
            // printf("output target : %s\n", bf);
			total_ms += time_elapsed_ms(start);
			
            bg_unused(output);

        }
        printf("loop %llu\n", _overloop);
    }

    for_array(_i, images) {
        stbi_image_free(images[_i]);
    }
    arrfree(&images);
    printf("Custom implementation : total %.3fms\nav time %.3fms\nbytes processed : %.4fMB, throughput : %.4fMB/ms", total_ms, total_ms / process_image_count, total_bytes_processed / (1024.0 * 1024.0), (double)total_bytes_processed /(1024.0 * 1024.0) / total_ms);
    
    return 0;

}
   

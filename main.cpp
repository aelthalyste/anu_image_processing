#include <stdio.h>
#include "anu.h"

#define BG_IMPLEMENTATION
#include "bg.hpp"

#include "anu_stb_include.h"

struct Linear_Allocator;
struct Image_Load_Queue;
struct Image_Load_Result;
struct Image_Write_Order;

void  init_image_write_order(Image_Write_Order *order, void *image, s32 w, s32 h, s32 channels);
void  init_image_load_queue(Image_Load_Queue *queue, u64 queue_len, Linear_Allocator *linear_allocator);
void  free_image(Image_Load_Result *r);
void  load_image_background(Image_Load_Result *result, const char *path);
DWORD image_load_thread(void *arg);

// some debug stuff
static inline int64_t start_prof();
static inline int64_t get_perf_freq();

// time elapsed in ms
static inline double time_elapsed_ms(int64_t start);


struct Image_Load_Result {
    void *image;
    volatile s64 done;
    s32 w, h;
    s32 channels;
};

struct Image_Load_Thread_Params {
    char const **paths;
    Image_Load_Result **queue;

    u64 len;

    volatile s64 begin;
    volatile s64 end;

} global_image_loader;

struct Image_Write_Order {
    void *image;
    s32 w, h;
    s32 channels;
    s32 result;
};

struct Image_Writer_Thread {
    char const **paths;
    Image_Write_Order **queue;
    s64 begin, end;
    s64 len;
} global_image_writer;

struct Image_Load_Queue {
    Image_Load_Result *arr;
    s64 len;

    s64 begin;
    s64 end;
    s64 rem;

    void queue(char *path);
    Image_Load_Result consume();
};



void Image_Load_Queue::queue(char *path) {
    
    // wait till queue empties
    while(rem==0);
    --rem;
    
    load_image_background(&arr[end], path);
    end = (end + 1) % len;
}

Image_Load_Result Image_Load_Queue::consume() {
    while (arr[begin].done == 0);

    Image_Load_Result result = arr[begin];

    begin = (begin + 1) % len;
    ++rem;

    return result;
}


void init_image_load_queue(Image_Load_Queue *queue, u64 queue_len, Linear_Allocator *linear_allocator) {
    queue->len = (s64)queue_len;
    queue->arr = (Image_Load_Result *)linear_allocator->allocate(queue_len * sizeof(queue->arr[0]));
    queue->begin = 0;
    queue->end   = 0;
    queue->rem   = queue->len;
}  


void free_image(Image_Load_Result *r) {
    stbi_image_free(r->image);
}

void load_image_background(Image_Load_Result *result, const char *path) {
    zero_memory(result, sizeof(*result));

    s64 indx = global_image_loader.end;
    

    global_image_loader.paths[indx] = path;
    global_image_loader.queue[indx] = result;
    global_image_loader.end = (indx + 1) % global_image_loader.len;
}

DWORD image_load_thread(void *arg) {
    (void)(arg);

    memset(&global_image_loader, 0, sizeof(global_image_loader));

    global_image_loader.len   = 1024;
    global_image_loader.paths = (char const **)bg_calloc(global_image_loader.len, 8);
    global_image_loader.queue = (Image_Load_Result **)bg_calloc(global_image_loader.len, 8);
    
    BG_ASSERT(global_image_loader.queue);
    BG_ASSERT(global_image_loader.paths);

    for (;;) {
        
        if (global_image_loader.begin != global_image_loader.end) {
            s64 indx = global_image_loader.begin;
            Image_Load_Result *result = global_image_loader.queue[indx];
            s32 n;
            result->image = stbi_load(global_image_loader.paths[indx], &result->w, &result->h, &n, 4);
            result->channels = 4;
            result->done = 1;
            global_image_loader.begin = (indx + 1) % global_image_loader.len;
        } else {
            Sleep(1);
        }

    }

}

#if 1

template<typename T>
struct Dequeue {
    T *arr;
    s64 len;
    s64 begin, end;
    s64 rem;
};

template<typename T>
void init_dequeue(Dequeue<T> *o, s64 len, Linear_Allocator *linear_allocator) {
    o->begin = 0;
    o->end   = 0;
    o->rem   = len;
    o->len   = len;
    o->arr   = (T *)linear_allocator->allocate(o->len * sizeof(o->arr[0]));
}


void save_image_background(const char *path, Image_Write_Order *order) {
    s64 indx = global_image_writer.end;     
    
    global_image_writer.paths[indx] = path;
    global_image_writer.queue[indx] = order;
    global_image_writer.end = (indx + 1) % global_image_writer.len;
}


DWORD image_write_thread(void *arg) {

    global_image_writer.len   = 1024;
    global_image_writer.queue = (Image_Write_Order **)bg_calloc(global_image_writer.len, 8);
    global_image_writer.paths = (char const **)bg_calloc(global_image_writer.len, 8);

    for(;;) {
        
        if (global_image_writer.begin != global_image_writer.end) {
            s64 indx = global_image_writer.begin;
            Image_Write_Order *order = global_image_writer.queue[indx];
            const char * target_path = global_image_writer.paths[indx];
            stbi_write_bmp(target_path, order->w, order->h, order->channels, order->image);
            global_image_writer.begin = (indx + 1) % global_image_writer.len;
            order->result = 1;
        } else {
            Sleep(1);
        }

    }

}
#endif


void init_image_write_order(Image_Write_Order *order, void *image, s32 w, s32 h, s32 channels) {
    order->image = image;
    order->w = w;
    order->h = h;
    order->channels = channels;
    order->result = 0;
}


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
    HANDLE thread_handle = CreateThread(0, Megabyte(1), image_load_thread, NULL, 0, NULL);
    BG_ASSERT(thread_handle != INVALID_HANDLE_VALUE);
    CloseHandle(thread_handle);

    const char *input_path   = "C:\\W\\anu_tests\\data\\ttt\\";
    const char *output_path  = "C:\\W\\anu_tests\\output\\";
	bg_unused(output_path);
	bg_unused(input_path);

    auto files = get_file_paths_in_directory(input_path);

    u64 memory_size = 1024 * 1024 * 128;
    u8 *memory = (u8 *)_aligned_malloc(1024 * 1024 * 128, 64);
    Linear_Allocator allocator = init_linear_allocator(memory, memory_size, 16);


    double total_ms = 0;
    u64 process_image_count = 0;
    u64 total_bytes_processed = 0;

    Dequeue<Image_Load_Result> dq;
    init_dequeue(&dq, 1024, &allocator);
    
    Image_Load_Queue q;
    init_image_load_queue(&q, 8, &allocator);

    u64 fi = 0;
    u64 image_count = files.len;

    while (image_count) {
        
        s64 c = q.rem;
        while(c-- && fi < files.len) q.queue(files[fi++]);
        
        image_count--;

        Image_Load_Result r = q.consume();
        defer({free_image(&r);});
        
        s64 start = start_prof();
        
        Allocator_Mark mark = allocator.mark();
        defer({allocator.restore(mark);});

        Process_Image_Result output = process_image((u8 *)r.image, r.w, r.h, r.channels, &allocator);

        process_image_count   += 1;
        total_bytes_processed += (u64)r.w * r.h * r.channels;
        total_ms += time_elapsed_ms(start);

        char bf[1024];
        memset(bf, 0, sizeof(bf));
        
        //sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
        //stbi_write_bmp(bf, output.pure_mask.w, output.pure_mask.h, 1, output.pure_mask.mask);

        // sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
        // stbi_write_bmp(bf, output.pure_mask.w, output.pure_mask.h, 1, output.pure_mask.mask);

        sprintf(bf, "C:\\w\\anu_tests\\output\\output_rect_%llu.bmp", image_count);
        stbi_write_bmp(bf, r.w, r.h, r.channels, r.image);
    }

#if 0
    for_n (_i, files.len) {

        Image_Load_Result r = q.consume();
        defer({free_image(&r);});
        
		s64 start = start_prof();
        

	    Allocator_Mark mark = allocator.mark();
        defer({allocator.restore(mark);});

        Process_Image_Result output = process_image((u8 *)r.image, r.w, r.h, r.channels, &allocator);

        process_image_count   += 1;
        total_bytes_processed += (u64)r.w * r.h * r.channels;
        total_ms += time_elapsed_ms(start);

        #if 1
        char bf[1024];
        memset(bf, 0, sizeof(bf));
        
        //sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
        //stbi_write_bmp(bf, output.pure_mask.w, output.pure_mask.h, 1, output.pure_mask.mask);

        // sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
        // stbi_write_bmp(bf, output.pure_mask.w, output.pure_mask.h, 1, output.pure_mask.mask);

        sprintf(bf, "C:\\w\\anu_tests\\output\\output_rect_%s", path_file_name(files[_i]));
        stbi_write_bmp(bf, r.w, r.h, r.channels, r.image);
        #endif

        bg_unused(output);
    }
#endif

    printf("Custom implementation : total %.3fms\nav time %.3fms\nbytes processed : %.4fMB, throughput : %.4fMB/ms", total_ms, total_ms / process_image_count, total_bytes_processed / (1024.0 * 1024.0), (double)total_bytes_processed /(1024.0 * 1024.0) / total_ms);
    
    return 0;

}
   

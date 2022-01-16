#define Kilobyte(val) ((val)*1024ll)
#define Megabyte(val) (Kilobyte(val)*1024ll)
#define Gigabyte(val) (Megabyte(val)*1024ll)

#include <stdio.h>

#include "anu.h"
#include "profile.h"
#include "memory.h"
#include "job.h"

#include "anu_stb_include.h"

#include <Windows.h>
#include <intrin.h>

#define SDL_MAIN_HANDLED 1
#include "SDL.h"

#pragma comment(lib, "mincore")
#pragma comment(lib, "User32")



struct Image_Load_Result {
    void *image;
    volatile int64_t done;
    int32_t w, h;
    int32_t channels;
};
typedef struct Image_Load_Result Image_Load_Result;


struct Image_Load_Thread_Params {
    char const **paths;
    Image_Load_Result **queue;

    uint64_t len;

    volatile int64_t begin;
    volatile int64_t end;

} ;
typedef struct Image_Load_Thread_Params Image_Load_Thread_Params;
Image_Load_Thread_Params global_image_loader;


struct Image_Write_Order {
    void *image;
    int32_t w, h;
    int32_t channels;
    int32_t result;
};
typedef struct Image_Write_Order Image_Write_Order;

struct Image_Writer_Thread {
    char const **paths;
    Image_Write_Order **queue;
    int64_t begin, end;
    int64_t len;
};
typedef struct Image_Writer_Thread Image_Writer_Thread;
Image_Writer_Thread global_image_writer;


struct Image_Load_Queue {
    Image_Load_Result *arr;
    int64_t len;

    int64_t begin;
    int64_t end;
    int64_t rem;
};
typedef struct Image_Load_Queue Image_Load_Queue; 



void  init_image_write_order(Image_Write_Order *order, void *image, int32_t w, int32_t h, int32_t channels);
void  init_image_load_queue(Image_Load_Queue *queue, uint64_t queue_len, Linear_Allocator *linear_allocator);
void  free_image(Image_Load_Result *r);
void  load_image_background(Image_Load_Result *result, const char *path);

static int image_load_thread(void *arg);
int image_write_thread(void *arg);



static inline void image_load_queue_queue(Image_Load_Queue *q, char *path) {
    // wait till queue empties
    while(q->rem==0);
    
    --q->rem;
    
    load_image_background(&q->arr[q->end], path);
    q->end = (q->end + 1) % q->len;
}

static inline void image_load_queue_consume(Image_Load_Queue *q, Image_Load_Result *result) {
    while (q->arr[q->begin].done == 0);

    *result = q->arr[q->begin];

    q->begin = (q->begin + 1) % q->len;
    ++q->rem;
}

static inline int64_t image_load_queue_test_consume(Image_Load_Queue *q, Image_Load_Result *result) {
    if (q->arr[q->begin].done == 0)
        return 0;

    Assert(q->begin < q->len);
    *result = q->arr[q->begin];
    q->begin = (q->begin + 1) % q->len;
    ++q->rem;
    return 1; 
}


void init_image_load_queue(Image_Load_Queue *q, uint64_t queue_len, Linear_Allocator *linear_allocator) {
    q->len = (int64_t)queue_len;
    q->arr = (Image_Load_Result *)linear_allocate(linear_allocator, queue_len * sizeof(q->arr[0]));
    q->begin = 0;
    q->end   = 0;
    q->rem   = q->len;
}  


void free_image(Image_Load_Result *r) {
    stbi_image_free(r->image);
}

void load_image_background(Image_Load_Result *result, const char *path) {
    memset(result, 0, sizeof(*result));

    int64_t indx = global_image_loader.end;
    

    global_image_loader.paths[indx] = path;
    global_image_loader.queue[indx] = result;
    global_image_loader.end = (indx + 1) % global_image_loader.len;
}

static int image_load_thread(void *arg) {

    (void)(arg);

    memset(&global_image_loader, 0, sizeof(global_image_loader));

    global_image_loader.len   = 1024;
    global_image_loader.paths = (char const **)calloc(global_image_loader.len, 8);
    global_image_loader.queue = (Image_Load_Result **)calloc(global_image_loader.len, 8);
    
    Assert(global_image_loader.queue);
    Assert(global_image_loader.paths);

    for (;;) {
        
        if (global_image_loader.begin != global_image_loader.end) {
            int64_t indx = global_image_loader.begin;
            Image_Load_Result *result = global_image_loader.queue[indx];
            int32_t n;
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


void save_image_background(const char *path, Image_Write_Order *order) {
    int64_t indx = global_image_writer.end;     
    
    global_image_writer.paths[indx] = path;
    global_image_writer.queue[indx] = order;
    global_image_writer.end = (indx + 1) % global_image_writer.len;
}


int image_write_thread(void *arg) {
    
    (void)(arg);

    global_image_writer.len   = 1024;
    global_image_writer.queue = (Image_Write_Order **)calloc(global_image_writer.len, 8);
    global_image_writer.paths = (char const **)calloc(global_image_writer.len, 8);

    for(;;) {
        
        if (global_image_writer.begin != global_image_writer.end) {
            int64_t indx = global_image_writer.begin;
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


void init_image_write_order(Image_Write_Order *order, void *image, int32_t w, int32_t h, int32_t channels) {
    order->image = image;
    order->w = w;
    order->h = h;
    order->channels = channels;
    order->result = 0;
}


static inline uint64_t string_length(const char *s) {
    const char *b = s;
    for (; *s != 0; s++);
    return s - b;
}

static inline uint64_t string_append(char *str, const char *app) {
    uint64_t strlen = string_length(str);
    uint64_t applen = string_length(app);
    memcpy(str + strlen, app, applen);
    str[strlen + applen] = 0;
    return strlen + applen;
}


static inline char ** get_file_paths_in_directory(const char *dir, uint64_t *file_count, Linear_Allocator *allocator) {
    
    char wildcard_dir[280];
    memset(wildcard_dir, 0, sizeof(wildcard_dir));
    uint64_t dirlen = string_length(dir);
    dirlen++; // null termination

    Assert(sizeof(wildcard_dir) > dirlen);

    uint64_t file_cap = 1024 * 4;
    char ** result = linear_allocate(allocator, file_cap * 8);

    string_append(wildcard_dir, dir);
    string_append(wildcard_dir, "\\*");


    WIN32_FIND_DATAA FDATA;
    HANDLE FileIterator = FindFirstFileA(wildcard_dir, &FDATA);
    uint64_t fc = 0;

    Assert(FileIterator != INVALID_HANDLE_VALUE);

    if (FileIterator != INVALID_HANDLE_VALUE) {

        while (FindNextFileA(FileIterator, &FDATA) != 0) {

            //@NOTE(Batuhan): Do not search for sub-directories, skip folders.
            if (FDATA.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                continue;
            }

            // @Stupid +5 just to be sure there is enough room for null termination
            uint64_t fn_len = string_length(FDATA.cFileName);
            fn_len += dirlen + 5;

            char *fnbuffer = (char *)linear_allocate_zero(allocator, fn_len);

            string_append(fnbuffer, dir);
            string_append(fnbuffer, "\\");
            string_append(fnbuffer, FDATA.cFileName);

            result[fc++] = fnbuffer;
            Assert(fc < file_cap);
        }


    }
    else {
        printf("Cant iterate directory\n");
    }


    FindClose(FileIterator);

    *file_count = fc;

    return result;
}

int g_mousex;
int g_mousey;
int g_mouseleftclick;
int g_mouserightclick;



int main(int argc, char* argv[]) {
    (void)(argc);
    (void)(argv);

    uint64_t memory_size = 1024 * 1024 * 128;
    uint8_t *memory = VirtualAlloc(NULL, memory_size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    Linear_Allocator allocator;
    init_linear_allocator(&allocator, memory, memory_size, 16);

    //stbi_set_flip_vertically_on_load(1);
    //stbi_flip_vertically_on_write(true);

    SDL_Init(SDL_INIT_EVERYTHING);
    //SDL_Init(SDL_INIT_VIDEO);
    if (NULL == SDL_CreateThread(image_load_thread, "image load thread", (void *)NULL)) {
        fprintf(stderr, "Unable to create image load thread, reason : %s", SDL_GetError());
        return 0;
    }


    const char *input_path   = "C:\\W\\anu_tests\\data\\";
    const char *output_path  = "C:\\W\\anu_tests\\output\\";
    (void)(output_path);
    

    uint64_t file_count = 0;
    char **files = get_file_paths_in_directory(input_path, &file_count, &allocator);



    Image_Load_Queue q;
    init_image_load_queue(&q, 32, &allocator);

    uint64_t fi = 0;

    uint32_t window_width  = 1280;
    uint32_t window_height = 720;
    (void)(window_width);
    (void)(window_height);

#if 1
    SDL_Window* window  = SDL_CreateWindow("SDL pixels", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, window_width, window_height, SDL_WINDOW_SHOWN);
    SDL_Surface* screen = SDL_GetWindowSurface(window);
    SDL_Surface* pixels = SDL_CreateRGBSurfaceWithFormat(0, window_width, window_height, 32, SDL_PIXELFORMAT_RGBA32);


    Image_Load_Result loaded_images[4096];
    memset(loaded_images, 0, sizeof(loaded_images));
    int64_t lic = 0;

    uint64_t default_image_index = 0;
    uint8_t *default_image = linear_allocate_zero(&allocator, window_width * window_height * 4 * 2);
    memset(default_image, 0xbb, window_width * window_height * 4 * 2);
    loaded_images[default_image_index].image = default_image;
    loaded_images[default_image_index].w = window_width;
    loaded_images[default_image_index].w = window_height;
    loaded_images[default_image_index].done = 1;
    loaded_images[default_image_index].channels = 4;

    int64_t image_to_be_rendered = 0;
    int64_t save_this_image = 0;
    int64_t process_this_image = 0;

    uint64_t subsample_window_size = 32;

    Allocator_Mark allocator_mark;
    memset(&allocator_mark, 0, sizeof(allocator_mark));

    Process_Image_Result pi_result;
    memset(&pi_result, 0, sizeof(pi_result));

    for(;;) {
        
        linear_allocator_mark(&allocator, &allocator_mark);

        SDL_Event ev;
        while (SDL_PollEvent(&ev))
        {
            if (ev.type == SDL_QUIT)
                return 0;

            if (ev.type == SDL_KEYDOWN) {
                
                if (ev.key.keysym.scancode == SDL_SCANCODE_ESCAPE)
                    return 0;

                if (ev.key.keysym.scancode == SDL_SCANCODE_E) {
                    image_to_be_rendered++;
                    memset(&pi_result, 0, sizeof(pi_result));
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_Q) {
                    image_to_be_rendered--;
                    memset(&pi_result, 0, sizeof(pi_result));
                }

                if (ev.key.keysym.scancode == SDL_SCANCODE_S)
                    save_this_image = 1;
                if (ev.key.keysym.scancode == SDL_SCANCODE_W)
                    process_this_image = 1;
                
                if (ev.key.keysym.scancode == SDL_SCANCODE_1) {
                    process_this_image = 1;
                    if (subsample_window_size > 8)
                        subsample_window_size -= 8;
                }
                if (ev.key.keysym.scancode == SDL_SCANCODE_2) {
                    process_this_image = 1;
                    subsample_window_size += 8;
                }

            }


        }


        // queue new images and try to consume from queue
        {
            int64_t c = q.rem;
            while(c-- && fi < file_count)
                image_load_queue_queue(&q, files[fi++]);

            if (image_load_queue_test_consume(&q, &loaded_images[lic+1]))
                ++lic;
        }

    
        if (image_to_be_rendered > lic)
            image_to_be_rendered = image_to_be_rendered % lic;
        if (image_to_be_rendered < 0)
            image_to_be_rendered = lic;

        Image_Load_Result *li = &loaded_images[image_to_be_rendered];
        
        unsigned char* framebuffer = linear_allocate(&allocator, li->w * li->h * li->channels);
        int w = li->w;
        int h = li->h;
        
        memcpy(framebuffer, li->image, li->w * li->h * li->channels);

        if (process_this_image) {
            process_this_image = 0;

            Assert(subsample_window_size != 0);
                        
            pi_result= process_image((uint8_t *)framebuffer, li->w, li->h, li->channels, subsample_window_size, subsample_window_size, &allocator);            
        }


        for (int i = 0; i < pi_result.rc; ++i)
            draw_rectangle(framebuffer, w, h, 4, pi_result.rectangles[i]);


        SDL_LockSurface(pixels);
        {
            SDL_FillRect(pixels, NULL, 0);

            int yc = window_height>h ? h : (int)window_height;
            Assert(li->image);
            Assert(framebuffer);

            memset(pixels->pixels, 0, window_width * window_height * 4);           
            for(int y=0;y<yc;++y) {
                char *input  = &((char *)framebuffer)[y*w*4];
                char *output = &((char *)pixels->pixels)[y*pixels->pitch];
                
                int copy_size = w*4;
                if (copy_size > pixels->pitch)
                    copy_size = window_width*4;
                memcpy(output, input, copy_size);                
            }

        }
        SDL_UnlockSurface(pixels);



        SDL_BlitSurface(pixels, NULL, screen, NULL);
        SDL_UpdateWindowSurface(window);

        if (save_this_image) {
            stbi_write_bmp("actual_image.bmp", li->w, li->h, 4, li->image);
            stbi_write_bmp("frame_buffer.bmp", window_height, window_width, 4, framebuffer);
            save_this_image = 0;
        }

        linear_allocator_restore(&allocator, allocator_mark);
    }
#endif


#if 0
    uint64_t image_count = file_count;

    uint64_t process_image_count = 0;
    uint64_t total_bytes_processed = 0;

    double total_ms = 0;

    while (image_count) {
        
        int64_t c = q.rem;
        while(c-- && fi < file_count) image_load_queue_queue(&q, files[fi++]);
        
        image_count--;

        Image_Load_Result r;
        image_load_queue_consume(&q, &r);
        
        int64_t start = start_prof();
        
        Allocator_Mark mark;
        linear_allocator_mark(&allocator, &mark);

        Process_Image_Result output = process_image((uint8_t *)r.image, r.w, r.h, r.channels, &allocator);
        (void)(output);

        process_image_count   += 1;
        total_bytes_processed += (uint64_t)r.w * r.h * r.channels;
        total_ms += time_elapsed_ms(start, start_prof());

        char bf[1024];
        memset(bf, 0, sizeof(bf));
        
        //sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
        //stbi_write_bmp(bf, output.pure_mask.w, output.pure_mask.h, 1, output.pure_mask.mask);

        // sprintf(bf, "C:\\w\\anu_tests\\output\\output_%s", path_file_name(files[_i]));
        // stbi_write_bmp(bf, output.pure_mask.w, output.pure_mask.h, 1, output.pure_mask.mask);

        sprintf(bf, "C:\\w\\anu_tests\\output\\output_rect_%llu.bmp", image_count);
        stbi_write_bmp(bf, r.w, r.h, r.channels, r.image);

        free_image(&r);
        linear_allocator_restore(&allocator, mark);

    }
    printf("Custom implementation : total %.3fms\nav time %.3fms\nbytes processed : %.4fMB, throughput : %.4fMB/ms", total_ms, total_ms / process_image_count, total_bytes_processed / (1024.0 * 1024.0), (double)total_bytes_processed /(1024.0 * 1024.0) / total_ms);

#endif

    
    return 0;

}
   

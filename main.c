#define Kilobyte(val) ((val)*1024ll)
#define Megabyte(val) (Kilobyte(val)*1024ll)
#define Gigabyte(val) (Megabyte(val)*1024ll)

#include <stdio.h>
#include "anu.h"

#include "anu_stb_include.h"
#include <Windows.h>
#include <intrin.h>

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
DWORD image_load_thread(void *arg);

// some debug stuff
static inline int64_t start_prof(void);
static inline int64_t get_perf_freq(void);

// time elapsed in ms
static inline double time_elapsed_ms(int64_t start);



void init_job_list(Job_List *jl, uint64_t queue_cap, Linear_Allocator *allocator) {
    jl->jobs = linear_allocate(allocator, sizeof(jl->jobs[0]) * queue_cap);
    jl->len = queue_cap;
    jl->rem = jl->len;
    jl->begin = 0;
    jl->end   = 0;
    jl->terminate_all = 0;
    InitializeCriticalSection(&jl->add_job_section);
    InitializeCriticalSection(&jl->consume_job_section);
}



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

DWORD image_load_thread(void *arg) {

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


DWORD image_write_thread(void *arg) {
    
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


static inline int64_t start_prof(void){
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return li.QuadPart;
}

static inline int64_t get_perf_freq(void){
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
    
    char wildcard_dir[280] ={};
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


void queue_job(Job_List *jl, DWORD (*proc)(void *arg), void *argument) {
    
    uint64_t st = start_prof();
    EnterCriticalSection(&jl->add_job_section);

    Assert(proc);

    while (jl->rem == 0)
        yield_execution();

    
    int64_t indx = jl->end;
    jl->jobs[indx].proc     = proc;
    jl->jobs[indx].argument = argument;
    jl->end = (indx + 1) % jl->len;
    
    InterlockedDecrement64(&jl->rem);
    LeaveCriticalSection(&jl->add_job_section);    
    double elapsed = time_elapsed_ms(st);
    // printf("q took :  %.4f\n", elapsed);

}

void synchronize_jobs(Job_List *jl) {
    while (jl->rem != jl->len) {
        yield_execution();
    }
    Assert(jl->rem == jl->len);
}

DWORD job_poll_thread(void *arg) {

    Job_List *jl = (Job_List *)arg;
    for (;;) {

        uint32_t consumed = 0;
        Job job;

        if (jl->terminate_all)
            break;

        uint64_t consume_st = start_prof();

        EnterCriticalSection(&jl->consume_job_section);
        if (jl->rem != jl->len) {
            consumed = 1; 
            job = jl->jobs[jl->begin];
            Assert(job.proc != NULL);

            jl->begin = (jl->begin + 1) % jl->len;
            InterlockedIncrement64(&jl->rem);
        }
        LeaveCriticalSection(&jl->consume_job_section);

        if (consumed != 0) {

            uint64_t st = start_prof();
            double elapsed_consume = time_elapsed_ms(consume_st);

            job.proc(job.argument);
            double elapsed_proc = time_elapsed_ms(st);

            // if a task runs for lower than 100microseconds, reconsider queuing it.
            Assert(elapsed > 0.1);
            // printf("Consume took %.4f, proc took %.4f\n", elapsed_consume, elapsed_proc);
        }
        else 
            yield_execution();

    }

    return 0;

}

DWORD mystupidthread(void *arg) {
    // Sleep(1);
    return 0;
}

int main() {
    
    uint64_t memory_size = 1024 * 1024 * 128;
    uint8_t *memory = (uint8_t *)_aligned_malloc(1024 * 1024 * 128, 64);
    Linear_Allocator allocator;
    init_linear_allocator(&allocator, memory, memory_size, 16);

    //stbi_set_flip_vertically_on_load(true);
    //stbi_flip_vertically_on_write(true);
    HANDLE thread_handle = CreateThread(0, Megabyte(1), image_load_thread, NULL, 0, NULL);
    Assert(thread_handle != INVALID_HANDLE_VALUE);
    CloseHandle(thread_handle);


    Job_List *jl = linear_allocate(&allocator, sizeof(*jl));
    init_job_list(jl, 1024 * 4, &allocator);


    for (uint64_t i = 0; i < 16; ++i) {
        HANDLE thread_handle = CreateThread(0, Megabyte(1), job_poll_thread, jl, 0, NULL);
        Assert(thread_handle != INVALID_HANDLE_VALUE);
        CloseHandle(thread_handle);        
    }

    uint64_t st = start_prof();
    for (uint64_t i = 0; i < 1024 * 4; ++i) {
        queue_job(jl, mystupidthread, NULL);
    }

    double elapsed = time_elapsed_ms(st);

    synchronize_jobs(jl);
    
    printf("elapsed %.4f ms, us per q : %f\n", elapsed, elapsed / 1024 * 1000.0);
    return 0;

    const char *input_path   = "C:\\W\\anu_tests\\customtest\\";
    const char *output_path  = "C:\\W\\anu_tests\\output\\";


    uint64_t file_count = 0;
    char **files = get_file_paths_in_directory(input_path, &file_count, &allocator);



    double total_ms = 0;
    uint64_t process_image_count = 0;
    uint64_t total_bytes_processed = 0;

    
    Image_Load_Queue q;
    init_image_load_queue(&q, 8, &allocator);

    uint64_t fi = 0;
    uint64_t image_count = file_count;

    while (image_count) {
        
        int64_t c = q.rem;
        while(c-- && fi < file_count) image_load_queue_queue(&q, files[fi++]);
        
        image_count--;

        Image_Load_Result r;
        image_load_queue_consume(&q, &r);
        
        int64_t start = start_prof();
        
        Allocator_Mark mark;
        linear_allocator_mark(&allocator, &mark);

        Process_Image_Result output = process_image(jl, (uint8_t *)r.image, r.w, r.h, r.channels, &allocator);

        process_image_count   += 1;
        total_bytes_processed += (uint64_t)r.w * r.h * r.channels;
        total_ms += time_elapsed_ms(start);

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
    
    return 0;

}
   

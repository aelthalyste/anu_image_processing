#include "job.h"
#include "memory.h"
#include "profile.h"


DWORD job_poll_thread(void *arg) {

    Job_List *jl = (Job_List *)arg;
    for (;;) {

        uint32_t consumed = 0;
        Job job;
        memset(&job, 0, sizeof(job));
        if (jl->terminate_all)
            break;

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
            Assert(job.proc);

            uint64_t st = start_prof();
            job.proc(job.argument);
            double elapsed_proc = time_elapsed_ms(st, start_prof());

            // if a task runs for lower than 50microseconds, reconsider queuing it.
            Assert(elapsed_proc > 0.05);
            // printf("proc took %.4f\n", elapsed_proc);
        }
        // else 
        //     yield_execution();

    }

    return 0;

}


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


void queue_job(Job_List *jl, DWORD (*proc)(void *arg), void *argument) {
    Assert(proc);

    EnterCriticalSection(&jl->add_job_section);
    while (jl->rem == 0);
        // yield_execution();

    int64_t indx = jl->end;
    jl->jobs[indx].proc     = proc;
    jl->jobs[indx].argument = argument;
    jl->end = (indx + 1) % jl->len;
    
    InterlockedDecrement64(&jl->rem);
    LeaveCriticalSection(&jl->add_job_section);    
}

void synchronize_jobs(Job_List *jl) {
    while (jl->rem != jl->len)
        yield_execution();
    
    Assert(jl->rem == jl->len);
}

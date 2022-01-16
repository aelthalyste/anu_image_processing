#pragma once

#include <Windows.h>
#include "stdint.h"
#include "memory.h"
#include "util.h"

#define yield_execution() SwitchToThread();


typedef struct 
{
    DWORD (*proc)(void *arg);
    void *argument;
} Job;


typedef struct {
    Job *jobs;
    volatile int64_t len;
    volatile int64_t begin;
    volatile int64_t end;
    volatile int64_t rem;
    volatile int32_t terminate_all;
    CRITICAL_SECTION add_job_section;
    CRITICAL_SECTION consume_job_section;
} Job_List;


void init_job_list(Job_List *jl, uint64_t queue_cap, Linear_Allocator *allocator);
void queue_job(Job_List *jl, DWORD (*proc)(void *arg), void *argument);
void synchronize_jobs(Job_List *jl);
DWORD job_poll_thread(void *arg);


DWORD job_poll_thread(void *arg);


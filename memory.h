#pragma once

#include <string.h>
#include <malloc.h>
#include "stdint.h"
#include "util.h"


typedef struct 
{
    void *memory;
    uint64_t size;
    uint64_t used;
    uint64_t aligment;
} Linear_Allocator;

typedef struct 
{
    uint64_t internal_mark;
} Allocator_Mark;



static inline void init_linear_allocator(Linear_Allocator *allocator, void *memory, uint64_t memory_size, uint64_t aligment) {
    allocator->memory   = memory;
    allocator->size     = memory_size;
    allocator->aligment = aligment ? aligment : 16;
    allocator->used     = 0;
}

static inline void * linear_allocate_aligned(Linear_Allocator *allocator, uint64_t s, uint64_t al) {
    uint64_t al_bonus = (uint64_t)((uint8_t *)allocator->memory + allocator->used) % al;
    if (al_bonus)
        al_bonus = al - al_bonus;
    void *result = NULL;
    if (s + al_bonus + allocator->used < allocator->size) {
        result = (uint8_t *)allocator->memory + allocator->used + al_bonus;
        allocator->used += al_bonus + s;
    }
    Assert(((int64_t)result % al) == 0);
    return result;
}

static inline void * linear_allocate(Linear_Allocator *allocator, uint64_t s) {
    return linear_allocate_aligned(allocator, s, allocator->aligment);   
}

static inline void * linear_allocate_zero_aligned(Linear_Allocator *allocator, uint64_t s, uint64_t al) {
    void *result = linear_allocate_aligned(allocator, s, al);
    memset(result, 0, s);
    return result;
}

static inline void * linear_allocate_zero(Linear_Allocator *allocator, uint64_t s) {
    void *result = linear_allocate(allocator, s);
    memset(result, 0, s);
    return result;         
}


static inline void linear_allocator_mark(Linear_Allocator *allocator, Allocator_Mark *mark) {
    mark->internal_mark = allocator->used;
}

static inline void linear_allocator_restore(Linear_Allocator *allocator, Allocator_Mark mark) {
    allocator->used = mark.internal_mark;
}

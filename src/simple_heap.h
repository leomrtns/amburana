#ifndef SIMPLE_HEAP_H
#define SIMPLE_HEAP_H
#include <biomcmc.h>

/* max heap structure for storing smallest hash values 
 * code inspired by https://gist.github.com/vgoel30/5d81e6abf9464930c1e126dab04d5be3  */

typedef struct heap64_struct* heap64;

struct heap64_struct {
  uint64_t *hash;
  int heap_size, n; // size allocated to heap vector, and n=currently existing elements 
};

heap64 new_heap64 (int heap_size);
void del_heap64 (heap64 h64);
uint64_t heap64_get_maximum (heap64 h64);
uint64_t heap64_remove_maximum (heap64 h64);
void heap64_insert (heap64 h64, uint64_t h); 
void heap64_bubble_down (heap64 h64, int p);
void heap64_bubble_up (heap64 h64, int index);
void heap64_finalise_heap_pop (heap64 h64); // DO NOT USE: qsort is much faster
void heap64_finalise_heap_qsort (heap64 h64);

#endif

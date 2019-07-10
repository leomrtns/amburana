/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file simple_heap.h 
 *  \brief Max heap structure for storing smallest hash values, and for creating graph neighbourhood 
 * code inspired by https://gist.github.com/vgoel30/5d81e6abf9464930c1e126dab04d5be3  
 */
#ifndef _simple_heap_h_ 
#define _simple_heap_h_ 
#include <biomcmc.h>

typedef struct heap64_struct* heap64;
typedef struct heap_pqueue_struct* heap_pqueue;

struct heap64_struct {
  uint64_t *hash;
  int heap_size, n; // size allocated to heap vector, and n=currently existing elements 
};

typedef struct {
  int id; 
  double priority; 
  void *v; // extra information, currently NULL 
} hpq_item_struct;

struct heap_pqueue_struct {
  hpq_item_struct *item; 
  int heap_size, n;
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

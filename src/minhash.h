/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file minhash.h 
 *  \brief heap minhash, b-bit (k-partition) minhash sketches for DNA sequences
 */
#ifndef _amburana_minhash_h_
#define _amburana_minhash_h_

#include "simple_heap.h" 

typedef struct sketch_set_struct* sketch_set;  // collection of minhashes etc for one kmerhash 
typedef struct heap_minhash_sketch_struct* heap_minhash_sketch;  // only one that needs maxheap
typedef struct bbit_minhash_sketch_struct* bbit_minhash_sketch;  // also can use distinct hashes instead of x/mod

struct sketch_set_struct
{
  uint8_t n_heap_mh, n_bbit_mh, n_distances;
  uint8_t i_heap[16], i_bbit[16]; // id of hash elements
  double *dist;
  int heap_mh_size, bbit_mh_bits;
  heap_minhash_sketch *heap_mh;
  bbit_minhash_sketch *bbit_mh;
  kmerhash kmer;
};

struct heap_minhash_sketch_struct
{
  int sketch_size;
  heap64 sketch;
  kmerhash kmer; // we need some constants defined here (kmer sizes)
};

struct bbit_minhash_sketch_struct
{ 
  int sketch_size, n_bits;
  uint64_t suffix_mask, *sketch; // each vector will have the min value over all from bin
  kmerhash kmer;
};
/*! \brief do not allocate vectors, just set sizes. Usefull if other functions need to know them */
sketch_set new_sketch_set_bare_numbers_only (kmerhash kmer, int heap_mh_size, int bbit_mh_bits);
sketch_set new_sketch_set (kmerhash kmer, int heap_mh_size, int bbit_mh_bits);
void del_sketch_set (sketch_set sset);
sketch_set new_sketch_set_from_dna (char *dna, size_t dna_length, kmerhash kmer, int heap_mh_size, int bbit_mh_bits);
void sketch_set_compare (sketch_set ss1, sketch_set ss2);

#endif
/*
struct cm_sketch_struct
{
  int size, count;
  uint32_t mod;
  int **freq;
};
*/

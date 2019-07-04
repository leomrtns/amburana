/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file minhash.h 
 *  \brief MinHash, b-bit (k-partition) sketches for DNA sequences
 */
#ifndef _amburana_minhash_h_
#define _amburana_minhash_h_

#include "simple_heap.h" 

typedef struct minhash_struct* minhash;   // only one that needs maxheap
typedef struct onephash_struct* onephash; // also can use distinct hashes instead of x/mod

struct minhash_struct
{
  int sketch_size, n_sketch;
  heap64 *sketch;
  kmerhash kmer; // we need some constants defined here (kmer sizes)
};

struct onephash_struct
{ 
  int sketch_size, n_sketch, n_bits;
  uint64_t suffix_mask, **sketch; // each vector will have the min value over all from bin
  uint64_t precision_mask;
  kmerhash kmer;
};

minhash new_minhash_from_dna (char *dna, size_t dna_length, int sketch_size, bool dense);
void del_minhash (minhash mh);
void compare_minhash (minhash mh1, minhash mh2, double *distance);

onephash new_onephash_from_dna (char *dna, size_t dna_length, int n_bits, bool dense);
void del_onephash (onephash oph);
void compare_onephash (onephash oh1, onephash oh2, double *distance);

#endif
/*
struct cm_sketch_struct
{
  int size, count;
  uint32_t mod;
  int **freq;
};
*/

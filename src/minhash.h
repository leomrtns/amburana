/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file minhash.h 
 *  \brief MinHash sketches for DNA sequences
 */
#ifndef _amburana_minhash_h_
#define _amburana_minhash_h_

#include "simple_heap.h" 

typedef struct minhash_struct* minhash;

struct minhash_struct
{
  int sketch_size;
  heap64 *sketch;
};

minhash new_minhash (int sketch_size);
void del_minhash (minhash mh);
minhash new_minhash_from_dna (char *dna, int dna_length, int sketch_size);
void compare_minhashes (minhash mh1, minhash mh2, double *result);

#endif

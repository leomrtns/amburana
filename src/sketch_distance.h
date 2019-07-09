/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file sketch_distance.h 
 *  \brief feed sketch set distances into distance_generator_struct:: 
 */
#ifndef _amburana_sketch_distance_h_
#define _amburana_sketch_distance_h_

#include "minhash.h" 

typedef struct sketch_distance_gen_struct* sketch_distance_gen;

struct sketch_distance_gen_struct
{
  alignment aln;
  sketch_set *sset, sset_zero;  /*!< \brief sset[] are actual kmer sets, sset_zero will have sketch sizes */
  double secs_sketches, secs_distances;
  distance_generator generator;
  int ref_counter;
};

sketch_distance_gen new_sketch_distance_gen (alignment aln, kmerhash kmer, int heap_mh_size, int bbit_mh_bits);
void del_sketch_distance_gen (sketch_distance_gen sd);
distance_generator new_distance_generator_from_sketch_set (sketch_distance_gen sd);

#endif

/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "sketch_distance.h"

void sketch_distance_gen_wrapper (void *data, int sample1, int sample2, double *result);

sketch_distance_gen
new_sketch_distance_gen (alignment aln, kmerhash kmer, int heap_mh_size, int bbit_mh_bits)
{ 
  int i;
  sketch_distance_gen sd = (sketch_distance_gen) biomcmc_malloc (sizeof (struct sketch_distance_gen_struct));
  sd->aln = aln; sd->aln->ref_counter++;
  sd->sset_zero = new_sketch_set_bare_numbers_only (kmer, heap_mh_size, bbit_mh_bits);
  sd->sset = (sketch_set*) biomcmc_malloc (aln->ntax * sizeof (sketch_set));
  for (i=0; i < aln->ntax; i++) sd->sset[i] = NULL;
  sd->generator = NULL; // this function is called before new_distance_generator(); link is done below; 
  sd->secs_sketches = sd->secs_distances = 0.;
  return sd;
}

void
del_sketch_distance_gen (sketch_distance_gen sd)
{
  int i;
  if (!sd) return;
  if (--sd->ref_counter) return;
  for (i = sd->aln->ntax-1; i >= 0; i--) del_sketch_set (sd->sset[i]);
  del_sketch_set (sd->sset_zero);
  del_alignment (sd->aln);
  del_distance_generator (sd->generator);
  free (sd);
}

distance_generator
new_distance_generator_from_sketch_set (sketch_distance_gen sd)
{
  distance_generator dg = new_distance_generator (sd->aln->ntax, sd->sset_zero->n_distances); 
  distance_generator_set_function_data (dg, &sketch_distance_gen_wrapper, (void*) sd); // sd will have extra data needed to calc distances
  sd->generator = dg; dg->ref_counter++;
  return dg;
}

void
sketch_distance_gen_wrapper (void *data, int sample1, int sample2, double *result)
{
  int k;
  clock_t time0, time1;
  sketch_distance_gen sd = (sketch_distance_gen) data;
  /* we only extract kmer set when distance is needed for the first time */
  time0 = clock ();
#pragma omp critical
   {
    if (!sd->sset[sample1]) sd->sset[sample1] = new_sketch_set_from_dna (sd->aln->character->string[sample1], sd->aln->character->nchars[sample1], 
                                                                         sd->sset_zero->kmer, sd->sset_zero->heap_mh_size, sd->sset_zero->bbit_mh_bits);
    if (!sd->sset[sample2]) sd->sset[sample2] = new_sketch_set_from_dna (sd->aln->character->string[sample2], sd->aln->character->nchars[sample2],
                                                                         sd->sset_zero->kmer, sd->sset_zero->heap_mh_size, sd->sset_zero->bbit_mh_bits);
   }
  time1 = clock (); sd->secs_sketches += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;

  sketch_set_compare (sd->sset[sample1], sd->sset[sample2]); // will store distances in its own vector
  time1 = clock (); sd->secs_distances += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); 

  for (k = 0; k < sd->generator->n_distances; k++) result[k] = sd->sset[sample1]->dist[k]; // we hope distance_generator calling is same as sd->generator...
  return;
}




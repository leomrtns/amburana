/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "minhash.h"
// TODO: exact kmer counter and countmin sketches
#define xFULL64 0xffffffffffffffffUL // max possible uint64_t

heap_minhash_sketch new_heap_minhash_sketch (kmerhash kmer, int sketch_size);
void del_heap_minhash_sketch (heap_minhash_sketch minh);
void heap_minhash_sketch_insert (heap_minhash_sketch minh, uint64_t hash, size_t location);
void heap_minhash_sketch_finalise (heap_minhash_sketch minh);
void heap_minhash_sketch_distance (heap_minhash_sketch mh1, heap_minhash_sketch mh2, int idx_in_kmer_params, double *w_dist, double *u_dist);

bbit_minhash_sketch new_bbit_minhash_sketch (kmerhash kmer, int n_bits);
void del_bbit_minhash_sketch (bbit_minhash_sketch oph);
void bbit_minhash_sketch_insert (bbit_minhash_sketch oph, uint64_t hash);
double bbit_minhash_sketch_distance (bbit_minhash_sketch oh1, bbit_minhash_sketch oh2, int idx_in_kmer_params);

sketch_set  // do not allocate vectors, usefull if other functions need to know sizes
new_sketch_set_bare_numbers_only (kmerhash kmer, int heap_mh_size, int bbit_mh_bits)
{
  sketch_set sset = (sketch_set) biomcmc_malloc (sizeof (struct sketch_set_struct));
  sset->kmer = kmer; sset->kmer->ref_counter++;
  sset->n_heap_mh = kmer->n_hash/2; // this division is hardcoded, I could change it in the future
  sset->n_bbit_mh = kmer->n_hash - sset->n_heap_mh;
  sset->heap_mh_size = heap_mh_size;
  sset->bbit_mh_bits = bbit_mh_bits;
  sset->heap_mh = NULL;
  sset->bbit_mh = NULL; 
  sset->n_distances = sset->n_heap_mh + sset->n_bbit_mh; // vector size, but we can store 2x if we remember that we have i and j from d(i,j)
  sset->dist = NULL; // (cont. above) and thus i has n_distances and j as well. Very ugly hack, but as long as sketch_distance_gen understands...
  return sset;
}

sketch_set
new_sketch_set (kmerhash kmer, int heap_mh_size, int bbit_mh_bits)
{
  int i;
  sketch_set sset = new_sketch_set_bare_numbers_only (kmer, heap_mh_size, bbit_mh_bits);
  sset->heap_mh = (heap_minhash_sketch*) biomcmc_malloc (sset->n_heap_mh * sizeof (heap_minhash_sketch));
  sset->bbit_mh = (bbit_minhash_sketch*) biomcmc_malloc (sset->n_bbit_mh * sizeof (bbit_minhash_sketch));
  for (i=0; i < sset->n_heap_mh; i++) {
    sset->heap_mh[i] = new_heap_minhash_sketch (kmer, heap_mh_size);
    sset->i_heap[i] = i;
  }
  for (i=0; i < sset->n_bbit_mh; i++) {
    sset->bbit_mh[i] = new_bbit_minhash_sketch (kmer, bbit_mh_bits);
    sset->i_bbit[i] = i + sset->n_heap_mh;
  }
  sset->dist = (double*) biomcmc_malloc (sset->n_distances * sizeof (double));
  return sset;
}

void
del_sketch_set (sketch_set sset)
{
  int i;
  if (!sset) return;
  for (i = sset->n_heap_mh-1; i >= 0; i--) del_heap_minhash_sketch (sset->heap_mh[i]);
  for (i = sset->n_bbit_mh-1; i >= 0; i--) del_bbit_minhash_sketch (sset->bbit_mh[i]);
  if (sset->heap_mh) free (sset->heap_mh);
  if (sset->bbit_mh) free (sset->bbit_mh);
  if (sset->dist) free (sset->dist);
  del_kmerhash (sset->kmer);
  free (sset);
}

sketch_set
new_sketch_set_from_dna (char *dna, size_t dna_length, kmerhash kmer, int heap_mh_size, int bbit_mh_bits)
{
  int i;
  sketch_set sset = new_sketch_set (kmer, heap_mh_size, bbit_mh_bits);
  link_kmerhash_to_dna_sequence (sset->kmer, dna, dna_length); 

  while ( kmerhash_iterator (kmer)) { // all magic happens here
    for (i = 0; i < sset->n_heap_mh; i++) heap_minhash_sketch_insert (sset->heap_mh[i], kmer->hash[ sset->i_heap[i] ], kmer->i);
    for (i = 0; i < sset->n_bbit_mh; i++) bbit_minhash_sketch_insert (sset->bbit_mh[i], kmer->hash[ sset->i_bbit[i] ]);
  }

  for (i = 0; i < sset->n_heap_mh; i++) heap_minhash_sketch_finalise (sset->heap_mh[i]);
  return sset;
}

void
sketch_set_compare (sketch_set ss1, sketch_set ss2)
{
  int i;
  if (ss1->kmer != ss2->kmer) biomcmc_error ("Sketch sets were initialised from different kmerhash_struct::");
  if (ss1->n_heap_mh != ss2->n_heap_mh) biomcmc_error ("Sketch sets with distinct number of heap_minhash_sketch sizes");
  if (ss1->n_bbit_mh != ss2->n_bbit_mh) biomcmc_error ("Sketch sets with distinct number of bbit_minhash_sketch sizes");
  
  for (i = 0; i < ss1->n_heap_mh; i++) // store 2 dists, weighted and unweighted, into two vectors
    heap_minhash_sketch_distance (ss1->heap_mh[i],ss2->heap_mh[i], ss1->i_heap[i], &(ss1->dist[i]), &(ss2->dist[i]));
  for (i = 0; i < ss1->n_bbit_mh; i++) 
    ss1->dist[i + ss1->n_heap_mh] = ss2->dist[i + ss2->n_heap_mh] = bbit_minhash_sketch_distance (ss1->bbit_mh[i],ss2->bbit_mh[i], ss1->i_bbit[i]);
  return;
}

heap_minhash_sketch
new_heap_minhash_sketch (kmerhash kmer, int sketch_size)
{
  heap_minhash_sketch minh = (heap_minhash_sketch) biomcmc_malloc (sizeof (struct heap_minhash_sketch_struct));
  if (sketch_size < 8)    sketch_size = 8; // ideal 16~128 
  if (sketch_size > 8192) sketch_size = 8192; 
  minh->sketch_size = sketch_size;
  minh->sketch = new_heap_hash64 (sketch_size);
  minh->kmer = kmer;
  minh->kmer->ref_counter++;
  return minh;
}

void
del_heap_minhash_sketch (heap_minhash_sketch minh)
{
  if (!minh) return;
  del_heap_hash64 (minh->sketch);
  del_kmerhash (minh->kmer);
  free (minh);
}

void // void is to allow for overloading (function pointers)
heap_minhash_sketch_insert (heap_minhash_sketch minh, uint64_t hash, size_t location)
{
  hpq_item item = {.id = (int) location, .freq = 1, .hash = hash};
  heap_hash64_insert (minh->sketch, item);
}

void
heap_minhash_sketch_finalise (heap_minhash_sketch minh)
{
  heap_hash64_finalise_heap_qsort (minh->sketch);
}

void
heap_minhash_sketch_distance (heap_minhash_sketch mh1, heap_minhash_sketch mh2, int idx_in_kmer_params, double *w_dist, double *u_dist)
{
  int j1, j2, n1, n2, common = 0, compared = 0, cosine = 0;
  hpq_item *h1, *h2;
  double r_u, r_w;

  if (mh1->sketch->heap_size != mh2->sketch->heap_size) 
    biomcmc_error ("Distinct size sketches: comparable in theory but more likely a programmer mistake");

  h1 = mh1->sketch->item; n1 = mh1->sketch->n;
  h2 = mh2->sketch->item; n2 = mh2->sketch->n;
  for (j1 = j2 = 0; (j1 < n1) && (j2 < n2) && compared < mh1->sketch->heap_size; compared++) { // both are in decreasing order
    if      (h1[j1].hash > h2[j2].hash) j1++; 
    else if (h1[j1].hash < h2[j2].hash) j2++; 
    else { cosine += (h1[j1].freq * h2[j2].freq); common++; j1++; j2++; } // same hash in both
  }
  if (compared < mh1->sketch->heap_size) compared += (n1 - j1 + n2 - j2); // try to complete the union operation (following Mash)
  if (compared > mh1->sketch->heap_size) compared = mh1->sketch_size;
  if (common == compared) *u_dist = 0.;
  if (common == 0) *u_dist = 1.;
  r_u = (double) (common) / (double) (compared); // unweighted (Jaccard similarity)
  *u_dist = -log (2. * r_u/(1. + r_u)) / (double) (mh1->kmer->p->size[idx_in_kmer_params]);
  r_w = (double)(cosine) / (mh1->sketch->sqrt_sum * mh2->sketch->sqrt_sum); // cosine similarity (cannot be < 0 since we use freq)
  *w_dist = -log (2. * r_w/(1. + r_w)) / (double) (mh1->kmer->p->size[idx_in_kmer_params]);
}

bbit_minhash_sketch
new_bbit_minhash_sketch (kmerhash kmer, int n_bits)
{
  int i, j;
  bbit_minhash_sketch oph = (bbit_minhash_sketch) biomcmc_malloc (sizeof (struct bbit_minhash_sketch_struct));
  if (n_bits < 4)  n_bits = 4; // good values seem to be  8~10 
  if (n_bits > 13) n_bits = 13;
  oph->n_bits = n_bits;
  oph->suffix_mask = 0UL;
  for (i = 0; i < oph->n_bits; i++) oph->suffix_mask |= (1UL << i); 

  oph->sketch_size = 1UL << oph->n_bits; // vector size will be 2^n_bits
  oph->sketch = (uint64_t*) biomcmc_malloc (oph->sketch_size * sizeof (uint64_t));
  for (j = 0; j < oph->sketch_size; j++) oph->sketch[j] = xFULL64; // max possible uint64_t

  oph->kmer = kmer; oph->kmer->ref_counter++;
  return oph;
}

void
del_bbit_minhash_sketch (bbit_minhash_sketch oph)
{
  if (!oph) return;
  if (oph->sketch) free (oph->sketch); 
  del_kmerhash (oph->kmer);
  free (oph);
}

void
bbit_minhash_sketch_insert (bbit_minhash_sketch oph, uint64_t hash)
{
  uint64_t prefix = hash >> oph->n_bits;      // hash with (64 - n_bits) bits
  uint32_t suffix = hash &  oph->suffix_mask; // int value of least sig n_bits of hash
    if (oph->sketch[suffix] > prefix) oph->sketch[suffix] = prefix;
}

double
bbit_minhash_sketch_distance (bbit_minhash_sketch oh1, bbit_minhash_sketch oh2, int idx_in_kmer_params)
{
  int j, compared = 0, common = 0; 
  double r;
  if (oh1->sketch_size != oh2->sketch_size) biomcmc_error ("can't compare one-perm hashes of distinct sketch sizes");
  for (j = 0; j < oh1->sketch_size; j++) 
    if ((oh1->sketch[j] != xFULL64) && (oh2->sketch[j] != xFULL64)) {
      compared++;
      if (oh1->sketch[j] == oh2->sketch[j]) common++;
    }
  if (common == compared) return 0.;
  if (common == 0) return 1.;
  r = (double) (common) / (double) (compared);
  return -log (2. * r/(1. + r)) / (double) (oh1->kmer->p->size[idx_in_kmer_params]);
}


/*
heap_minhash
new_heap_minhash_from_dna (char *dna, size_t dna_length, int sketch_size, bool dense)
{
  int i;
  kmerhash kmer = new_kmerhash (dense); // false = 4 bits per site (o.w. 2 bits) 
  link_kmerhash_to_dna_sequence (kmer, dna, dna_length); // not ideal, since we can create ONE kmerhash for whole analysis 
  heap_minhash_sketch minh = new_heap_minhash_sketch (sketch_size, kmer->n_hash); // use all hashes available by kmerhash
  minh->kmer = kmer;

  while ( kmerhash_iterator (minh->kmer)) {
    for (i = 0; i < minh->n_sketch; i++) heap64_insert (minh->sketch[i], minh->kmer->hash[i]);
  }

  for (i = 0; i < minh->n_sketch; i++) heap64_finalise_heap_qsort (minh->sketch[i]);
  return minh;
}
bbit_minhash
new_bbit_minhash_from_dna (char *dna, size_t dna_length, int n_bits, bool dense)
{
  int i;
  uint64_t prefix = 0UL;
  uint32_t suffix = 0UL;

  kmerhash kmer = new_kmerhash (dense); // false = 4 bits per site (o.w. 2 bits) 
  link_kmerhash_to_dna_sequence (kmer, dna, dna_length); // not ideal, since we can create ONE kmerhash for whole analysis 
  bbit_minhash_sketch oph = new_bbit_minhash_sketch (n_bits, kmer->n_hash);
  oph->kmer = kmer; 

  while ( kmerhash_iterator (oph->kmer)) for (i = 0; i < oph->n_sketch; i++) {
    prefix = oph->kmer->hash[i] >> oph->n_bits; // hash with (64 - n_bits) bits
    suffix = oph->kmer->hash[i] &  oph->suffix_mask;   // int value of least sig n_bits of hash
    if (oph->sketch[i][suffix] > prefix) oph->sketch[i][suffix] = prefix;
  }
  return oph;
}

new_cm_sketch (int max_vector_size)
{ // this may not be a locally-sensitive hashing (LSH) since similar inputs go to distinct buckets
  int i, j;
  cm_sketch cm = (cm_sketch) biomcmc_malloc (sizeof (struct cm_sketch_struct));
  if (max_vector_size < 16) max_vector_size = 16;
  cm->size = max_vector_size;
  cm->mod = 0xffffffff / (max_vector_size + 1); // plusone for case hash == MAX 
  cm->count = 0; 
  cm->freq = (int**) biomcmc_malloc (8 * sizeof (int*));
  for (i = 0; i < 8; i++) {
    cm->freq[i] = (int*) biomcmc_malloc (cm->size * sizeof (int));
    for (j = 0; j < cm->size; j++) cm->freq[i][j] = 0;
  }
// obs: when adding, can check if hash > max_size (only then hash/(max/nbuck)
  return cm;
}
  for (i=0; i < 8; i++) cm->freq[i][ (int) (h32[i]/cm->mod) ]++;  
  cm->count++;

void
compare_cm_sketches (cm_sketch cm1, cm_sketch cm2, double *result)
{
  int i, j;
  double x, frac = (double)(cm1->count) / (double)(cm2->count); 
  for (i=0; i<8; i++) result[i] = 0.;
  if (cm1->size != cm2->size) biomcmc_error ("can't compare sketches of distinct sizes");
  for (i=0; i<8; i++) {
    for (j = 0; j < cm1->size; j++) {
      // a/m - b/n = (na - mb)/mn = dividing both terms by n = (na/n - mb/n)/ (mn/n) = (a - m/n x b)/m
      x = (double)(cm1->freq[i][j]) - frac * (double)(cm2->freq[i][j]);
      result[i] += x * x;
    }
    result[i] /= (double)(cm1->count * cm1->count); 
  }
}

*/

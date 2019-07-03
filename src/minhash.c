/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "minhash.h"

#define xFULL64 0xffffffffffffffffUL // max possible uint64_t

minhash new_minhash (int sketch_size, int n_sketch);
onephash new_onephash (int n_bits, int n_sketch);


minhash
new_minhash (int sketch_size, int n_sketch)
{
  int i;
  minhash minh = (minhash) biomcmc_malloc (sizeof (struct minhash_struct));
  if (sketch_size < 8)    sketch_size = 8; // ideal 16~128 
  if (sketch_size > 2048) sketch_size = 2048; // should be small, since longer kmer have bigger sizes
  minh->n_sketch = n_sketch; // should match hashes available at kmerhash (or added by hand from those)
  minh->sketch = (heap64*) biomcmc_malloc (minh->n_sketch * sizeof (heap64*));
  for (i = 0; i < minh->n_sketch; i++) minh->sketch[i] = new_heap64 (sketch_size * (i+1)); // longer kmers have longer sketches 
  return minh;
}

void
del_minhash (minhash minh)
{
  int i;
  if (!minh) return;
  if (minh->sketch) {
    for (i=minh->n_sketch-1; i >=0; i--) del_heap64 (minh->sketch[i]);
    free (minh->sketch);
  }
  free (minh);
}

minhash
new_minhash_from_dna (char *dna, size_t dna_length, int sketch_size)
{
  int i;
  kmerhash kmer = new_kmerhash_from_dna_sequence (dna, dna_length, true); // false = 4 bits per site (o.w. 2 bits) 
  minhash minh = new_minhash (sketch_size, kmer->n_hash); // use all hashes available by kmerhash

  while ( kmerhash_iterator (kmer)) {
    for (i = 0; i < minh->n_sketch; i++) heap64_insert (minh->sketch[i], kmer->hash[i]);
  }

  for (i = 0; i < minh->n_sketch; i++) heap64_finalise_heap_qsort (minh->sketch[i]);
  del_kmerhash (kmer);
  return minh;
}

void
compare_minhash (minhash mh1, minhash mh2, double *result)
{
  int i, j1, j2, n1, n2, common, compared;
  uint64_t *h1, *h2;
  if (mh1->n_sketch != mh2->n_sketch) biomcmc_error ("can't map between minhash sketches");
  for (i = 0; i < mh1->n_sketch; i++) {
    if (mh1->sketch[i]->heap_size != mh2->sketch[i]->heap_size) biomcmc_error ("Distinct size sketches: comparable in theory but likely a mistake");
    h1 = mh1->sketch[i]->hash;    n1 = mh1->sketch[i]->n;
    h2 = mh2->sketch[i]->hash;    n2 = mh2->sketch[i]->n;
    common = compared = 0;
    for (j1 = j2 = 0; (j1 < n1) && (j2 < n2) && compared < mh1->sketch[i]->heap_size; compared++) { // both are in decreasing order
      if (h1[j1] > h2[j2]) j1++; 
      else if (h1[j1] < h2[j2]) j2++; 
      else { common++; j1++; j2++; } // same hash in both
    }
    if (compared < mh1->sketch[i]->heap_size) compared += (n1 - j1 + n2 - j2); // try to complete the union operation (following Mash)
    if (compared > mh1->sketch[i]->heap_size) compared = mh1->sketch_size;
    result[i] = (double) (common) / (double) (compared);
    // dist = -log(2.*result/(1.+result)) / kmersize;
    // pvalue =  gsl_cdf_binomial_Q(x - 1, r, sketchSize); r=pX*pY/(pX+ pY- pX*pY); pX = 1./(1.+kmerSpace/n2); pY = 1./(1.+kmerSpace/n1); kmerspace=4**kmersize
  }
}

onephash
new_onephash (int n_bits, int n_sketch)
{
  int i, j;
  onephash oph = (onephash) biomcmc_malloc (sizeof (struct onephash_struct));
  if (n_bits < 3)  n_bits = 3; // good value 8~10 
  if (n_bits > 14) n_bits = 14;
  oph->n_bits = n_bits;
  oph->suffix_mask = oph->precision_mask = 0UL;
  for (i = 0; i < oph->n_bits; i++) oph->suffix_mask |= (1UL << i); 
  for (i = 0; i < (64 - oph->n_bits - 2); i++) oph->precision_mask |= (1UL << i);

  oph->sketch_size = 1UL << oph->n_bits; // each vector size will be 2^n_bits

  oph->n_sketch = n_sketch; // should match hashes available at kmerhash (or added by hand from those)
  oph->sketch = (uint64_t**) biomcmc_malloc (oph->n_sketch * sizeof (uint64_t*));
  for (i = 0; i < oph->n_sketch; i++) {
    oph->sketch[i] = (uint64_t*) biomcmc_malloc (oph->sketch_size * sizeof (uint64_t));
    for (j = 0; j < oph->sketch_size; j++) oph->sketch[i][j] = xFULL64; // max possible uint64_t
  } 
  return oph;
}

void
del_onephash (onephash oph)
{
  int i;
  if (!oph) return;
  if (oph->sketch) {
    for (i = oph->n_sketch - 1; i >= 0; i--) if (oph->sketch[i]) free (oph->sketch[i]); 
    free (oph->sketch);
  }
  free (oph);
}

onephash
new_onephash_from_dna (char *dna, size_t dna_length, int n_bits)
{
  int i;
  uint64_t prefix = 0UL;
  uint32_t suffix = 0UL;

  kmerhash kmer = new_kmerhash_from_dna_sequence (dna, dna_length, true); // false = 4 bits per site (o.w. 2 bits) 
  onephash oph = new_onephash (n_bits, kmer->n_hash);

  while ( kmerhash_iterator (kmer)) for (i = 0; i < oph->n_sketch; i++) {
    prefix = kmer->hash[i] >> oph->n_bits; // hash with (64 - n_bits) bits
    suffix = kmer->hash[i] &  oph->suffix_mask;   // int value of least sig n_bits of hash
    if (oph->sketch[i][suffix] > prefix) oph->sketch[i][suffix] = prefix;
  }
  del_kmerhash (kmer);
  return oph;
}

void
compare_onephash (onephash oh1, onephash oh2, double *result)
{
  int i, j, compared, common; 
  if (oh1->sketch_size != oh2->sketch_size) biomcmc_error ("can't compare one-perm hashes of distinct sketch sizes");
  if (oh1->n_sketch != oh2->n_sketch) biomcmc_error ("can't compare one-perm hashes with different number of sketches");
  for (i = 0; i < oh1->n_sketch; i++) {
    compared = common = 0;
    for (j = 0; j < oh1->sketch_size; j++) 
      if ((oh1->sketch[i][j] != xFULL64) && (oh2->sketch[i][j] != xFULL64)) {
        compared++;
        if ((oh1->sketch[i][j] & oh1->precision_mask) == (oh2->sketch[i][j] & oh2->precision_mask)) common++;
    }
    result[i] = (double) (common) / (double) (compared);
  }
}

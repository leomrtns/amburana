/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "minhash.h"

/* almost identical to dna4bits[] from alignment[], but has forward and reverse */
uint8_t  dna4bits[256][2] = {{0xff}}; /* DNA base to bitpattern translation, with 1st element set to arbitrary value */

void update_minhash_from_fixedhash (minhash mh, uint64_t hash_f, uint64_t hash_r);
static void fixedhash_values_from_16mer (int dnachar, uint64_t *hf, uint64_t *hr);
void initialize_dna_to_bit_table (void);
void biomcmc_hashint64_to_vector (uint64_t x, uint32_t *out);
void update_cm_sketch_from_fixedhash (cm_sketch cm, uint64_t hash_f, uint64_t hash_r);

cm_sketch
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

  return cm;
}

void
del_cm_sketch (cm_sketch cm)
{
  int i;
  if (!cm) return;
  if (cm->freq) {
    for (i=7; i >=0; i--) if (cm->freq[i]) free (cm->freq[i]);
    free (cm->freq);
  }
  free (cm);
}

cm_sketch 
new_fixedhash_sketch_from_dna (char *dna, int dna_length, int sketch_size)
{
  int i;
  cm_sketch cm = new_cm_sketch (sketch_size);
  uint64_t hash_f = 0UL, hash_r = 0UL;
  if (dna4bits[0][0] == 0xff) initialize_dna_to_bit_table ();
  for (i = 0; i < 15; i++) fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r); 
  for (i = 15; i < dna_length; i++) {
    fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r);
    update_cm_sketch_from_fixedhash (cm, hash_f, hash_r);
    //for (j = 0; j < 16; j++) printf ("%d ", (int)((hash_f >> 4*j) & 15LL)); for (j = 0; j < 16; j++) printf ("%c", dna[i-j]);
  }
  return cm;
}

void
update_cm_sketch_from_fixedhash (cm_sketch cm, uint64_t hash_f, uint64_t hash_r)
{
  uint32_t h32[8]; // int = (-x,x); uint = (0,2x) 
  uint64_t small_h=0;  
  int i; 
  if (hash_f < hash_r) small_h = hash_f;
  else small_h = hash_r; 

  biomcmc_hashint64_to_vector (small_h, h32); // first four elements of h32[] are filled
  small_h = biomcmc_hashint64_salted (small_h, 4); // avalanche used in xxhash 
  biomcmc_hashint64_to_vector (small_h, h32 + 4); // last four elements of h32[] are filled
  for (i=0; i < 8; i++) cm->freq[i][ (int) (h32[i]/cm->mod) ]++;  
  cm->count++;
}

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

minhash
new_minhash (int sketch_size)
{
  int i;
  minhash mh = (minhash) biomcmc_malloc (sizeof (struct minhash_struct));
  if (sketch_size < 16) sketch_size = 16;
  mh->sketch_size = sketch_size;
  mh->sketch = (heap64*) biomcmc_malloc (4 * sizeof (heap64*));
  for (i = 0; i < 4; i++) mh->sketch[i] = new_heap64 (sketch_size); 
  return mh;
}

void
del_minhash (minhash mh)
{
  int i;
  if (!mh) return;
  if (mh->sketch) {
    for (i=3; i >=0; i--) del_heap64 (mh->sketch[i]);
    free (mh->sketch);
  }
  free (mh);
}

minhash
new_minhash_from_dna (char *dna, int dna_length, int sketch_size)
{
  int i;
  minhash mh = new_minhash (sketch_size);
  uint64_t hash_f = 0UL, hash_r = 0UL;
  if (dna4bits[0][0] == 0xff) initialize_dna_to_bit_table ();
  for (i = 0; i < 15; i++) fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r); 
  for (i = 15; i < dna_length; i++) {
    fixedhash_values_from_16mer ((int) dna[i], &hash_f, &hash_r);
    update_minhash_from_fixedhash (mh, hash_f, hash_r);
    //for (j = 0; j < 16; j++) printf ("%d ", (int)((hash_f >> 4*j) & 15LL)); for (j = 0; j < 16; j++) printf ("%c", dna[i-j]);
  }
  for (i=0;i<4;i++) heap64_finalise_heap_qsort (mh->sketch[i]);
  return mh;
}

void
update_minhash_from_fixedhash (minhash mh, uint64_t hash_f, uint64_t hash_r)
{
  uint64_t small_h = 0UL, out[2];  
  if (hash_f < hash_r) small_h = hash_f;
  else small_h = hash_r;

  biomcmc_murmurhash3_128bits (&small_h, 8, 57, out); 
  heap64_insert (mh->sketch[0], out[0]);
  heap64_insert (mh->sketch[1], out[1]);
  heap64_insert (mh->sketch[2], biomcmc_xxh64 (&small_h, 4, 3)); // 4 means only first 32 bits
  small_h = small_h >> 32UL; // last 32 bits
  heap64_insert (mh->sketch[3], biomcmc_xxh64 (&small_h, 4, 3));  // 4 means only first 32 bits
}

void
compare_minhashes (minhash mh1, minhash mh2, double *result)
{ // obs: minhash calcs distance (not similarity), counts how many k-mers are in common until intersection is sketch_size
  int i, j1, j2, n1, n2, common, compared;
  uint64_t *h1, *h2;
  if (mh1->sketch_size != mh2->sketch_size) biomcmc_error ("can't compare sketches of distinct sizes");
  for (i=0; i<4; i++) { 
    h1 = mh1->sketch[i]->hash;    n1 = mh1->sketch[i]->n;
    h2 = mh2->sketch[i]->hash;    n2 = mh2->sketch[i]->n;
    common = compared = 0;
    for (j1 = j2 = 0; (j1 < n1) && (j2 < n2) && compared < mh1->sketch_size; compared++) { // both are in decreasing order
      if (h1[j1] > h2[j2]) j1++; 
      else if (h1[j1] < h2[j2]) j2++; 
      else { common++; j1++; j2++; } // same hash in both
    }
    if (compared < mh1->sketch_size) compared += (n1 - j1 + n2 - j2); // try to complete the union operation (following Mash)
    if (compared > mh1->sketch_size) compared = mh1->sketch_size;
    result[i] = (double) (common) / (double) (compared);
    // dist = -log(2.*result/(1.+result)) / kmersize;
    // pvalue =  gsl_cdf_binomial_Q(x - 1, r, sketchSize); r=pX*pY/(pX+ pY- pX*pY); pX = 1./(1.+kmerSpace/n2); pY = 1./(1.+kmerSpace/n1); kmerspace=4**kmersize
  }
}

static void
fixedhash_values_from_16mer (int dnachar, uint64_t *hf, uint64_t *hr)
{
  *hf = *hf << 4 | dna4bits[dnachar][0]; // forward
  *hr = *hf >> 4 | ((uint64_t)(dna4bits[dnachar][1]) << 60UL); // reverse
}

void
initialize_dna_to_bit_table (void)
{
  int i;
  for (i = 0; i < 256; i++) dna4bits[i][0] = dna4bits[i][1] = 0;
  /* The ACGT is PAUP convention (and maybe DNAml, fastDNAml); PAML uses TCAG ordering */
  dna4bits['A'][0] = 1;   dna4bits['A'][1] = 8;  /* .   A */ /* 0001 */ /* reverse is 'T'    = 8  */
  dna4bits['B'][0] = 14;  dna4bits['B'][1] = 7;  /* .TGC  */ /* 1110 */ /* reverse is 'ACG'  = 7  */
  dna4bits['C'][0] = 2;   dna4bits['C'][1] = 4;  /* .  C  */ /* 0010 */ /* reverse is 'G'    = 4  */
  dna4bits['D'][0] = 13;  dna4bits['D'][1] = 11; /* .TG A */ /* 1101 */ /* reverse is 'TCA'  = 11 */
  dna4bits['G'][0] = 4;   dna4bits['G'][1] = 2;  /* . G   */ /* 0100 */ /* reverse is 'C'    = 2  */
  dna4bits['H'][0] = 11;  dna4bits['H'][1] = 13; /* .T CA */ /* 1011 */ /* reverse is 'TGA'  = 13 */
  dna4bits['K'][0] = 12;  dna4bits['K'][1] = 3;  /* .TG   */ /* 1100 */ /* reverse is 'AC'   = 3  */
  dna4bits['M'][0] = 3;   dna4bits['M'][1] = 12; /* .  CA */ /* 0011 */ /* reverse is 'TG'   = 12 */
  dna4bits['N'][0] = 15;  dna4bits['N'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['O'][0] = 15;  dna4bits['O'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['R'][0] = 5;   dna4bits['R'][1] = 10; /* . G A */ /* 0101 */ /* reverse is 'TC'   = 10 */
  dna4bits['S'][0] = 6;   dna4bits['S'][1] = 6;  /* . GC  */ /* 0110 */ /* reverse is 'GC'   = 6  */
  dna4bits['T'][0] = 8;   dna4bits['T'][1] = 1;  /* .T    */ /* 1000 */ /* reverse is 'A'    = 1  */
  dna4bits['U'][0] = 8;   dna4bits['U'][1] = 1;  /* .T    */ /* 1000 */ /* reverse is 'A'    = 1  */
  dna4bits['V'][0] = 7;   dna4bits['V'][1] = 14; /* . GCA */ /* 0111 */ /* reverse is 'TGC'  = 14 */
  dna4bits['W'][0] = 9;   dna4bits['W'][1] = 9;  /* .T  A */ /* 1001 */ /* reverse is 'TA'   = 9  */
  dna4bits['X'][0] = 15;  dna4bits['X'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['Y'][0] = 10;  dna4bits['Y'][1] = 5;  /* .T C  */ /* 1010 */ /* reverse is 'GA'   =  5 */
  dna4bits['?'][0] = 15;  dna4bits['?'][1] = 15; /* .TGCA */ /* 1111 */ /* reverse is 'TGCA' = 15 */
  dna4bits['-'][0] = 0;   dna4bits['-'][1] = 0;  /* .TGCA */ /* fifth state */
}
// trick : 0xAAAAAAAAAA 0x5555555555 (A=1010 5=0101) for 32bits

void
biomcmc_hashint64_to_vector (uint64_t x, uint32_t *out) /* 64 bits, splits into two 32bits blocks */
{
  uint32_t low = x;
  uint32_t high = x >> 32UL;
  out[0] =biomcmc_hashint_salted (low, 1);
  out[1] =biomcmc_hashint_salted (high, 2);
  x =biomcmc_hashint64_salted (x, 5);
  out[2] = (uint32_t) x;
  out[3] = (uint32_t) (x >> 32UL);
}



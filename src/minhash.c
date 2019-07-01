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
void fixedhash_values_from_16mer (int dnachar, uint64_t *hf, uint64_t *hr);
void initialize_dna_to_bit_table (void);

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
  uint64_t small_h = 0UL;  
  if (hash_f < hash_r) small_h = hash_f;
  else small_h = hash_r;

  heap64_insert (mh->sketch[0], biomcmc_xxh64 (&small_h, 8, 3));
  heap64_insert (mh->sketch[1], biomcmc_murmurhash3 (&small_h, 8, 57, NULL)); // NULL since I dont need out[]
  heap64_insert (mh->sketch[2], biomcmc_hashint64_salted (small_h, 5));
  heap64_insert (mh->sketch[3], biomcmc_hashint64_salted (small_h, 10));
}

void
compare_minhashes (minhash mh1, minhash mh2, double *result)
{
  int i, j1, j2, n1, n2, count;
  uint64_t *h1, *h2;
  if (mh1->sketch_size != mh2->sketch_size) biomcmc_error ("can't compare sketches of distinct sizes");
  for (i=0; i<4; i++) { 
    count = 0;
    h1 = mh1->sketch[i]->hash;    n1 = mh1->sketch[i]->n;
    h2 = mh2->sketch[i]->hash;    n2 = mh2->sketch[i]->n;

    for (j1 = j2 = 0; (j1 < n1) && (j2 < n2); ) { // both are in decreasing order
      if (h1[j1] > h2[j2]) j1++; 
      else if (h1[j1] < h2[j2]) j2++; 
      else { count++; j1++; j2++; } // same hash in both
    }
    result[i] = (double) (count) / (double) (n1 > n2 ? n1 : n2);
  }
}

void
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


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
  double scale = 1.;
  minhash minh = (minhash) biomcmc_malloc (sizeof (struct minhash_struct));
  if (sketch_size < 8)    sketch_size = 8; // ideal 16~128 
  if (sketch_size > 2048) sketch_size = 2048; // should be small, since longer kmer have bigger sizes
  minh->n_sketch = n_sketch; // should match hashes available at kmerhash (or added by hand from those)
  minh->sketch = (heap64*) biomcmc_malloc (minh->n_sketch * sizeof (heap64*));
  for (i = 0; i < minh->n_sketch; i++) {
    minh->sketch[i] = new_heap64 ((int)(scale * (double)(sketch_size))); // longer kmers have longer sketches 
    scale *= 1.2; // doubles every ~4 iterations
  }
  minh->kmer = NULL;
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
  del_kmerhash (minh->kmer);
  free (minh);
}

minhash
new_minhash_from_dna (char *dna, size_t dna_length, int sketch_size, bool dense)
{
  int i;
  kmerhash kmer = new_kmerhash_from_dna_sequence (dna, dna_length, dense); // false = 4 bits per site (o.w. 2 bits) 
  minhash minh = new_minhash (sketch_size, kmer->n_hash); // use all hashes available by kmerhash
  minh->kmer = kmer;

  while ( kmerhash_iterator (minh->kmer)) {
    for (i = 0; i < minh->n_sketch; i++) heap64_insert (minh->sketch[i], minh->kmer->hash[i]);
  }

  for (i = 0; i < minh->n_sketch; i++) heap64_finalise_heap_qsort (minh->sketch[i]);
  return minh;
}

void
compare_minhash (minhash mh1, minhash mh2, double *distance)
{
  int i, j1, j2, n1, n2, common, compared;
  uint64_t *h1, *h2;
  double r;
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
    if (common == compared) distance[i] = 0.;
    else if (common == 0) distance[i] = 1.;
    else {
      r = (double) (common) / (double) (compared);
      distance[i] = -log (2. * r/(1. + r)) / (double) (mh1->kmer->nsites_kmer[i]);
    }
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
  oph->kmer = NULL;
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
  del_kmerhash (oph->kmer);
  free (oph);
}

onephash
new_onephash_from_dna (char *dna, size_t dna_length, int n_bits, bool dense)
{
  int i;
  uint64_t prefix = 0UL;
  uint32_t suffix = 0UL;

  kmerhash kmer = new_kmerhash_from_dna_sequence (dna, dna_length, dense); // false = 4 bits per site (o.w. 2 bits) 
  onephash oph = new_onephash (n_bits, kmer->n_hash);
  oph->kmer = kmer; 

  while ( kmerhash_iterator (oph->kmer)) for (i = 0; i < oph->n_sketch; i++) {
    prefix = oph->kmer->hash[i] >> oph->n_bits; // hash with (64 - n_bits) bits
    suffix = oph->kmer->hash[i] &  oph->suffix_mask;   // int value of least sig n_bits of hash
    if (oph->sketch[i][suffix] > prefix) oph->sketch[i][suffix] = prefix;
  }
  return oph;
}

void
compare_onephash (onephash oh1, onephash oh2, double *distance)
{
  int i, j, compared, common; 
  double r;
  if (oh1->sketch_size != oh2->sketch_size) biomcmc_error ("can't compare one-perm hashes of distinct sketch sizes");
  if (oh1->n_sketch != oh2->n_sketch) biomcmc_error ("can't compare one-perm hashes with different number of sketches");
  for (i = 0; i < oh1->n_sketch; i++) {
    compared = common = 0;
    for (j = 0; j < oh1->sketch_size; j++) 
      if ((oh1->sketch[i][j] != xFULL64) && (oh2->sketch[i][j] != xFULL64)) {
        compared++;
        if ((oh1->sketch[i][j] & oh1->precision_mask) == (oh2->sketch[i][j] & oh2->precision_mask)) common++;
    }
    if (common == compared) distance[i] = 0.;
    else if (common == 0) distance[i] = 1.;
    else {
      r = (double) (common) / (double) (compared);
      distance[i] = -log (2. * r/(1. + r)) / (double) (oh1->kmer->nsites_kmer[i]);
    }
  }
}
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
#include <amburana.h>
#include "simple_heap.c"
#include "minhash.c"
//#include "parasail/matrices/nuc44.h"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_file *fasta;
  struct arg_int  *sketch;
  struct arg_int  *nbits;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .fasta   = arg_filen(NULL, NULL, NULL, 1, 1, "fasta file"),
    .sketch  = arg_int0("s", "size", "<n>", "sketch size for minhash"),
    .nbits   = arg_int0("b", "bits", "<n>", "number of bits for suffix of one-permutation minhash"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.fasta, params.sketch, params.nbits, params.end};
  params.argtable = argtable; 
  params.sketch->ival[0] = 64; // default (must be before parsing)
  params.nbits->ival[0] = 5; // default (must be before parsing)
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  arg_parse (argc, argv, params.argtable); // returns >0 if errors were found, but this info also on params.end->count
  print_usage (params, argv[0]);
  return params;
}

void
del_arg_parameters (arg_parameters params)
{
  if (params.help) free (params.help);
  if (params.version) free (params.version);
  if (params.fasta) free (params.fasta);
  if (params.sketch) free (params.sketch);
  if (params.nbits) free (params.nbits);
  if (params.end) free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return;

  if (params.end->count) {  // params.end holds error messages
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("Error when reading arguments from command line\n\n");
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("Just testing at the moment. Move along, nothing to see here!\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf (" And here is where a nice explanation would be.\n\n");
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i,j,k;
  clock_t time0, time1;
  sketch_set *sset;
  //parasail_result_t* nwresult;
  alignment aln;

  arg_parameters params = get_parameters_from_argv (argc, argv);

  time0 = clock (); 
  aln = read_alignment_from_file ((char*) params.fasta->filename[0]);
  sset = (sketch_set*) biomcmc_malloc (aln->ntax * sizeof (sketch_set));
  time1 = clock (); fprintf (stderr, "  time to read alignment and set sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=0; i < aln->ntax; i++) sset[i] = new_onephash_from_dna (aln->character->string[i], aln->character->nchars[i], params.nbits->ival[0], true);
  time1 = clock (); fprintf (stderr, "  time to calculate sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) {
    printf ("\n%38s %38s ", aln->taxlabel->string[j], aln->taxlabel->string[i]);

    //nwresult = parasail_nw_banded (aln->character->string[i], aln->character->nchars[i], 
    //                               aln->character->string[j], aln->character->nchars[j], 1, 4, 16, &parasail_nuc44); 
    //printf (" %14d    ", parasail_result_get_score (nwresult));

    compare_onephash (oh[i], oh[j], dist);
    compare_minhash (mh[i], mh[j], dist+7);
    for (k=0;k<14;k++) printf ("%9.7lf ", dist[k]);
    //parasail_result_free(nwresult);
  }

  for (i= aln->ntax -1; i >- 0; i--) { del_onephash (oh[i]); del_minhash (mh[i]); }
  if (oh) free (oh);
  if (mh) free (mh);
  del_alignment (aln);

  del_arg_parameters (params);
  return EXIT_SUCCESS;
}


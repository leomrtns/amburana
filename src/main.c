#include <amburana.h>
#include "simple_heap.c"
#include "minhash.c"
#include "parasail/matrices/nuc44.h"

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
  onephash *oh;
  minhash *mh;
  double dist[8];
  //parasail_result_t* nwresult;
  alignment aln;

  arg_parameters params = get_parameters_from_argv (argc, argv);

  time0 = clock (); 
  aln = read_alignment_from_file ((char*) params.fasta->filename[0]);
  oh = (onephash*) biomcmc_malloc (aln->ntax * sizeof (onephash));
  mh = (minhash*)  biomcmc_malloc (aln->ntax * sizeof (minhash));
  time1 = clock (); printf ("  time to read alignment and set sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=0; i < aln->ntax; i++) oh[i] = new_onephash_from_dna (aln->character->string[i], aln->character->nchars[i], params.nbits->ival[0]);
  time1 = clock (); printf ("  time to calculate sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=0; i < aln->ntax; i++) mh[i] = new_minhash_from_dna (aln->character->string[i], aln->character->nchars[i], params.sketch->ival[0]);
  time1 = clock (); printf ("  time to calculate minhashes: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) {
    printf ("\n%38s %38s ", aln->taxlabel->string[j], aln->taxlabel->string[i]);

    //nwresult = parasail_nw_banded (aln->character->string[i], aln->character->nchars[i], 
    //                               aln->character->string[j], aln->character->nchars[j], 1, 4, 16, &parasail_nuc44); 
    //printf (" %14d    ", parasail_result_get_score (nwresult));

    compare_onephash (oh[i], oh[j], dist);
    compare_minhash (mh[i], mh[j], dist+4);
    for (k=0;k<8;k++) printf ("%9.7lf ", dist[k]);
    //parasail_result_free(nwresult);
  }

  for (i= aln->ntax -1; i >- 0; i--) { del_onephash (oh[i]); del_minhash (mh[i]); }
  if (oh) free (oh);
  if (mh) free (mh);
  del_alignment (aln);

  del_arg_parameters (params);
  return EXIT_SUCCESS;
}


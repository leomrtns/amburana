#include "amburana.h" 
#include "simple_heap.c"
#include "minhash.c"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_file *spname;
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
    .help = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .spname = arg_file0("s","species", "<file name>", "file "),
    .end  = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.spname, params.end};
  params.argtable = argtable; 
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
  if (params.spname) free (params.spname);
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
  cm_sketch *cm;
  minhash *mh;
  double dist[8];
  alignment aln;

  arg_parameters params = get_parameters_from_argv (argc, argv);

  time0 = clock (); 
  aln = read_alignment_from_file ((char*) params.spname->filename[0]);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  cm = (cm_sketch*) biomcmc_malloc (aln->ntax * sizeof (cm_sketch));
  mh = (minhash*) biomcmc_malloc (aln->ntax * sizeof (minhash));

  for (i=0; i < aln->ntax; i++) cm[i] = new_fixedhash_sketch_from_dna (aln->character->string[i], aln->character->nchars[i], 128);
  time1 = clock (); printf ("  time to calculate sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i=0; i < aln->ntax; i++) mh[i] = new_minhash_from_dna (aln->character->string[i], aln->character->nchars[i], 128);
  time1 = clock (); printf ("  time to calculate minhashes: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) {
    printf ("\n%38s %38s ", aln->taxlabel->string[j], aln->taxlabel->string[i]); 
    compare_cm_sketches (cm[i], cm[j], dist);
    for (k=0;k<8;k++) printf ("%8.6lf ", dist[k]);
    compare_minhashes (mh[i], mh[j], dist);
    for (k=0;k<4;k++) printf ("%8.6lf ", dist[k]);
  }

  for (i= aln->ntax -1; i >- 0; i--) { del_cm_sketch (cm[i]); del_minhash (mh[i]); }
  if (cm) free (cm);
  if (mh) free (mh);
  del_alignment (aln);

  del_arg_parameters (params);
  return EXIT_SUCCESS;
}


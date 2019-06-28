#include "simple_heap.h" 

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
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
    .end  = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.end};
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
  int i,j, nreps = 400, heapsize=500000;
  uint64_t rad;
  double elapsed;
  clock_t time0, time1;
  heap64 h64;

  biomcmc_random_number_init (0);
  arg_parameters params = get_parameters_from_argv (argc, argv);

  time0 = clock (); elapsed = 0.;
  for (i=0; i < nreps; i++) {
    h64 = new_heap64 (heapsize);
    for (j=0; j < heapsize+1; j++) {
      rad = (uint64_t) biomcmc_rng_get_32(); 
      heap64_insert (h64, rad);
    }
    heap64_finalise_heap_pop (h64);
    time1 = clock (); 
    if (!(i%50)) printf ("heap: %.8lf secs\n", (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); 
    elapsed += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;
    del_heap64(h64);
  }
  printf ("average heap : %lf secs\n", elapsed/(double)(nreps));


  time0 = clock (); elapsed = 0.;
  for (i=0; i < nreps; i++) {
    h64 = new_heap64 (heapsize);
    for (j=0; j < heapsize+1; j++) {
      rad = (uint64_t) biomcmc_rng_get_32(); 
      heap64_insert (h64, rad);
    }
    heap64_finalise_heap_qsort (h64);
    time1 = clock (); 
    if (!(i%50)) printf ("qsort: %.8lf secs\n", (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); 
    elapsed += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;
    del_heap64(h64);
  }
  printf ("average heap : %lf secs\n", elapsed/(double)(nreps));


  biomcmc_random_number_finalize ();
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}


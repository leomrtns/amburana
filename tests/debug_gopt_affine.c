#include <amburana.h>
#include <zlib.h>

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_file *fasta;
  struct arg_int  *sketch;
  struct arg_int  *nbits;
  struct arg_int  *kmerset;
  struct arg_int  *minsamp;
  struct arg_dbl  *epsilon;
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
    .kmerset = arg_int0("k", "kmer", "{0,1,2,3}", "code for set of kmers"),
    .minsamp = arg_int0("m", "min", "{2,3,...}", "min samples to be considered core (OPTICS)"),
    .epsilon = arg_dbl0("e", "epsilon", "<float>", "max dist to be considered neighbour (OPTICS)"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.fasta, params.sketch, params.nbits, params.kmerset, params.minsamp, params.epsilon, params.end};
  params.argtable = argtable; 
  params.sketch->ival[0] = 256; // default (must be before parsing)
  params.nbits->ival[0] = 10; // default (must be before parsing)
  params.kmerset->ival[0] = 3;
  params.minsamp->ival[0] = 2;
  params.epsilon->dval[0] = 1.;
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
  if (params.kmerset) free (params.kmerset);
  if (params.minsamp) free (params.minsamp);
  if (params.epsilon) free (params.epsilon);
  if (params.end) free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  int i;
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
    printf ("The choices for the sets of kmers are:\n");
    for (i=0;i<4;i++) printf ("%d\t for %s analysis\n", i, biomcmc_kmer_class_string[i]);
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j, k;
  sketch_distance_gen sd;
  distance_generator dg;
  kmerhash kmer;
  alignment aln;
  goptics_cluster gop;
  affineprop_cluster ap;
  clock_t time0, time1;
  double timing1 = 0., timing2 = 0.;
  char *iscore[]={"     ","core "};

  arg_parameters params = get_parameters_from_argv (argc, argv);
  kmer = new_kmerhash (params.kmerset->ival[0]); 

  aln = read_alignment_from_file ((char*) params.fasta->filename[0]);
  sd = new_sketch_distance_gen (aln, kmer, params.sketch->ival[0], params.nbits->ival[0]);
  dg = new_distance_generator_from_sketch_set (sd);
  ap = new_affineprop_cluster (dg);

  for (k=0; k < dg->n_distances; k++) {
    distance_generator_set_which_distance (dg, k);
    time0 = clock ();

    gop = new_goptics_cluster_run (dg, params.minsamp->ival[0], params.epsilon->dval[0]);
    assign_goptics_clusters (gop, 0.1 * params.epsilon->dval[0]);
    time1 = clock (); timing1 += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;

    for (i=-1; i < gop->n_clusters; i++) {
      printf ("\t\t\t\t k=%d  GOP  cluster %d:\n", k, i);
      for (j=0; j < gop->d->n_samples; j++) if (gop->cluster[j] == i) 
        printf (" %s %-40s \t %7.5lf    %7.5lf\n",iscore[(int)gop->core[j]], aln->taxlabel->string[j], gop->reach_distance[j], gop->core_distance[j]);
    }
    del_goptics_cluster (gop); // recreated every time

    affineprop_run (ap, 200, 0.6, 0.5); // n_iter, preference (in quantile), damping (weight of previous iteration)
    time1 = clock (); timing2 += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;

    for (i=0; i < ap->n_clusters; i++) {
      printf ("\t\t\t\t k=%2d Affine  cluster %2d: convergence: %3d exemplar :: %s\n", k, i, ap->n_converging, aln->taxlabel->string[ ap->exemplars[i] ]);
      for (j=0; j < ap->d->n_samples; j++) if (ap->cluster[j] == i) 
        printf ("%3d %-40s \n", j, aln->taxlabel->string[j]);
    }
  }
  fprintf (stderr, "calculated sketches in %lf secs and pairwise distances in %lf secs\n", sd->secs_sketches, sd->secs_distances);
  fprintf (stderr, "average time for gOPTICS: %lf and for affine prop %lf \n", timing1/(double)(dg->n_distances), timing2/(double)(dg->n_distances));
  fprintf (stderr, "timing for gOPTICS  %lf and affine prop %lf as given internally\n",  gop->timing_secs/(double)(dg->n_distances), ap->timing_secs/(double)(dg->n_distances));


  del_affineprop_cluster (ap);
  del_alignment (aln);
  del_sketch_distance_gen (sd);
  del_distance_generator (dg);
  del_kmerhash (kmer);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}


#include <amburana.h>
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

#ifndef TEST_FILE_DIR
#define TEST_FILE_DIR "./files/"
#endif
// use it like memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
char filename[2048] = TEST_FILE_DIR; // now we can memcpy() file names _after_ prefix_size
size_t prefix_size = strlen(TEST_FILE_DIR); // all modifications to filename[] come after prefix_size

START_TEST(heap_timing_function)
{
  int i,j, nreps = 10, heapsize=1000;
  uint64_t rad;
  double t_pop, t_qsort;
  clock_t time0, time1;
  heap64 h64;

  biomcmc_random_number_init (0);

  time0 = clock (); t_pop = 0.;
  for (i=0; i < nreps; i++) {
    h64 = new_heap64 (heapsize);
    for (j=0; j < heapsize+1; j++) {
      rad = (uint64_t) biomcmc_rng_get_32(); 
      heap64_insert (h64, rad);
    }
    heap64_finalise_heap_pop (h64);
    time1 = clock (); 
    t_pop+= (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;
    del_heap64(h64);
  }
  printf ("average heap : %lf secs\n", t_pop/(double)(nreps));

  time0 = clock (); t_qsort = 0.;
  for (i=0; i < nreps; i++) {
    h64 = new_heap64 (heapsize);
    for (j=0; j < heapsize+1; j++) {
      rad = (uint64_t) biomcmc_rng_get_32(); 
      heap64_insert (h64, rad);
    }
    heap64_finalise_heap_qsort (h64);
    time1 = clock (); 
    t_qsort+= (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1;
    del_heap64(h64);
  }
  printf ("average heap : %lf secs\n", t_qsort/(double)(nreps));
  biomcmc_random_number_finalize ();
  if (t_qsort > 2 * t_pop) ck_abort_msg ("qsort implementation too slow");
}
END_TEST

START_TEST(fixedhash_sketch_alignment_function)
{
  int i,j;
  cm_sketch *cm;
  double dist[8];
  clock_t time0, time1;
  alignment aln;

  time0 = clock ();
  memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
  aln = read_alignment_from_file (filename);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  cm = (cm_sketch*) biomcmc_malloc (aln->ntax * sizeof (cm_sketch));
  for (i=0; i < aln->ntax; i++) cm[i] = new_fixedhash_sketch_from_dna (aln->character->string[i], aln->character->nchars[i], 64);
  time1 = clock (); printf ("  time to calculate sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) compare_cm_sketches (cm[i], cm[j], dist);
  time1 = clock (); printf ("  time to compare sketches: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i= aln->ntax -1; i >- 0; i--) del_cm_sketch (cm[i]); 
  if (cm) free (cm);
  del_alignment (aln);
  if (false) ck_abort_msg ("dummy");
}
END_TEST

START_TEST(minhash_alignment_function)
{
  int i,j;
  minhash *mh;
  double dist[4];
  clock_t time0, time1;
  alignment aln;

  time0 = clock ();
  memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
  aln = read_alignment_from_file (filename);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  mh = (minhash*) biomcmc_malloc (aln->ntax * sizeof (minhash));
  for (i=0; i < aln->ntax; i++) mh[i] = new_minhash_from_dna (aln->character->string[i], aln->character->nchars[i], 64);
  time1 = clock (); printf ("  time to calculate minhashes: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i=1; i < aln->ntax; i++) for (j=0; j < i; j++) compare_minhashes (mh[i], mh[j], dist);
  time1 = clock (); printf ("  time to compare minhashes: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  for (i= aln->ntax -1; i >- 0; i--) del_minhash (mh[i]); 
  if (mh) free (mh);
  del_alignment (aln);
  if (false) ck_abort_msg ("dummy");
}
END_TEST


Suite * simpleheap_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("Simple MaxHeap");
  tc_case = tcase_create("max heap finalising timings");
  tcase_add_test(tc_case, heap_timing_function);
  tcase_add_test(tc_case, fixedhash_sketch_alignment_function);
  tcase_add_test(tc_case, minhash_alignment_function);
  suite_add_tcase(s, tc_case);
  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  sr = srunner_create (simpleheap_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed > 0) ? TEST_FAILURE:TEST_SUCCESS;
}

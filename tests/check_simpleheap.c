#include <simple_heap.h>
#include <minhash.h>
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

Suite * simpleheap_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("Simple MaxHeap");
  tc_case = tcase_create("max heap finalising timings");
  tcase_add_test(tc_case, heap_timing_function);
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

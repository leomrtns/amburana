#include "simple_heap.h" 

/* max heap structure for storing smallest hash values 
 * code inspired by https://gist.github.com/vgoel30/5d81e6abf9464930c1e126dab04d5be3  */

uint64_t heap64_get_maximum (heap64 h64) { return h64->hash[1]; }

heap64
new_heap64 (int heap_size)
{
  int i;
  heap64 h64 = (heap64) biomcmc_malloc (sizeof (struct heap64_struct)); 
  h64->n = 0;
  h64->heap_size = (heap_size < 4 ? 4:heap_size);
  h64->hash = (uint64_t*) biomcmc_malloc ((h64->heap_size + 1) * sizeof (uint64_t));
  for(i = 0; i < h64->heap_size + 1; i++) h64->hash[i] = 0UL;
  return h64;
}

void
del_heap64 (heap64 h64)
{
  if (h64) {
    if (h64->hash) free (h64->hash);
    free (h64);
  }
}

uint64_t 
heap64_remove_maximum (heap64 h64) // a.k.a. pop() 
{
	uint64_t max = h64->hash[1]; //the root is the maximum element
  if (h64->n == 0) return 0UL;
  h64->hash[1] = h64->hash[h64->n]; // replace the element at the top with last element in the heap
  h64->n -= 1; //reduce the total elements in the heap
  heap64_bubble_down (h64, 1);
  return max;
}

void 
heap64_insert (heap64 h64, uint64_t h) 
{
  if (h64->n == h64->heap_size) { // replace if smaller than maximum
    if (h < h64->hash[1]) { h64->hash[1] = h; heap64_bubble_down (h64, 1); }
  }
  else { // regular insert (a.k.a. push)
    h64->n += 1;
    h64->hash[h64->n] = h;
		heap64_bubble_up (h64, h64->n);
	}
}

void
heap64_bubble_down (heap64 h64, int p)
{//bubbles down an element into it's proper place in the heap
	int i, c = 2*p, max_index = p; //c= index of younger child, max_index = index of maximum child
  for(i = 0; i < 2; i++) if(c + i <= h64->n) {
    //check to see if the data at min_index is smaller than the data at the child
    if(h64->hash[max_index] < h64->hash[c+i])	max_index = c + i;
  }
	if (max_index != p) {
    uint64_t tmp = h64->hash[p];
    h64->hash[p] = h64->hash[max_index];
    h64->hash[max_index] = tmp;
		heap64_bubble_down (h64, max_index);
	}
}

void 
heap64_bubble_up (heap64 h64, int index)
{ //bubbles up the last element of the heap to maintain heap structure
  int parent = (index == 1 ? -1 : (int)(index/2));
  if (parent == -1) return;	//if we are at the root of the heap, no parent
  //if the parent node has a smaller data value, we need to bubble up
  if(h64->hash[parent] < h64->hash[index]) {
    uint64_t tmp = h64->hash[parent];
    h64->hash[parent] = h64->hash[index];
    h64->hash[index] = tmp;
    heap64_bubble_up (h64, parent);
  }
}

void
heap64_finalise_heap_pop (heap64 h64) // MUCH SLOWER THAN QSORT
{ // notice that this will be in _decreasing_ order
  int i, new_size = h64->n; // we need new_size since remove_max() will decrease h64->n
  uint64_t *newhash;

  newhash = (uint64_t*) biomcmc_malloc (new_size * sizeof (uint64_t));
  for (i=0; i < new_size; i++) newhash[i] = heap64_remove_maximum (h64);
  if (h64->hash) free (h64->hash);
  h64->hash = newhash;
  h64->heap_size = h64->n = new_size; // in case n < heap_size
}

void
heap64_finalise_heap_qsort (heap64 h64)
{
  h64->hash[0] = h64->hash[h64->n]; // heap goes from [1...n] and we want [0...n-1]
  qsort (h64->hash, h64->n, sizeof (uint64_t), compare_uint64_decreasing); // to be consistent with pop()
  h64->hash = (uint64_t*) biomcmc_realloc ((uint64_t*) h64->hash, h64->n * sizeof (uint64_t));
  h64->heap_size = h64->n;
}

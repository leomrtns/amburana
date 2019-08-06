#include "simple_heap.h" 

/* max heap structure for storing smallest hash values 
 * code inspired by https://gist.github.com/vgoel30/5d81e6abf9464930c1e126dab04d5be3  
 * and also for creating kNN of points for clustering (max heap as hash queue storing distances) */

static void heap64_bubble_down (heap64 h64, int p);
static void heap64_bubble_up (heap64 h64, int index);
static bool heap_hash64_item_already_here (heap_hash64 pq, hpq_item item, int p);
static void heap_hash64_bubble_down (heap_hash64 pq, int p);
static void heap_hash64_bubble_up (heap_hash64 pq, int index);
static int compare_hpq_item_increasing (const void *a, const void *b);


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
  if (h64) { if (h64->hash) free (h64->hash); free (h64); }
}

uint64_t heap64_get_maximum (heap64 h64) { return h64->hash[1]; }

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

static void
heap64_bubble_down (heap64 h64, int p)
{//bubbles down an element into it's proper place in the heap
	int i, c = 2*p, max_index = p; //c= index of younger child, max_index = index of maximum child
  for (i = 0; i < 2; i++) if(c + i <= h64->n) {
    //check to see if the data at min_index is smaller than the data at the child
    if (h64->hash[max_index] < h64->hash[c+i])	max_index = c + i;
  }
	if (max_index != p) {
    uint64_t tmp = h64->hash[p];
    h64->hash[p] = h64->hash[max_index];
    h64->hash[max_index] = tmp;
		heap64_bubble_down (h64, max_index);
	}
}

static void 
heap64_bubble_up (heap64 h64, int index)
{ //bubbles up the last element of the heap to maintain heap structure
  int parent = (index == 1 ? -1 : (int)(index/2));
  if (parent == -1) return;	//if we are at the root of the heap, no parent
  //if the parent node has a smaller data value, we need to bubble up
  if (h64->hash[parent] < h64->hash[index]) {
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
  qsort (h64->hash, h64->n, sizeof (uint64_t), compare_uint64_increasing); // to be consistent with pop() we should use "decreasing"
  h64->hash = (uint64_t*) biomcmc_realloc ((uint64_t*) h64->hash, h64->n * sizeof (uint64_t));
  h64->heap_size = h64->n;
}

heap_hash64 
new_heap_hash64 (int heap_size)
{
  int i;
  heap_hash64 pq = (heap_hash64) biomcmc_malloc (sizeof (struct heap64_struct)); 
  pq->n = 0;
  pq->heap_size = (heap_size < 2 ? 2: heap_size);
  pq->item = (hpq_item*) biomcmc_malloc ((pq->heap_size + 1) * sizeof (hpq_item));
  for (i = 0; i < pq->heap_size + 1; i++) { pq->item[i].freq = pq->item[i].id = 0; pq->item[i].hash = 0ULL; }
  return pq;
}

void
del_heap_hash64 (heap_hash64 pq)
{
  if (!pq) return;
  if (pq->item) free (pq->item); 
  free (pq); 
}

hpq_item
heap_hash64_get_maximum (heap_hash64 pq) { return pq->item[1]; }

hpq_item
heap_hash64_remove_maximum (heap_hash64 pq) // a.k.a. pop() 
{
  hpq_item max = pq->item[1]; //the root is the maximum element
  if (pq->n == 0) return (hpq_item) {0,0,0};
  pq->item[1] = pq->item[pq->n]; // replace the element at the top with last element in the heap
  pq->n--; 
  heap_hash64_bubble_down (pq, 1);
  return max;
}

void 
heap_hash64_insert (heap_hash64 pq, hpq_item item) // struct, not a pointer
{
  if (heap_hash64_item_already_here (pq, item, 1)) return; // if exists, then increase count  
  if (pq->n == pq->heap_size) { // replace if smaller than maximum
    if (item.hash < pq->item[1].hash) { 
      pq->item[1].id   = item.id; 
      pq->item[1].freq = item.freq;
      pq->item[1].hash = item.hash;
      heap_hash64_bubble_down (pq, 1); 
    }
  }
  else { // regular insert (a.k.a. push)
    pq->n++;
    pq->item[pq->n].id   = item.id;
    pq->item[pq->n].freq = item.freq; // should be one
    pq->item[pq->n].hash = item.hash;
    heap_hash64_bubble_up (pq, pq->n);
  }
}

static bool
heap_hash64_item_already_here (heap_hash64 pq, hpq_item item, int p)
{ // must traverse both children since MaxHeap is not ordered like BST
  if (p > pq->n) return false; // it can be equal, since we have one extra element (item[1...n] )
  if (item.hash == pq->item[p].hash) { pq->item[p].freq++; return true; }
  if (item.hash  > pq->item[p].hash) return false;
  bool below = heap_hash64_item_already_here (pq, item, 2 * p); // 'left' child 
  if (!below) below = heap_hash64_item_already_here (pq, item, 2 * p + 1); // traverse right child if not found yet
  return below;
}

static void
heap_hash64_bubble_down (heap_hash64 pq, int p)
{//bubbles down an element into it's proper place in the heap
	int i, c = 2*p, max_index = p; //c= index of younger child, max_index = index of maximum child
  for (i = 0; i < 2; i++) if(c + i <= pq->n) {
    //check to see if the data at min_index is smaller than the data at the child
    if (pq->item[max_index].hash < pq->item[c+i].hash)	max_index = c + i;
  }
	if (max_index != p) {
    hpq_item tmp = pq->item[p];
    pq->item[p] = pq->item[max_index];
    pq->item[max_index] = tmp;
    heap_hash64_bubble_down (pq, max_index);
  }
}

static void 
heap_hash64_bubble_up (heap_hash64 pq, int index)
{ //bubbles up the last element of the heap to maintain heap structure
  int parent = (index == 1 ? -1 : (int)(index/2));
  if (parent == -1) return;	//if we are at the root of the heap, no parent
  //if the parent node has a smaller data value, we need to bubble up
  if (pq->item[parent].hash < pq->item[index].hash) {
    hpq_item tmp = pq->item[parent];
    pq->item[parent] = pq->item[index];
    pq->item[index] = tmp;
    heap_hash64_bubble_up (pq, parent);
  }
}

void
heap_hash64_finalise_heap_qsort (heap_hash64 pq)
{
  hpq_item tmp = pq->item[0];  pq->item[0] = pq->item[pq->n];  pq->item[pq->n] = tmp;   // heap goes from [1...n] and we want [0...n-1]
//  for (int i=pq->n-4;i<pq->n-1;i++) printf ("%20lu %6d %3d > ", pq->item[i].hash, pq->item[i].id, pq->item[i].freq);
  qsort (pq->item, pq->n, sizeof (hpq_item), compare_hpq_item_increasing);
//  for (int i=pq->n-4;i<pq->n-1;i++) printf (" < %20lu %6d %3d", pq->item[i].hash, pq->item[i].id, pq->item[i].freq);
//  printf (":: DEBUG\n");
  if (pq->n < pq->heap_size - 1) { 
    pq->item = (hpq_item*) biomcmc_realloc ((hpq_item*) pq->item, pq->n * sizeof (hpq_item));
    pq->heap_size = pq->n; 
  }
  int sqrt_sum = 0;
  for (int i = 0; i < pq->n; i++) sqrt_sum += (pq->item[i].freq * pq->item[i].freq);
  pq->sqrt_sum = sqrt ((double)(sqrt_sum));
}

static int
compare_hpq_item_increasing (const void *a, const void *b)
{
  if ( ((hpq_item *)b)->hash < ((hpq_item *)a)->hash ) return 1;
  if ( ((hpq_item *)b)->hash > ((hpq_item *)a)->hash ) return -1;
  return 0;
}

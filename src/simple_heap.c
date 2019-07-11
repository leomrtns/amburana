#include "simple_heap.h" 

/* max heap structure for storing smallest hash values 
 * code inspired by https://gist.github.com/vgoel30/5d81e6abf9464930c1e126dab04d5be3  
 * and also for creating kNN of points for clustering (max heap as priority queue storing distances) */

static void heap64_bubble_down (heap64 h64, int p);
static void heap64_bubble_up (heap64 h64, int index);
static void heap_pqueue_bubble_down (heap_pqueue pq, int p);
static void heap_pqueue_bubble_up (heap_pqueue pq, int index);
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

heap_pqueue 
new_heap_pqueue (int heap_size)
{
  int i;
  heap_pqueue pq = (heap_pqueue) biomcmc_malloc (sizeof (struct heap64_struct)); 
  pq->n = 0;
  pq->heap_size = (heap_size < 2 ? 2: heap_size);
  pq->item = (hpq_item*) biomcmc_malloc ((pq->heap_size + 1) * sizeof (struct hpq_item_struct));
  for (i = 0; i < pq->heap_size + 1; i++) { pq->item[i]->id = 0; pq->item[i]->priority = 0.; pq->item[i]->v = NULL; }
  return pq;
}

void
del_heap_pqueue (heap_pqueue pq)
{
  if (pq) { if (pq->item) free (pq->item); free (pq); }
}

hpq_item
heap_pqueue_get_maximum (heap_pqueue pq) { return pq->item[1]; }

hpq_item
heap_pqueue_remove_maximum (heap_pqueue pq) // a.k.a. pop() 
{
  hpq_item max = pq->item[1]; //the root is the maximum element
  if (pq->n == 0) return NULL;
  pq->item[1] = pq->item[pq->n]; // replace the element at the top with last element in the heap
  pq->n--; 
  heap_pqueue_bubble_down (pq, 1);
  return max;
}

void 
heap_pqueue_insert (heap_pqueue pq, struct hpq_item_struct item) // struct, not a pointer
{
  if (pq->n == pq->heap_size) { // replace if smaller than maximum
    if (item.priority < pq->item[1]->priority) { 
      pq->item[1]->id = item.id; 
      pq->item[1]->priority = item.priority; 
      pq->item[1]->v = item.v; 
      heap_pqueue_bubble_down (pq, 1); 
    }
  }
  else { // regular insert (a.k.a. push)
    pq->n++;
    pq->item[pq->n]->id = item.id;
    pq->item[pq->n]->priority = item.priority;
    pq->item[pq->n]->v = item.v;
    heap_pqueue_bubble_up (pq, pq->n);
  }
}

static void
heap_pqueue_bubble_down (heap_pqueue pq, int p)
{//bubbles down an element into it's proper place in the heap
	int i, c = 2*p, max_index = p; //c= index of younger child, max_index = index of maximum child
  for (i = 0; i < 2; i++) if(c + i <= pq->n) {
    //check to see if the data at min_index is smaller than the data at the child
    if (pq->item[max_index]->priority < pq->item[c+i]->priority)	max_index = c + i;
  }
	if (max_index != p) {
    hpq_item tmp = pq->item[p];
    pq->item[p] = pq->item[max_index];
    pq->item[max_index] = tmp;
    heap_pqueue_bubble_down (pq, max_index);
  }
}

static void 
heap_pqueue_bubble_up (heap_pqueue pq, int index)
{ //bubbles up the last element of the heap to maintain heap structure
  int parent = (index == 1 ? -1 : (int)(index/2));
  if (parent == -1) return;	//if we are at the root of the heap, no parent
  //if the parent node has a smaller data value, we need to bubble up
  if (pq->item[parent]->priority < pq->item[index]->priority) {
    hpq_item tmp = pq->item[parent];
    pq->item[parent] = pq->item[index];
    pq->item[index] = tmp;
    heap_pqueue_bubble_up (pq, parent);
  }
}

void
heap_pqueue_finalise_heap_qsort (heap_pqueue pq)
{
  hpq_item tmp = pq->item[0];  pq->item[0] = pq->item[pq->n];  pq->item[pq->n] = tmp;   // heap goes from [1...n] and we want [0...n-1]
  qsort (pq->item, pq->n, sizeof (struct hpq_item_struct), compare_hpq_item_increasing);
  if (pq->n < pq->heap_size - 1)
    pq->item = (hpq_item*) biomcmc_realloc ((hpq_item) pq->item, pq->n * sizeof (struct hpq_item_struct));
  pq->heap_size = pq->n;
}

static int
compare_hpq_item_increasing (const void *a, const void *b)
{
  if ( ((hpq_item)b)->priority < ((hpq_item)a)->priority ) return 1;
  if ( ((hpq_item)b)->priority > ((hpq_item)a)->priority ) return -1;
  return 0;
}

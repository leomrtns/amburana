/**
 * Copyright 2014 Gagarine Yaikhom (MIT License)
 * https://github.com/gyaikhom/agglomerative-hierarchical-clustering
 * Implements Agglomerative Hierarchical Clustering algorithm.
 **/

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */

#define alloc_mem(N, T) (T *) calloc(N, sizeof(T))
#define alloc_fail(M) fprintf(stderr,                                   \
                              "Failed to allocate memory for %s.\n", M)
#define read_fail(M) fprintf(stderr, "Failed to read %s from file.\n", M)
#define invalid_node(I) fprintf(stderr,                                 \
                                "Invalid cluster node at index %d.\n", I)

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct neighbour_s neighbour_t;

float (*distance_fptr) (float **, const int *, const int *, int, int);

struct cluster_s {
  int num_items; /* number of items that was clustered */
  int num_clusters; /* current number of root clusters */
  int num_nodes; /* number of leaf and merged clusters */
  cluster_node_t *nodes; /* leaf and merged clusters */
  float **distances; /* distance between leaves */
};

struct cluster_node_s {
  int type; /* type of the cluster node */
  int is_root; /* true if cluster hasn't merged with another */
  int height; /* height of node from the bottom */
  coord_t centroid; /* centroid of this cluster */
  char *label; /* label of a leaf node */
  int *merged; /* indexes of root clusters merged */
  int num_items; /* number of leaf nodes inside new cluster */
  int *items; /* array of leaf nodes indices inside merged clusters */
  neighbour_t *neighbours; /* sorted linked list of distances to roots */
};

struct neighbour_s {
  int target; /* the index of cluster node representing neighbour */
  float distance; /* distance between the nodes */
  neighbour_t *next, *prev; /* linked list entries */
};

//typedef struct coord_s { float x, y; } coord_t;
//struct item_s { coord_t coord; char label[MAX_LABEL_LEN]; };


void
agglomerative_cluster_run ()
{
  cluster_t *cluster = agglomerate(num_items, items);
  get_k_clusters(cluster, k);
}


float **generate_distance_matrix(int num_items, const item_t items[])
{
  float **matrix = alloc_mem(num_items, float *);
  if (matrix) {
    for (int i = 0; i < num_items; ++i) {
      matrix[i] = alloc_mem(num_items, float);
      if (!matrix[i]) {
        alloc_fail("distance matrix row");
        num_items = i;
        for (i = 0; i < num_items; ++i) free(matrix[i]);
        free(matrix);
        matrix = NULL;
        break;
      }
    }
    if (matrix)  fill_euclidean_distances(matrix, num_items, items);
  } else alloc_fail("distance matrix");
  return matrix;
}

float single_linkage (float **distances, const int a[], const int b[], int m, int n)
{
  // LEO : should receive cluster_node_t (with centroid info etc) instead of numbers
  float min = FLT_MAX, d;
  for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) {
      d = distances[a[i]][b[j]];
      if (d < min)  min = d;
    }
  return min;
}

float complete_linkage (float **distances, const int a[], const int b[], int m, int n)
{
  float d, max = 0.0 /* assuming distances are positive */;
  for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) {
      d = distances[a[i]][b[j]];
      if (d > max)  max = d;
    }
  return max;
}
// Ward d(A,B) = sum(d(a,b))/(A+B) - sum(d(a1,a2))/A - sum(db1,b2))/B (see nearestneighbor chain on wiki)

float average_linkage (float **distances, const int a[], const int b[], int m, int n)
{
  float total = 0.0;
  for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) total += distances[a[i]][b[j]];
  return total / (m * n);
}

float get_distance (cluster_t *cluster, int index, int target)
{
  /* if both are leaves, just use the distances matrix */
  if (index < cluster->num_items && target < cluster->num_items)  return cluster->distances[index][target];
  else {
    cluster_node_t *a = &(cluster->nodes[index]);
    cluster_node_t *b = &(cluster->nodes[target]);
    return distance_fptr(cluster->distances, a->items, b->items, a->num_items, b->num_items);
  }
}

void free_neighbours (neighbour_t *node)
{
  neighbour_t *t;
  while (node) {
    t = node->next;
    free(node);
    node = t;
  }
}

void free_cluster_nodes(cluster_t *cluster)
{
  for (int i = 0; i < cluster->num_nodes; ++i) {
    cluster_node_t *node = &(cluster->nodes[i]);
    if (node->label) free(node->label);
    if (node->merged) free(node->merged);
    if (node->items) free(node->items);
    if (node->neighbours) free_neighbours(node->neighbours);
  }
  free(cluster->nodes);
}

void free_cluster(cluster_t * cluster)
{
  if (cluster) {
    if (cluster->nodes) free_cluster_nodes(cluster);
    if (cluster->distances) {
      for (int i = 0; i < cluster->num_items; ++i) free(cluster->distances[i]);
      free(cluster->distances);
    }
    free(cluster);
  }
}

void insert_before(neighbour_t *current, neighbour_t *neighbours, cluster_node_t *node)
{
  neighbours->next = current;
  if (current->prev) {
    current->prev->next = neighbours;
    neighbours->prev = current->prev;
  } 
  else node->neighbours = neighbours;
  current->prev = neighbours;
}

void insert_after(neighbour_t *current, neighbour_t *neighbours)
{
  neighbours->prev = current;
  current->next = neighbours;
}

void insert_sorted(cluster_node_t *node, neighbour_t *neighbours)
{
  neighbour_t *temp = node->neighbours;
  while (temp->next) {
    if (temp->distance >= neighbours->distance) {
      insert_before(temp, neighbours, node);
      return;
    }
    temp = temp->next;
  }
  if (neighbours->distance < temp->distance)  insert_before(temp, neighbours, node);
  else  insert_after(temp, neighbours);
}

neighbour_t *add_neighbour(cluster_t *cluster, int index, int target)
{
  neighbour_t *neighbour = alloc_mem(1, neighbour_t);
  if (neighbour) {
    neighbour->target = target;
    neighbour->distance = get_distance(cluster, index, target);
    cluster_node_t *node = &(cluster->nodes[index]);
    if (node->neighbours) insert_sorted(node, neighbour);
    else node->neighbours = neighbour;
  } else alloc_fail("neighbour node");
  return neighbour;
}

cluster_t *update_neighbours(cluster_t *cluster, int index)
{
  cluster_node_t *node = &(cluster->nodes[index]);
  if (node->type == NOT_USED) {
    invalid_node(index);
    cluster = NULL;
  } else {
    int root_clusters_seen = 1, target = index;
    while (root_clusters_seen < cluster->num_clusters) {
      cluster_node_t *temp = &(cluster->nodes[--target]);
      if (temp->type == NOT_USED) {
        invalid_node(index);
        cluster = NULL;
        break;
      }
      if (temp->is_root) {
        ++root_clusters_seen;
        add_neighbour(cluster, index, target);
      }
    }
  }
  return cluster;
}

#define init_leaf(cluster, node, item, len)             \
  do {                                            \
    strncpy(node->label, item->label, len); \
    node->centroid = item->coord;           \
    node->type = LEAF_NODE;                 \
    node->is_root = 1;                      \
    node->height = 0;                       \
    node->num_items = 1;                    \
    node->items[0] = cluster->num_nodes++;  \
  } while (0)                                     \

cluster_node_t *add_leaf(cluster_t *cluster, const item_t *item)
{
  cluster_node_t *leaf = &(cluster->nodes[cluster->num_nodes]);
  int len = strlen(item->label) + 1;
  leaf->label = alloc_mem(len, char);
  if (leaf->label) {
    leaf->items = alloc_mem(1, int);
    if (leaf->items) {
      init_leaf(cluster, leaf, item, len);
      cluster->num_clusters++;
    } else {
      alloc_fail("node items");
      free(leaf->label);
      leaf = NULL;
    }
  } else {
    alloc_fail("node label");
    leaf = NULL;
  }
  return leaf;
}

#undef init_leaf

cluster_t *add_leaves(cluster_t *cluster, item_t *items)
{
  for (int i = 0; i < cluster->num_items; ++i) {
    if (add_leaf(cluster, &items[i])) update_neighbours(cluster, i);
    else { cluster = NULL; break; }
  }
  return cluster;
}

void print_cluster_items(cluster_t *cluster, int index)
{
  cluster_node_t *node = &(cluster->nodes[index]);
  fprintf(stdout, "Items: ");
  if (node->num_items > 0) {
    fprintf(stdout, "%s", cluster->nodes[node->items[0]].label);
    for (int i = 1; i < node->num_items; ++i) fprintf(stdout, ", %s", cluster->nodes[node->items[i]].label);
  }
  fprintf(stdout, "\n");
}

void print_cluster_node(cluster_t *cluster, int index)
{
  cluster_node_t *node = &(cluster->nodes[index]); 
  fprintf(stdout, "Node %d - height: %d, centroid: (%5.3f, %5.3f)\n", index, node->height, node->centroid.x, node->centroid.y);
  if (node->label) fprintf(stdout, "\tLeaf: %s\n\t", node->label);
  else fprintf(stdout, "\tMerged: %d, %d\n\t", node->merged[0], node->merged[1]);
  print_cluster_items(cluster, index);
  fprintf(stdout, "\tNeighbours: ");
  neighbour_t *t = node->neighbours;
  while (t) {
    fprintf(stdout, "\n\t\t%2d: %5.3f", t->target, t->distance);
    t = t->next;
  }
  fprintf(stdout, "\n");
}

void merge_items(cluster_t *cluster, cluster_node_t *node, cluster_node_t **to_merge)
{
  node->type = A_MERGER;
  node->is_root = 1;
  node->height = -1;
// LEO: main function to modify (replace centroid by ward's avge distance^2 between points etc) 
  /* copy leaf indexes from merged clusters */
  int k = 0, idx;
  coord_t centroid = { .x = 0.0, .y = 0.0 };
  for (int i = 0; i < 2; ++i) {
    cluster_node_t *t = to_merge[i];
    t->is_root = 0; /* no longer root: merged */
    if (node->height == -1 || node->height < t->height) node->height = t->height;
    for (int j = 0; j < t->num_items; ++j) {
      idx = t->items[j];
      node->items[k++] = idx;
    }
    centroid.x += t->num_items * t->centroid.x;
    centroid.y += t->num_items * t->centroid.y;
  }
  /* calculate centroid */
  node->centroid.x = centroid.x / k;
  node->centroid.y = centroid.y / k;
  node->height++;
}

#define merge_to_one(cluster, to_merge, node, node_idx)         \
  do {                                                    \
    node->num_items = to_merge[0]->num_items +      \
    to_merge[1]->num_items;                 \
    node->items = alloc_mem(node->num_items, int);  \
    if (node->items) {                              \
      merge_items(cluster, node, to_merge);   \
      cluster->num_nodes++;                   \
      cluster->num_clusters--;                \
      update_neighbours(cluster, node_idx);   \
    } else {                                        \
      alloc_fail("array of merged items");    \
      free(node->merged);                     \
      node = NULL;                            \
    }                                               \
  } while(0)                                              \

cluster_node_t *merge(cluster_t *cluster, int first, int second)
{
  int new_idx = cluster->num_nodes;
  cluster_node_t *node = &(cluster->nodes[new_idx]);
  node->merged = alloc_mem(2, int);
  if (node->merged) {
    cluster_node_t *to_merge[2] = {
      &(cluster->nodes[first]),
      &(cluster->nodes[second])
    };
    node->merged[0] = first;
    node->merged[1] = second;
    merge_to_one(cluster, to_merge, node, new_idx);
  } else {
    alloc_fail("array of merged nodes");
    node = NULL;
  }
  return node;
}

#undef merge_to_one

void find_best_distance_neighbour(cluster_node_t *nodes, int node_idx, neighbour_t *neighbour, float *best_distance, int *first, int *second)
{
  while (neighbour) {
    if (nodes[neighbour->target].is_root) {
      if (*first == -1 || neighbour->distance < *best_distance) {
        *first = node_idx;
        *second = neighbour->target;
        *best_distance = neighbour->distance;
      }
      break;
    }
    neighbour = neighbour->next;
  }
}


int find_clusters_to_merge(cluster_t *cluster, int *first, int *second)
{
  float best_distance = 0.0;
  int root_clusters_seen = 0;
  int j = cluster->num_nodes; /* traverse hierarchy top-down */
  *first = -1;
  while (root_clusters_seen < cluster->num_clusters) {
    cluster_node_t *node = &(cluster->nodes[--j]);
    if (node->type == NOT_USED || !node->is_root)  continue;
    ++root_clusters_seen;
    find_best_distance_neighbour(cluster->nodes, j, node->neighbours, &best_distance, first, second);
  }
  return *first;
}

cluster_t *merge_clusters(cluster_t *cluster)
{
  int first, second;
  while (cluster->num_clusters > 1) {
    if (find_clusters_to_merge(cluster, &first, &second) != -1) merge(cluster, first, second);
  }
  return cluster;
}

#define init_cluster(cluster, num_items, items)                         \
  do {                                                            \
    cluster->distances =                                    \
    generate_distance_matrix(num_items, items);     \
    if (!cluster->distances)                                \
    goto cleanup;                                   \
    cluster->num_items = num_items;                         \
    cluster->num_nodes = 0;                                 \
    cluster->num_clusters = 0;                              \
    if (add_leaves(cluster, items))                         \
    merge_clusters(cluster);                        \
    else                                                    \
    goto cleanup;                                   \
  } while (0)                                                     \

cluster_t *agglomerate(int num_items, item_t *items)
{
  cluster_t *cluster = alloc_mem(1, cluster_t);
  if (cluster) {
    cluster->nodes = alloc_mem(2 * num_items - 1, cluster_node_t);
    if (cluster->nodes) init_cluster(cluster, num_items, items);
    else {
      alloc_fail("cluster nodes");
      goto cleanup;
    }
  } else alloc_fail("cluster");
  goto done;

cleanup:
  free_cluster(cluster);
  cluster = NULL;

done:
  return cluster;
}

#undef init_cluster

int print_root_children(cluster_t *cluster, int i, int nodes_to_discard)
{
  cluster_node_t *node = &(cluster->nodes[i]);
  int roots_found = 0;
  if (node->type == A_MERGER) {
    for (int j = 0; j < 2; ++j) {
      int t = node->merged[j];
      if (t < nodes_to_discard) {
        print_cluster_items(cluster, t);
        ++roots_found;
      }
    }
  }
  return roots_found;
}

void get_k_clusters(cluster_t *cluster, int k)
{
  if (k < 1)
    return;
  if (k > cluster->num_items)
    k = cluster->num_items;

  int i = cluster->num_nodes - 1;
  int roots_found = 0;
  int nodes_to_discard = cluster->num_nodes - k + 1;
  while (k) {
    if (i < nodes_to_discard) {
      print_cluster_items(cluster, i);
      roots_found = 1;
    } else
      roots_found = print_root_children(cluster, i, nodes_to_discard);
    k -= roots_found;
    --i;
  }
}

void print_cluster(cluster_t *cluster)
{
  for (int i = 0; i < cluster->num_nodes; ++i) print_cluster_node(cluster, i);
}

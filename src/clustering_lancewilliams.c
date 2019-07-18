/**
 * Copyright 2014 Gagarine Yaikhom (MIT License)
 * https://github.com/gyaikhom/agglomerative-hierarchical-clustering
 * Implements Agglomerative Hierarchical Clustering algorithm.
 **/

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */

#define alloc_mem(N, T) (T *) calloc(N, sizeof(T))

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct neighbour_s neighbour_t;

struct cluster_s {
  int num_items; /* number of items that was clustered */
  int num_clusters; /* current number of root clusters */
  int num_nodes; /* number of leaf and merged clusters */
  cluster_node_t *nodes; /* leaf and merged clusters */
};

struct cluster_node_s {
  int type; /* type of the cluster node */
  int is_root; /* true if cluster hasn't merged with another */
  int height; /* height of node from the bottom */
  coord_t centroid; /* centroid of this cluster */
  char *label; /* label of a leaf node */
  int *merged; /* indexes of root clusters merged */
  int n_items; /* number of leaf nodes inside new cluster */
  int *item; /* array of leaf nodes indices inside merged clusters */
  neighbour_t *neighbours; /* sorted linked list of distances to roots */
};

struct neighbour_s {
  int target; /* the index of cluster node representing neighbour */
  double distance; /* distance between the nodes */
  double variance; /* LEO: store info about new cluster (Ward variance for instance)  */
  neighbour_t *next, *prev; /* linked list entries */
};

agglomerative_cluster
new_agglomerative_cluster ()
{
  ag->cluster = malloc(cluster_t);
}

void
agglomerative_cluster_run ()
{
  ag->cluster_t *cluster = agglomerate(num_items, items);
  get_k_clusters(cluster, k);
}

double get_distance (agglomerative_cluster ag, int index, int target)
{
  /* if both are leaves, just use the distances matrix */
  if (index < ag->cluster->num_items && target < ag->cluster->num_items) return distance_generator_get (ag->d, index, target);
  return ag->distance_fptr (ag, a->items, ag->cluster->nodes[index], ag->cluster->nodes[target]);
}

double single_linkage (agglomerative_cluster ag,  cluster_node_t *a, cluster_node_t *b)
{
  // LEO : should receive cluster_node_t (with centroid info etc) instead of numbers
  double min = DBL_MAX, d;
  for (int i = 0; i < a->n_items; ++i) for (int j = 0; j < b->n_items; ++j) {
    d = distance_generator_get (ag->d, a->item[i], b->item[j]);
    if (d < min)  min = d;
  }
  return min;
}

double complete_linkage (agglomerative_cluster ag,  cluster_node_t *a, cluster_node_t *b)
{
  double d, max = -DBL_MAX;
  for (int i = 0; i < a->n_items; ++i) for (int j = 0; j < b->n_items; ++j) {
    d = distance_generator_get (ag->d, a->item[i], b->item[j]);
    if (d > max)  max = d;
  }
  return max;
}

double average_linkage (agglomerative_cluster ag,  cluster_node_t *a, cluster_node_t *b)
{
  double d = 0.0;
  for (int i = 0; i < a->n_items; ++i) for (int j = 0; j < b->n_items; ++j) d +=  distance_generator_get (ag->d, a->item[i], b->item[j]);
  return d/ (a->n_items * b->n_items);
}

double ward_linkage (agglomerative_cluster ag,  cluster_node_t *a, cluster_node_t *b)
{
  double tmp_d, d = 0.0;
  for (int i = 0; i < a->n_items; ++i) for (int j = 0; j < b->n_items; ++j) {
    tmp_d = distance_generator_get (ag->d, a->item[i], b->item[j]);
    d += (tmp_d * tmp_d);
  }
  return  d / (a->n_items + b->n_items) - a->variance/a->n_items - b->variance/b->n_items;
  // save d somewhere, if accepted we update it to d + a->var + b->var
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

void find_best_distance_neighbour(cluster_node_t *nodes, int node_idx, neighbour_t *neighbour, double *best_distance, int *first, int *second)
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
  double best_distance = 0.0;
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


/**
 * Copyright 2014 Gagarine Yaikhom (MIT License)
 * https://github.com/gyaikhom/agglomerative-hierarchical-clustering
 * Implements Agglomerative Hierarchical Clustering algorithm.
 **/

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
  bool is_leaf; /* leaf or merged */ 
  bool is_root; /* true if cluster hasn't merged with another */
  double height[2], length; /* height of node from the bottom (= distance) and to lower node */
  int merged[2]; /* indexes of root clusters merged */
  int n_items; /* number of leaf nodes inside new cluster */
  int *item; /* array of leaf nodes indices inside merged clusters */
  neighbour_t *neighbours; /* sorted linked list of distances to roots */
};

struct neighbour_s {
  int target; /* the index of cluster node representing neighbour */
  double distance[2]; /* distance between the nodes and extra info (variance for Ward's mehod) */
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

double get_distance_into_neighbour (agglomerative_cluster ag, neighbour_t neighbour, int index, int target)
{
  /* if both are leaves, just use the distances matrix */
  // TODO update distance[2] (not just return) 
  if (index < ag->cluster->num_items && target < ag->cluster->num_items) return distance_generator_get (ag->d, index, target);
  return ag->distance_fptr (ag, a->items, ag->cluster->nodes[index], ag->cluster->nodes[target]);
}

double single_linkage (agglomerative_cluster ag,  cluster_node_t *a, cluster_node_t *b)
{
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
    if (temp->distance[0] >= neighbours->distance[0]) {
      insert_before(temp, neighbours, node);
      return;
    }
    temp = temp->next;
  }
  if (neighbours->distance[0] < temp->distance[0]) insert_before(temp, neighbours, node);
  else  insert_after(temp, neighbours);
}

neighbour_t *add_neighbour(cluster_t *cluster, int index, int target)
{
  neighbour_t *neighbour = alloc_mem (1, neighbour_t);
  neighbour->target = target;
  get_distance_into_neighbour (cluster, neighbour, index, target);
  cluster_node_t *node = &(cluster->nodes[index]);
  if (node->neighbours) insert_sorted(node, neighbour);  // LEO: only called here
  else node->neighbours = neighbour;
  return neighbour;
}

cluster_t *update_neighbours(cluster_t *cluster, int index)
{
  cluster_node_t *node = &(cluster->nodes[index]);
  int root_clusters_seen = 1, target = index;
  while (root_clusters_seen < cluster->num_clusters) {
    cluster_node_t *temp = &(cluster->nodes[--target]);
    if (temp->is_root) {
      ++root_clusters_seen;
      add_neighbour(cluster, index, target);  // LEO: only called here
    }
  }
  return cluster;
}

void merge_items (cluster_t *cluster, cluster_node_t *node, int ids[2])
{
  int k = 0, idx;
  node->is_leaf = false;
  node->is_root = true;

  for (int i = 0; i < 2; ++i) {
    cluster_node_t *t = &(cluster->nodes[ids[i]]);
    t->is_root = false; /* no longer root: merged */
    for (int j = 0; j < t->num_items; ++j) { /* copy leaf indexes from merged clusters */
      idx = t->items[j];
      node->items[k++] = idx;
    }
    t->length = (node->height[0] - t->height[0]) / 2.; // branch length is half the difference in cluster distances
  }
}

cluster_node_t *merge (cluster_t *cluster, int ids[2], double best_distance[2])
{
  int new_idx = cluster->num_nodes;
  cluster_node_t *node = &(cluster->nodes[new_idx]);
  node->merged[0] = ids[0];
  node->merged[1] = ids[1];
  node->height[0] = best_distance[0]; // cluster distance 
  node->height[1] = best_distance[1]; // extra info (variance in Ward, or avge distance in upgma)
  node->num_items = cluster->nodes[ids[0]].num_items + cluster->nodes[ids[1]].num_items;
  node->items = alloc_mem(node->num_items, int);
  merge_items (cluster, node, ids);
  cluster->num_nodes++;
  cluster->num_clusters--;
  update_neighbours(cluster, node_idx);
  return node;
}

cluster_t *agglomerate(int num_items, item_t *items)
{
  int ids[2];
  double best_distance[2] = {DBL_MAX, 0.};
  cluster_t *cluster = alloc_mem(1, cluster_t);
  cluster->nodes = alloc_mem(2 * num_items - 1, cluster_node_t);
  cluster->distances =  generate_distance_matrix(num_items, items);
  cluster->num_items = num_items;
  cluster->num_nodes = 0;
  cluster->num_clusters = 0;

  for (int i = 0; i < cluster->num_items; ++i) if (add_leaf (cluster, &items[i])) update_neighbours(cluster, i); 

  while (cluster->num_clusters > 1) if (find_clusters_to_merge(cluster, ids, best_distance) != -1) merge (cluster, ids, best_distance);
  return cluster;
}

cluster_node_t *add_leaf (cluster_t *cluster, const item_t *item)
{
  cluster_node_t *leaf = &(cluster->nodes[cluster->num_nodes]);
  leaf->items = alloc_mem(1, int);
  node->label_id = cluster->num_nodes; // redundant, same as items[0] 
  node->is_leaf = node->is_root = true;
  node->height[0] = node->height[1] = node->length  = 0.;
  node->num_items = 1;
  node->items[0] = cluster->num_nodes++;
  cluster->num_clusters++;
  return leaf;
}

int find_clusters_to_merge (cluster_t *cluster, int ids[2], double best_distance[2])
{
  int root_clusters_seen = 0;
  int j = cluster->num_nodes; /* traverse hierarchy top-down */
  ids[0] = -1;
  best_distance[0] = DBL_MAX; // maybe zero is fine, due to "-1"
  while (root_clusters_seen < cluster->num_clusters) {
    cluster_node_t *node = &(cluster->nodes[--j]);
    if (!node->is_root)  continue;
    ++root_clusters_seen;
    find_best_distance_neighbour (cluster->nodes, j, node->neighbours, ids, best_distance);
  }
  return ids[0];
}

void find_best_distance_neighbour (cluster_node_t *nodes, int node_idx, neighbour_t *neighbour, int ids[2], double best_distance[2])
{
  while (neighbour) {
    if (nodes[neighbour->target].is_root) {
      if (ids[0] == -1 || neighbour->distance[0] < best_distance[0]) {
        ids[0] = node_idx;
        ids[1] = neighbour->target;
        best_distance[0] = neighbour->distance[0];
        best_distance[1] = neighbour->distance[1];
      }
      break;
    }
    neighbour = neighbour->next;
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

int print_root_children(cluster_t *cluster, int i, int nodes_to_discard)
{
  cluster_node_t *node = &(cluster->nodes[i]);
  int roots_found = 0;
  if (!node->is_leaf) {
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
    fprintf(stdout, "%3d", cluster->nodes[node->items[0]].label_id);
    for (int i = 1; i < node->num_items; ++i) fprintf(stdout, ", %3d", cluster->nodes[node->items[i]].label_id);
  }
  fprintf(stdout, "\n");
}

void print_cluster_node(cluster_t *cluster, int index)
{
  cluster_node_t *node = &(cluster->nodes[index]); 
  fprintf(stdout, "Node %d - height: %d, centroid: (%5.3f, %5.3f)\n", index, node->height, node->centroid.x, node->centroid.y);
  if (node->is_leaf) fprintf(stdout, "\tLeaf: %3d\n\t", node->label_id);
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


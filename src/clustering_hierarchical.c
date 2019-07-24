/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "clustering_hierarchical.h"

typedef struct hierarchical_cluster_struct* hierarchical_cluster;
typedef struct hc_link hc_link_t;
typedef struct hc_cluster hc_cluster_t;
typedef struct hc_level hc_level_t;

struct hierarchical_cluster_struct {
  distance_generator d;
  hc_cluster_t *clusters;
  hc_link_t *links;
  hc_level_t *levels;
  char linkage;
};

struct hc_link {
  int index_l;
  int count;
  double distance; // extra info from ward goes here
  hc_cluster_t *source;
  hc_cluster_t *target;
};

struct hc_cluster {
  int index_c;
  int sample_id;
  int offset;
  int final_level;
  int removed;
  hc_link_t *link;
  hc_link_t *minimum;
  hc_cluster_t *prev;
  hc_cluster_t *next;
  hc_cluster_t *cluster;
};

struct hc_level {
  double linkage;
  int id, source_id, target_id;
  hc_link_t *link;
  hc_level_t *next;
};

void hc_cluster_init_zero (hc_cluster_t *clusters, int length);
hierarchical_cluster new_hierarchical_cluster (distance_generator dg);
void del_hierarchical_cluster (hierarchical_cluster ac);

hc_link_t *hc_link_min (hc_link_t *a, hc_link_t* b);
double hc_average_linkage (hc_link_t *source, hc_link_t *target);
hc_link_t* hc_merge_link (hc_link_t *source, hc_link_t *target, hierarchical_cluster hac);
hc_link_t* hc_update (hc_cluster_t *cluster, hc_link_t *links);
hc_link_t* hc_next (hierarchical_cluster hac, hc_cluster_t *cluster, hc_link_t *links, hc_link_t *lastLink);
void hc_unlink (hc_cluster_t *cluster);
void hc_merge_cluster (hc_cluster_t *source, hc_cluster_t *target);

void
hc_cluster_init_zero (hc_cluster_t *clusters, int length) 
{
  hc_cluster_t *curr;
  hc_cluster_t *prev;

  curr = &clusters[0];
  curr->index_c = curr->offset = curr->final_level = curr->removed  = 0;
  curr->minimum = curr->link = NULL;
  curr->next = curr->prev = curr->cluster = NULL;
  curr->sample_id = -1;
  for (int i = 0; i < length; i++) {
    prev = &clusters[i];
    curr = &clusters[i + 1];
    curr->offset = curr->removed  = 0;
    curr->minimum = curr->link = NULL;
    curr->next = curr->cluster = NULL;
    curr->final_level = curr->index_c = curr->sample_id = i;
    curr->prev = prev;
    prev->next = curr;
  }
}
 
hierarchical_cluster
new_hierarchical_cluster (distance_generator dg)
{
  int i, length, offset = 0, total = 0;
  hc_cluster_t *curr, *prev;
  hc_link_t *thislink;

  hierarchical_cluster hac = (hierarchical_cluster) biomcmc_malloc (sizeof (struct hierarchical_cluster_struct));
  hac->d = dg; dg->ref_counter++;

  hac->clusters = (hc_cluster_t*) biomcmc_malloc ((dg->n_samples + 1) * sizeof(hc_cluster_t));
  hc_cluster_init_zero (hac->clusters, dg->n_samples);
  curr = hac->clusters;

  for (length = 0; (curr = curr->next); length++) total += length;
  hac->links  = (hc_link_t*) biomcmc_malloc (total * sizeof(hc_link_t));
  hac->levels = (hc_level_t*) biomcmc_malloc (hac->d->n_samples * sizeof (hc_level_t));

  for (i = 0; i < length; i++) {
    prev = &(hac->clusters)[i];
    curr = &(hac->clusters)[i + 1];
    curr->offset = offset;

    while (prev->prev != NULL) {
      thislink = &(hac->links)[curr->offset + prev->index_c];
      thislink->source = prev;
      thislink->target = curr;
      thislink->count = 1;
      thislink->distance = distance_generator_get (dg, prev->index_c, curr->index_c);
      curr->minimum = NULL;
      if (curr->minimum == NULL || thislink->distance < curr->minimum->distance) curr->minimum = thislink;
      prev = prev->prev;
      offset++;
    }
  }
  hac->linkage = 'u';
  return hac;
}

void
del_hierarchical_cluster (hierarchical_cluster hac)
{
  if (!hac) return;
  if (hac->clusters) free (hac->clusters);
  if (hac->links) free (hac->links);
  if (hac->levels) free (hac->levels);
  del_distance_generator (hac->d);
  free (hac);
}

topology
hierarchical_cluster_topology (distance_generator dg, char linkage)
{
  hierarchical_cluster hac = new_hierarchical_cluster (dg);
  hc_cluster_t *target = NULL, *source = NULL, *cluster = NULL;
  hc_link_t *thislink = NULL;
  hc_level_t *level, *prev = NULL;
  int n = hac->d->n_samples, i, offset = 0;
  topology tree = new_topology (dg->n_samples);

  switch (linkage) {
    case ('i'): case ('I'): hac->linkage = 'i'; break; // single linkage
    case ('a'): case ('A'): hac->linkage = 'a'; break; // complete linkage
    case ('w'): case ('W'): hac->linkage = 'w'; break; // weighted linkage (WPGMA)
    case ('u'): case ('U'): hac->linkage = 'u'; break; // average linkage (UPGMA)
    default:                hac->linkage = 'u'; break;
  };

  cluster = hac->clusters;
  for (i = 0; i < n; i++) tree->blength[i] = 0.; // initially holds height (a.k.a. distance, linkage)

  while ((thislink = hc_next (hac, hac->clusters, hac->links, thislink)) != NULL) {
    target = thislink->target;
    source = thislink->source;
    thislink->index_l = offset;
    level = &(hac->levels)[offset];
    level->id = hac->d->n_samples + offset++; // node number in tree
    if (prev != NULL) prev->next = level;
    level->next      = NULL;
    level->linkage   = thislink->distance;
    level->source_id = source->sample_id;
    level->target_id = target->sample_id;
    level->link      = thislink;
    if (source->prev == NULL) cluster = source->next;
    source->link = thislink;
    hc_merge_cluster (source, target);
    
    tree->blength[level->id] = level->linkage; // branch lengths provisionally have cluster distances (branch lengths are a difference in distances)
    prev = level;
  }

  /* store proper indexes into levels (in order) */
  for (level = hac->levels; level; level = level->next) {
    int s_i = level->source_id + 1; // + 1 since cluster[0] is not a leaf 
    int t_i = level->target_id + 1; 
    level->source_id = (&cluster[s_i])->final_level;
    level->target_id = (&cluster[t_i])->final_level;
    tree->blength[level->source_id] = (level->linkage  - tree->blength[level->source_id])/2.;
    tree->blength[level->target_id] = (level->linkage  - tree->blength[level->target_id])/2.;
    (&cluster[s_i])->final_level = level->id; // leaves "climb up" the tree, so that the store ID of current LCA 
    (&cluster[t_i])->final_level = level->id; 
    //printf("DEBUG 2 :: %2d %2d    %2d  %8.5lf\n", level->source_id, level->target_id, level->id,  level->linkage);
  }

  tree->blength[2*n-2] = 0.;  // still has overall distance
  /* topology just needs to updates pointers up, right, and left (sisters are updated later) */
  for (level = hac->levels; level; level = level->next) { // for each internal node
    i = level->id;
    tree->nodelist[i]->left  = tree->nodelist[ level->source_id ];
    tree->nodelist[i]->right = tree->nodelist[ level->target_id ];
    tree->nodelist[i]->right->up = tree->nodelist[i]->left->up = tree->nodelist[i];
  }
  
  update_topology_sisters (tree);
  update_topology_traversal (tree);
  del_hierarchical_cluster (hac);
  return tree;
}

hc_link_t 
*hc_link_min (hc_link_t *a, hc_link_t* b) 
{
  if (a == NULL) return b;
  if (b == NULL) return a;
  if (a->distance < b->distance) return a;
  return b;
}

double
hc_average_linkage (hc_link_t *source, hc_link_t *target) // UPGMA
{ // target is the one with new values (see += below...)
  target->count += source->count;
  return (double)((source->distance * source->count) + (target->distance * target->count)) / (double)(target->count);
}

double
hc_weighted_linkage (hc_link_t *source, hc_link_t *target) 
{ // looks less weighted than above but it means that _nodes_ end up being weighted (UPGMA gives same weight to all nodes)
  return (source->distance + target->distance) / 2.;
}

double
hc_minimum_linkage (hc_link_t *source, hc_link_t *target) 
{
  if (source->distance < target->distance) return source->distance;
  return target->distance;
}

double
hc_maximum_linkage (hc_link_t *source, hc_link_t *target) 
{
  if (source->distance > target->distance) return source->distance;
  return target->distance;
}

hc_link_t* 
hc_merge_link (hc_link_t *source, hc_link_t *target, hierarchical_cluster hac) 
{ 
  switch (hac->linkage) {
    case ('i'): target->distance = hc_minimum_linkage (source, target); break;
    case ('a'): target->distance = hc_maximum_linkage (source, target); break;
    case ('w'): target->distance = hc_weighted_linkage (source, target); break;
    case ('u'):
    default:    target->distance = hc_average_linkage (source, target); break;
  };
  return target; 
}

hc_link_t* 
hc_update (hc_cluster_t *cluster, hc_link_t *links) 
{
  hc_link_t *min = NULL, *link;
  int offset = cluster->offset;
  while ((cluster = cluster->prev) != NULL) {
    link = &links[offset + cluster->index_c];
    if (link->source->removed == 0 && (min == NULL || link->distance < min->distance)) min = link;
  }
  return min;
}

hc_link_t*
hc_next (hierarchical_cluster hac, hc_cluster_t *cluster, hc_link_t *links, hc_link_t *lastLink) 
{
  hc_link_t *link, *mergedlink = NULL, *minLink = NULL;

  while (cluster != NULL) {
    if (cluster->minimum != NULL && cluster->minimum->source->removed == 1) cluster->minimum = hc_update (cluster, links);
    if (lastLink != NULL) {
      hc_cluster_t *source = lastLink->source;
      hc_cluster_t *target = lastLink->target;
      if (cluster->index_c > target->index_c) {
        mergedlink = hc_merge_link (&links[cluster->offset + source->index_c], &links[cluster->offset + target->index_c], hac);
        cluster->minimum = hc_link_min (mergedlink, cluster->minimum);
      } else if (cluster->index_c < target->index_c) {
        link = &links[source->offset + cluster->index_c];
        if (cluster->index_c > source->index_c) link = &links[cluster->offset + source->index_c];
        mergedlink = hc_merge_link (link, &links[target->offset + cluster->index_c], hac);
        target->minimum = hc_link_min (mergedlink, target->minimum);
      }
    }
    if (cluster->minimum != NULL && (minLink == NULL || cluster->minimum->distance < minLink->distance))  minLink = cluster->minimum;
    cluster = cluster->next;
  }
  return minLink;
}

void 
hc_unlink (hc_cluster_t *cluster) 
{
  cluster->removed = 1;
  cluster->next->prev = cluster->prev;
  if (cluster->prev != NULL) cluster->prev->next = cluster->next;
}

void 
hc_merge_cluster (hc_cluster_t *source, hc_cluster_t *target) 
{
  while (target->cluster != NULL) target = target->cluster;
  target->cluster = source;
  hc_unlink (source);
}

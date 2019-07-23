/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "clustering_hierarchical.h"

typedef struct hc_link hc_link_t;
typedef struct hc_cluster hc_cluster_t;
typedef struct hc_level hc_level_t;

struct hc_link {
  int index_l;
  int count;
//  int *clusters;
  double distance; // extra info from ward goes here
  hc_cluster_t *source;
  hc_cluster_t *target;
};

struct hc_cluster {
  int index_c;
  int offset;
  int length;
  int removed;
  hc_link_t *link;
  hc_link_t *minimum;
  hc_cluster_t *prev;
  hc_cluster_t *next;
  hc_cluster_t *cluster;
};

struct hc_level {
  double linkage;
  int source_id;
  int target_id;
  hc_link_t *link;
  hc_level_t *next;
  int *clusters;
};

void hc_cluster_init_zero (hc_cluster_t *clusters, int length);
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
  curr->index_c = curr->offset = curr->length = curr->removed  = 0;
  curr->minimum = curr->link = NULL;
  curr->next = curr->prev = curr->cluster = NULL;
  for (int i = 0; i < length; i++) {
    prev = &clusters[i];
    curr = &clusters[i + 1];
    curr->offset = curr->length = curr->removed  = 0;
    curr->minimum = curr->link = NULL;
    curr->next = curr->cluster = NULL;
    curr->index_c   = i; 
    curr->prev      = prev;
    prev->next      = curr;
  }
}
 
hierarchical_cluster
new_hierarchical_cluster (distance_generator dg)
{
  int i, length, offset = 0, total = 0;
  hc_cluster_t *curr, *clusters, *prev;
  hc_link_t *thislink, *links;

  hierarchical_cluster hac = (hierarchical_cluster) biomcmc_malloc (sizeof (struct hierarchical_cluster_struct));
  hac->d = dg; dg->ref_counter++;

  clusters = (hc_cluster_t*) biomcmc_malloc((dg->n_samples + 1) * sizeof(hc_cluster_t));
  hc_cluster_init_zero (clusters, dg->n_samples);
  curr = clusters;

  for (length = 0; (curr = curr->next); length++) total += length;
  links = (hc_link_t*) biomcmc_malloc (total * sizeof(hc_link_t));

  for (i = 0; i < length; i++) {
    prev = &clusters[i];
    curr = &clusters[i + 1];
    curr->offset = offset;

    while (prev->prev != NULL) {
      thislink = &links[curr->offset + prev->index_c];
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
  hac->clusters = (void*) clusters;
  hac->links    = (void*) links;
  hac->levels   = NULL;
  hac->linkage = 'u';
  hac->timing_secs = 0.;
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
}

void
hierarchical_cluster_run (hierarchical_cluster hac, char linkage)
{
  hc_cluster_t *target = NULL, *source = NULL, *cluster;
  hc_link_t *thislink = NULL, *links;
  hc_level_t *levels = (hc_level_t*) biomcmc_malloc (hac->d->n_samples * sizeof(hc_level_t));
  hc_level_t *level, *prev = NULL;
  int offset = 0;
  clock_t time1, time0;

  cluster = (hc_cluster_t*) hac->clusters;
  links = (hc_link_t*) hac->links;
  switch (linkage) {
    case ('u'):
    case ('U'): hac->linkage = 'u'; break;
    default:    hac->linkage = 'u'; break;
  };
  time0 = clock ();
  while ((thislink = hc_next (hac, cluster, links, thislink)) != NULL) {
    target = thislink->target;
    source = thislink->source;
    thislink->index_l = offset;
    level = &levels[offset++];

    // link->clusters = (int*) malloc((length - offset) * sizeof(int) * 2);
    if (prev != NULL) prev->next = level;
    level->next      = NULL;
    level->linkage   = thislink->distance;
    level->source_id = source->index_c;
    level->target_id = target->index_c;
    // level->clusters = link->clusters;
    if (source->prev == NULL) cluster = source->next;
    source->link = thislink;
    hc_merge_cluster (source, target);
    prev = level;
  }
  time1 = clock (); hac->timing_secs += (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  hac->levels = levels;
  for (level = levels; level; level = level->next) {  
    printf("DEBUG:: %2d %2d %8.5lf\n", level->source_id, level->target_id, level->linkage);
//    for (target = level->link->target; target->prev; target = target->prev)  printf(" - %2d ", target->index_c);
  }
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
hc_average_linkage (hc_link_t *source, hc_link_t *target) 
{ // target is the one with new values (see += below...)
  return ((source->distance * source->count) + (target->distance * target->count)) / (target->count += source->count);
}

hc_link_t* 
hc_merge_link (hc_link_t *source, hc_link_t *target, hierarchical_cluster hac) 
{ 
  switch (hac->linkage) {
    case ('u'):
    default:
      target->distance = hc_average_linkage (source, target);
      break;
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
      // lastLink->clusters[index++] = cluster->index;
      // lastLink->clusters[index++] = cluster->length;
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
  target->length += source->length;
  while (target->cluster != NULL) target = target->cluster;
  target->cluster = source;
  hc_unlink (source);
}

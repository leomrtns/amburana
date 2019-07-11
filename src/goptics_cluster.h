/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file goptics_cluster.h 
 *  \brief OPTICS algorithm based on https://github.com/guineri/GOPTICS  
 */
#ifndef _amburana_goptics_cluster_h_
#define _amburana_goptics_cluster_h_

#include "sketch_distance.h" 

typedef struct goptics_cluster_struct* goptics_cluster;

struct goptics_cluster_struct
{
  int *Va_i, *Va_n; // Va_i[pts] where starts at Ea_ids and Ea_dist list; va_n[pts] = number of neighbours <both opaque>
  double epsilon; 
  int min_points, num_edges, n_clusters; // minpts from user, num_edges = number of dists < epsilon 
  int *cluster, *order, n_order; // samples sorted by reachability order, and resulting cluster number
  double *core_distance, *reach_distance, max_distance; // max_dist is a convenience number to replace DBL_MAX in output
  bool *core;
  void *Ea, *heap, *points; // void b/c I don't want to expose local structs
  distance_generator d; // d->n_samples 
};

goptics_cluster new_goptics_cluster (distance_generator dg, int min_points, double epsilon);
goptics_cluster new_goptics_cluster_run (distance_generator dg, int min_points, double epsilon);
void del_goptics_cluster (goptics_cluster gop);
void assign_goptics_clusters (goptics_cluster gop, double cluster_eps);

#endif

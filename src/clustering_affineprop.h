/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file clustering_affineprop.h 
 *  \brief Affinity propagation clustering  
 */
#ifndef _clustering_affineprop_h_
#define _clustering_affineprop_h_
#include "sketch_distance.h" 

typedef struct affineprop_cluster_struct* affineprop_cluster;

struct affineprop_cluster_struct
{
  double *A, *R, *S;
  double preference;
  int *cluster, n_clusters;
  distance_generator d;
};

affineprop_cluster new_affineprop_cluster (distance_generator dg);
void del_affineprop_cluster (affineprop_cluster ap);
void affineprop_run (affineprop_cluster ap, int iterations, double quantile);

#endif

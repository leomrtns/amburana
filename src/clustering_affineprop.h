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
  double *A, *R, *S, *tmp;
  double *pref;
  int *cluster, *exemplars, n_clusters; /*!< \brief cluster[n_samples], exemplars[n_clusters], with exemplars[ cluster[i] ] giving exemplar for sample i */
  int n_converging;
  distance_generator d;
  rng_tt800_struct tt800;
};

affineprop_cluster new_affineprop_cluster (distance_generator dg);
void del_affineprop_cluster (affineprop_cluster ap);
void affineprop_run (affineprop_cluster ap, int iterations, double quantile, double lambda);

#endif

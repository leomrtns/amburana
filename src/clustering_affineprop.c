/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/* Affinity Propagation, based on implementation by Sunghoon Heo
 * Reference : Brandan J Frey Science 2007 Feb  https://github.com/tahuh/ap_clust/
 */
#include "clustering_affineprop.h"

static void 
affineprop_update_iteration (affineprop_cluster ap);

affineprop_cluster
new_affineprop_cluster (distance_generator dg)
{
  affineprop_cluster ap = (affineprop_cluster) biomcmc_malloc (sizeof (struct affineprop_cluster_struct));
  ap->d = dg; dg->ref_counter++;
  
  ap->A = (double*) biomcmc_malloc (dg->n_samples * dg->n_samples * sizeof (double));
  ap->R = (double*) biomcmc_malloc (dg->n_samples * dg->n_samples * sizeof (double));
  ap->S = (double*) biomcmc_malloc (dg->n_samples * dg->n_samples * sizeof (double)); // use it as *(A + i*n + j)
  ap->cluster = (int*) biomcmc_malloc (dg->n_samples * sizeof (int));  
  ap->n_clusters = 0;
  ap->preference = 1.; // s[i][i], close to min d(i,j) is fewer classes, closer to max d[i,j] is more classes
  return ap;
}

void 
del_affineprop_cluster (affineprop_cluster ap)
{
  if (!ap) return;
  if (ap->A) free (ap->A);
  if (ap->R) free (ap->R);
  if (ap->S) free (ap->S);
  if (ap->cluster) free (ap->cluster);
  del_distance_generator (ap->d);
  free (ap);
}

void
affineprop_run (affineprop_cluster ap, int iterations, double quantile)
{ // quantile is value for sim[i][i], or input preference, usually median (=0.5)
  int i, j, loc, ns = ap->d->n_samples, n = 0;
  double minvec, *vec = (double*) biomcmc_malloc (ap->d->n_samples * sizeof (double));
  if (quantile < 0.01) quantile = 0.01;
  if (quantile > 0.99) quantile = 0.99;

  /* 1. Find preference based on similarities */  // A[] is used temporarily to store similarities (-1 x distance)
  for(j = 1; j < ns; j++) for (i = 0; i < j; i++) ap->A[n++] = - distance_generator_get (ap->d, i, j); // negative distance = similarity
  ap->preference = biomcmc_wirth_algorithm (ap->A, ap->d->n_samples, (int)(quantile * (double)(ap->d->n_samples)));
  /* 2. Fill similarity matrix, including S[i][i] preference */
  n = 0; // if we subsample, then S[i][i] can take into account previous iterations (i.e. previous exemplars have higher value)
  for (j = 1; j < ns; j++) for (i = 0; i < j; i++) *(ap->S + i * ns + j) = *(ap->S + j * ns + i) =  ap->A[n++];
  for (j = 0; j < ns; j++) *(ap->S + j * ns + j) = ap->preference; // how likely sample i is of becoming an exemplar
  for (j = 0; j < ns*ns; j++) ap->A[i] = ap->R[i] = 0.;
  /* 3. update Responsibility and Availability matrices */
  for (i = 0; i < iterations; i++) affineprop_update_iteration (ap);
  /* 4. Find exemplars */
  for (n = 0; n < ns; n++) {
    for (j = 0; j < ns; j++) vec[j] = *(ap->A + n * ns + j) + *(ap->R + n * ns + j); // A[i][k] + R[i][k];
    minvec = DBL_MAX; loc = 0; // original call max but stores large negative values
    for (j = 0; j < ns; j++) if (vec[j] < minvec) {minvec = vec[j]; loc = j; } 
    ap->cluster[n] = loc;
  }
  free (vec);
}

static void 
affineprop_update_iteration (affineprop_cluster ap)
{
  double sikp, a_ik;
  double max_rik = -1e-32, sum = 0.0;
  unsigned int i, ip, k, kp, ns = ap->d->n_samples;
  for (i = 0; i < ns; i++) for (k = 0; k < ns; k++) {
    for (kp = 0; kp < ns; kp++) if (kp != k) {
      sikp = *(ap->S + i * ns + kp) + *(ap->A + i * ns + kp);
      if (sikp > max_rik) max_rik = sikp;
    }
    *(ap->R + i * ns + k) =  *(ap->S + i * ns + k) - max_rik; // R[i][k] = S[i][k] - max_rik;
    a_ik = *(ap->A + i * ns + k);
    if (k == i) {
      for (ip = 0; ip < ns; ip++) if (ip != i) a_ik += (0 > *(ap->R + ip * ns + k)? 0 : *(ap->R + ip * ns + k));
      *(ap->A + i * ns + k) = a_ik;
      continue;	// force loop
    } else {
      sum = 0.0;
      for (ip = 0; ip < ns; ip++) if ((ip != k) && (ip != i)) sum += (0 > *(ap->R + ip * ns + k)? 0 : *(ap->R + ip * ns + k));
      a_ik = *(ap->R + k * ns + k) + sum;
      if (a_ik > 0.) a_ik = 0.;
      *(ap->A + i * ns + k) = a_ik;
      continue;	// force loop
    }
  }
}


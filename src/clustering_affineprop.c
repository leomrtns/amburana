/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */
/* Affinity Propagation, based on  https://github.com/tahuh/ap_clust/ and https://github.com/UBod/apcluster  */

#include "clustering_affineprop.h"
#define ITER_CONVERGE 20

static void affineprop_generate_similarity_matrix (affineprop_cluster ap, double quantile);
static bool affineprop_update_iteration (affineprop_cluster ap, double lambda);

affineprop_cluster
new_affineprop_cluster (distance_generator dg)
{
  int i;
  int64_t time_seed[2];
  affineprop_cluster ap = (affineprop_cluster) biomcmc_malloc (sizeof (struct affineprop_cluster_struct));
  ap->d = dg; dg->ref_counter++;
  
  ap->A = (double**) biomcmc_malloc (dg->n_samples * sizeof (double*));
  ap->R = (double**) biomcmc_malloc (dg->n_samples * sizeof (double*));
  ap->S = (double**) biomcmc_malloc (dg->n_samples * sizeof (double*));
  for (i = 0; i < dg->n_samples; i++) {
    ap->A[i] = (double*) biomcmc_malloc (dg->n_samples * sizeof (double));
    ap->R[i] = (double*) biomcmc_malloc (dg->n_samples * sizeof (double));
    ap->S[i] = (double*) biomcmc_malloc (dg->n_samples * sizeof (double));
  }
  ap->tmp  = (double*) biomcmc_malloc (dg->n_samples * sizeof (double)); // vector[n], not matrix[n][n] as above
  ap->pref = (double*) biomcmc_malloc (dg->n_samples * sizeof (double)); 
  ap->cluster = (int*) biomcmc_malloc (dg->n_samples * sizeof (int));  
  ap->exemplars = (int*) biomcmc_malloc (2 * dg->n_samples * sizeof (int)); // shoud b smaller than n_samples, but we also use for convergence 
  ap->n_clusters = 0;
  ap->n_converging = 0;
  ap->timing_secs = 0.;

  biomcmc_get_time (time_seed);
  rng_set_tt800 (&(ap->tt800), time_seed[0] | time_seed[1]);
  return ap;
}

void 
del_affineprop_cluster (affineprop_cluster ap)
{
  int i;
  if (!ap) return;
  if (ap->A) { for (i = ap->d->n_samples-1; i >= 0; i--) if (ap->A[i]) free (ap->A[i]);
    free (ap->A);
  }
  if (ap->R) { for (i = ap->d->n_samples-1; i >= 0; i--) if (ap->R[i]) free (ap->R[i]);
    free (ap->R);
  }
  if (ap->S) { for (i = ap->d->n_samples-1; i >= 0; i--) if (ap->S[i]) free (ap->S[i]);
    free (ap->S);
  }
  if (ap->tmp) free (ap->tmp);
  if (ap->cluster) free (ap->cluster);
  if (ap->exemplars) free (ap->exemplars);
  del_distance_generator (ap->d);
  free (ap);
}

void
affineprop_run (affineprop_cluster ap, int iterations, double quantile, double lambda)
{ // quantile is value for sim[i][i], or input preference, usually median (=0.5)
  int i, j, ns = ap->d->n_samples;
  double max_sim, *vec = (double*) biomcmc_malloc (ap->d->n_samples * sizeof (double));
  bool converged = false;
  clock_t time1, time0;
  if (quantile < 0.01) quantile = 0.01;
  if (quantile > 0.99) quantile = 0.99;
  if (lambda > 0.99) lambda = 0.99;
  if (lambda < 0.01) lambda = 0.01;
  time0 = clock ();

  /* 1. similarity (negative distance) and preferences similarity[i][i] into S[] */
  affineprop_generate_similarity_matrix (ap, quantile);

  for (i = 0; i < ns; i++) ap->cluster[i] = ap->exemplars[i] = ap->exemplars[ns + i] = -1; // exemplars[ 2 x n_samples ]

  /* 2. update Responsibility and Availability matrices R[] and A[] */
  for (i = 0, converged = false; (i < iterations) && (!converged); i++) converged = affineprop_update_iteration (ap, lambda);
  if (converged) ap->n_converging = i; // recycling n_converging, which was used to count how many iterations with stable exemplars 
  else ap->n_converging = 0; 

  /* 3. Fill cluster entries for exemplars (which are calculated at every iteration for convergence) */
  for (i = 0; i < ap->n_clusters; i++) ap->cluster[ ap->exemplars[i] ] = i; // cluster = index in n_samples of closest exemplar

  /* 4. Find exemplar closest to sample */ 
  for (i = 0; i < ns; i++) if (ap->cluster[i] < 0) { // only for non-exemplars
    max_sim = ap->S[i][ap->exemplars[0]]; 
    ap->cluster[i] = 0;
    for (j = 1; j < ap->n_clusters; j++) if (max_sim < ap->S[i][ap->exemplars[j]]) {
      max_sim = ap->S[i][ap->exemplars[j]]; 
      ap->cluster[i] = j;
    }
  }
  time1 = clock (); ap->timing_secs += (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  free (vec);
}

static void
affineprop_generate_similarity_matrix (affineprop_cluster ap, double quantile)
{
  int i, j, k, ns = ap->d->n_samples;
  double min_sim = DBL_MAX, noise;

  /* 1. Find preference based on similarities */  
  for(i = 0; i < ns; i++) {
    k = 0;
    for (j = 0; j < ns; j++) if (i != j) ap->tmp[k++] = - distance_generator_get (ap->d, i, j); // negative distance = similarity
    ap->S[i][i] = biomcmc_wirth_algorithm (ap->tmp, k, (int)(quantile * (double)(k)));
  }

  /* 3. Fill similarity matrix, with jitter (noise) to remove degeneracies */
  for (i = 0; i < ns; i++) {
    for (j = 0; j < i; j++) { 
      noise = 1.e-12 * ((double)(rng_get_tt800 (&(ap->tt800)))/(double)(0xffffffffU)); // unif(0,1) 
      ap->S[i][j] = - distance_generator_get (ap->d, i, j) + noise;
      if (ap->S[i][j] < min_sim) min_sim = ap->S[i][j];
    }
    noise = 1.e-12 * ((double)(rng_get_tt800 (&(ap->tt800)))/(double)(0xffffffffU)); // unif(0,1) 
    ap->S[i][i] += noise;
    if (ap->S[i][i] < min_sim) min_sim = ap->S[i][i];
  }

  /* 4. Rescale similarity to be positive (notice the j<=i , to include S[i][i])  */
  for (i = 0; i < ns; i++) for (j = 0; j <=i; j++) { ap->S[i][j] -= min_sim; ap->S[j][i] = ap->S[i][j]; }

  /* 5. Clean A[] (and also R[] just to be safe) */
  for (i = 0; i < ns; i++) for (j = 0; j <=i; j++)  ap->A[j][i] = ap->R[j][i] = ap->A[i][j] = ap->R[i][j] = 0.;
}

static bool 
affineprop_update_iteration (affineprop_cluster ap, double lambda)
{
  double aux_dbl, max1 = -1e-62, max2 = -1e-62, new_val;
  int i, j, k, ns = ap->d->n_samples;

  for (i = 0; i < ns; i++) { // responsibilities
    max1 = max2 = -1e-62; k = 0;
    for (j = 0; j < ns; j++) {
      aux_dbl=  ap->S[i][j] + ap->A[i][j];
      if (aux_dbl > max1) { max2 = max1; max1 = aux_dbl; k = j; }
      else if (aux_dbl > max2) max2 = aux_dbl;
    }
    for (j = 0; j < ns; j++) {
      new_val = ap->S[i][j] - (j == k ? max2 : max1);
      ap->R[i][j] = (1 - lambda) * new_val + lambda *  ap->R[i][j];
    }
  }

  for (i = 0; i < ns; i++) { // availabilities 
    aux_dbl = 0.;
    for (j = 0; j < ns; j++) {
      if ((ap->R[j][i] < 0.) && (j != i)) ap->tmp[j] = 0.;
      else ap->tmp[j] = ap->R[j][i];
      aux_dbl += ap->tmp[j];
    }
    for (j = 0; j < ns; j++) {
      new_val = aux_dbl - ap->tmp[j];
      if ((new_val > 0) && (j != i)) new_val = 0.;
      ap->A[j][i] = (1 - lambda) * new_val + lambda * ap->A[j][i];
    }
  }

  j = ns;  // check for convergence, using exemplars[ns + i] 
  for (i = 0; i < ns; i++) if ((ap->R[i][i] + ap->A[i][i]) > 0.) ap->exemplars[j++] = i;
  j -= ns; // now n has new number of clusters
  if ((j > 0) && (j == ap->n_clusters)) for (i = 0; (i < j) && (ap->exemplars[i] == ap->exemplars[ns + i]); i++);
  if (i == j) { 
    ap->n_converging++; // all exemplars are same as previous iteration
    if (ap->n_converging == ITER_CONVERGE) return true;
    return false;
  }
  ap->n_clusters = j; // not the same, so we must update:
  for (i = 0; i < j; i++) ap->exemplars[i] = ap->exemplars[ns+i];
  ap->n_converging = 0;
  return false;
}

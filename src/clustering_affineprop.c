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
  int time_seed[2];
  affineprop_cluster ap = (affineprop_cluster) biomcmc_malloc (sizeof (struct affineprop_cluster_struct));
  ap->d = dg; dg->ref_counter++;
  
  ap->A = (double*) biomcmc_malloc (dg->n_samples * dg->n_samples * sizeof (double));
  ap->R = (double*) biomcmc_malloc (dg->n_samples * dg->n_samples * sizeof (double));
  ap->S = (double*) biomcmc_malloc (dg->n_samples * dg->n_samples * sizeof (double)); // use it as *(A + i*n + j)
  ap->tmp  = (double*) biomcmc_malloc (dg->n_samples * sizeof (double)); // vector[n], not matrix[n][n] as above
  ap->pref = (double*) biomcmc_malloc (dg->n_samples * sizeof (double)); 
  ap->cluster = (int*) biomcmc_malloc (dg->n_samples * sizeof (int));  
  ap->exemplars = (int*) biomcmc_malloc (2 * dg->n_samples * sizeof (int)); // shoud b smaller than n_samples, but we also use for convergence 
  ap->n_clusters = 0;
  ap->n_converging = 0;

  biomcmc_get_time (time_seed);
  rng_set_tt800 (&(ap->tt800), time_seed[0] | time_seed[1]);
  return ap;
}

void 
del_affineprop_cluster (affineprop_cluster ap)
{
  if (!ap) return;
  if (ap->A) free (ap->A);
  if (ap->R) free (ap->R);
  if (ap->S) free (ap->S);
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
  if (quantile < 0.01) quantile = 0.01;
  if (quantile > 0.99) quantile = 0.99;
  if (lambda > 0.99) lambda = 0.99;
  if (lambda < 0.01) lambda = 0.01;

  affineprop_generate_similarity_matrix (ap, quantile);

  for (i = 0; i < ns; i++) ap->cluster[i] = ap->exemplars[i] = ap->exemplars[ns + i] = -1; // exemplars[ 2 x n_samples ]
//  for (j = 0; j < 10; j++) { for (i = 0; i < 10; i++) printf ("%7.5e ", *(ap->S + j * ns + i)); printf ("   DEBUG\n"); }

  /* 3. update Responsibility and Availability matrices */
  for (i = 0, converged = false; (i < iterations) && (!converged); i++) converged = affineprop_update_iteration (ap, lambda);
  if (converged) ap->n_converging = i; // recycling n_converging, which was used to count how many iterations with stable exemplars 
  else ap->n_converging = 0; 

  /* 4. Fill cluster entries for exemplars (which are calculated at every iteration for convergence) */
  for (i = 0; i < ap->n_clusters; i++) ap->cluster[ ap->exemplars[i] ] = i; // cluster = index in n_samples of closest exemplar

  /* 5. Find exemplar closest to sample */ 
  for (i = 0; i < ns; i++) if (ap->cluster[i] < 0) { // only for non-exemplars
    max_sim = *(ap->S + i * ns + ap->exemplars[0]); 
    ap->cluster[i] = 0;
    for (j = 1; j < ap->n_clusters; j++) if (max_sim < *(ap->S + i * ns + ap->exemplars[j])) {
      max_sim = *(ap->S + i * ns + ap->exemplars[j]); 
      ap->cluster[i] = j;
    }
  }
  free (vec);
}

static void
affineprop_generate_similarity_matrix (affineprop_cluster ap, double quantile)
{
  int i, j, ns = ap->d->n_samples;
  double min_sim;

  /* Find preference based on similarities */  // A[] is used temporarily to store similarities (-1 x distance)
  for(j = 0; j < ns; j++) { 
    for (i = 0; i < ns; i++) ap->A[i] = - distance_generator_get (ap->d, i, j); // negative distance = similarity
    ap->pref[j] = biomcmc_wirth_algorithm (ap->A, ns, (int)(quantile * (double)(ns)));
  }
  /* Fill similarity matrix, including S[i][i] preference, _temporarily_ into A[] (since wirth algo _modifies_ vector) */
  for (j = 1; j < ns; j++) for (i = 0; i < j; i++) 
    *(ap->A + i * ns + j) = *(ap->A + j * ns + i) = - distance_generator_get (ap->d, i, j);
  for (j = 0; j < ns; j++) *(ap->A + j * ns + j) = ap->pref[j]; 
  /* Jitter similarity to remove degeneracies */
  for (j = 0; j < ns*ns; j++) ap->A[j] += 1.e-12 * ((double)(rng_get_tt800 (&(ap->tt800)))/(double)(0xffffffffU)); // unif(0,1) 
  /* Rescale similarity to be positive */
  for (j = 0; j < ns*ns; j++) ap->S[j] = ap->A[j]; // first we copy A back to S
  min_sim = biomcmc_wirth_algorithm (ap->A, ns * ns, 0); // This function modifies ap->A . min_value is the most negative ;)
  for (j = 0; j < ns*ns; j++) ap->S[j] -= min_sim;
  /* clean A[] (and also R[] just to be safe) */
  for (j = 0; j < ns*ns; j++) ap->A[i] = ap->R[i] = 0.;
}

static bool 
affineprop_update_iteration (affineprop_cluster ap, double lambda)
{
  double aux_dbl, max1 = -1e-62, max2 = -1e-62, new_val;
  int i, j, k, ns = ap->d->n_samples;

  for (i = 0; i < ns; i++) { // responsibilities
    max1 = max2 = -1e-62; k = 0;
    for (j = 0; j < ns; j++) {
      aux_dbl= *(ap->S + i * ns + j) + *(ap->A + i * ns + j);
      if (aux_dbl > max1) {
        max2 = max1;
        max1 = aux_dbl;
        k = j;
      }
      else if (aux_dbl > max2) max2 = aux_dbl;
    }
    for (j = 0; j < ns; j++) {
      new_val = *(ap->S + i * ns + j) - (j == k ? max2 : max1);
      *(ap->R + i * ns + j) = (1 - lambda) * new_val + lambda *   *(ap->R + i * ns + j);
    }
  }

  for (i = 0; i < ns; i++) { // availabilities 
    aux_dbl = 0.;
    for (j = 0; j < ns; j++) {
      if ((*(ap->R + j * ns + i) < 0.) && (j != i)) ap->tmp[j] = 0.;
      else ap->tmp[j] = *(ap->R + j * ns + i);
      aux_dbl += ap->tmp[j];
    }
    for (j = 0; j < ns; j++) {
      new_val = aux_dbl - ap->tmp[j];
      if ((new_val > 0) && (j != i)) new_val = 0.;
      *(ap->A + j * ns + i) = (1 - lambda) * new_val + lambda *   *(ap->A + j * ns + i);
    }
  }

  j = ns;  // check for convergence, using exemplars[ns + i] 
  for (i = 0; i < ns; i++) if ((*(ap->R + i * ns + i) + *(ap->A + i * ns + i)) > 0.) ap->exemplars[j++] = i;
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

void 
affineprop_update_iteration_old (affineprop_cluster ap)
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
  }
  for (i = 0; i < ns; i++) for (k = 0; k < ns; k++) {
    a_ik = *(ap->A + i * ns + k); // maybe it's wrong?
    a_ik = 0.;
    if (k == i) {
      for (ip = 0; ip < ns; ip++) if (ip != i) a_ik += (0 > *(ap->R + ip * ns + k)? 0 : *(ap->R + ip * ns + k));
      *(ap->A + i * ns + k) = a_ik;
      continue;	// force loop
    } else {
      sum = 0.0;
      for (ip = 0; ip < ns; ip++) if ((ip != k) && (ip != i)) sum += (0 > *(ap->R + ip * ns + k)? 0 : *(ap->R + ip * ns + k));
      a_ik = *(ap->R + k * ns + k) + sum;
      if (a_ik > 0.) a_ik = 0.; // a[i,k] = MIN(0, a_ik)  [[ above it's MAX() ]]
      *(ap->A + i * ns + k) = a_ik;
      continue;	// force loop
    }
  }
}

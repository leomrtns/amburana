/* 
 * This file is part of amburana. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file clustering_hierarchical.h 
 *  \brief Hierarchical agglomerative clustering, based on https://github.com/rchipka/libhcluste 
 */
#ifndef _clustering_hierarchical_h_
#define _clustering_hierarchical_h_
#include "sketch_distance.h" 

topology hierarchical_cluster_topology (distance_generator dg, char linkage);

#endif

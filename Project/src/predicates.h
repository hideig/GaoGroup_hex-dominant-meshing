/*
 * HXT - Copyright (C) <2016-2018> <Université catholique de Louvain (UCL), Belgique>
 *
 * List of the contributors to the development of HXT: see AUTHORS file.
 * Description and complete License: see LICENSE file.
 * 	
 * This program (HXT) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU Lesser
 * General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING and COPYING.LESSER files).  If not, 
 * see <http://www.gnu.org/licenses/>.
 */

#ifndef _ROBUST_PREDICATES_H_
#define _ROBUST_PREDICATES_H_

#include "hxt_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

extern double splitter;
extern double o3dstaticfilter;
extern double o3derrboundA;
extern double ispstaticfilter;
extern double isperrboundA;

double exactinit(double maxx, double maxy, double maxz);
// double incircle(double *pa, double *pb, double *pc, double *pd);

/** \todo Please comment the freaking variable type you are using. [JP]
 */
double insphere(
  const double* const __restrict__ pa,
  const double* const __restrict__ pb,
  const double* const __restrict__ pc,
  const double* const __restrict__ pd,
  const double* const __restrict__ pe);

double orient3dd(
  const double* const __restrict__ pa,
  const double* const __restrict__ pb,
  const double* const __restrict__ pc,
  const double* const __restrict__ pd);

int grow_expansion(
  int elen, const double *e, double b, double *h);

int grow_expansion_zeroelim(
  int elen, const double *e, double b, double *h);

int fast_expansion_sum_zeroelim(
  const int elen, const double *e,
  const int flen, const double *f, double *h);

int fast_expansion_sum(
  int elen, const double *e, int flen, const double *f, double *h);

int fast_expansion_sum_zeroelim(
  const int elen, const double *e,
  const int flen, const double *f, double *h);

int scale_expansion(
  const int elen, const double *e, const double b, double *h);

int scale_expansion_zeroelim(
  const int elen, const double *e, const double b, double *h);

#ifdef __cplusplus
}
#endif

#endif

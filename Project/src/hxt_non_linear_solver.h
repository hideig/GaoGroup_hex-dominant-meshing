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

#ifndef HXT_NON_LINEAR_SOLVER_H
#define HXT_NON_LINEAR_SOLVER_H

#include "hxt_linear_system.h"
#include "hxt_api.h"

typedef HXTStatus HXTNonLinearSolverCallbackF(void *data, const double *solution, HXTLinearSystem *sys, double *f);
typedef HXTStatus HXTNonLinearSolverCallbackDF(void *data, const double *solution, HXTLinearSystem *sys);

HXTStatus hxtNewtonRaphson(HXTLinearSystem *nrSys, double *solution, int size, int maxiter, double tol, HXTNonLinearSolverCallbackF *fcb, HXTNonLinearSolverCallbackDF *dfcb, void *data);


#ifdef HXT_HAVE_PETSC
HXTStatus hxtNonLinearSolverPETSc(HXTLinearSystem *sys, double *solution, HXTNonLinearSolverCallbackF *f, HXTNonLinearSolverCallbackDF *df, void *data, const char *petscOptions);
typedef HXTStatus HXTOptimizationTaoCallbackF(void *data, const double *solution, double *f);
typedef HXTStatus HXTOptimizationTaoCallbackG(void *data, const double *solution, HXTLinearSystem *sys, double *f);
typedef HXTStatus HXTOptimizationTaoCallbackH(void *data, const double *solution, HXTLinearSystem *sys);
HXTStatus hxtOptimizationTao(HXTLinearSystem *sys, double *solution, HXTOptimizationTaoCallbackF *fcb, HXTOptimizationTaoCallbackG *gcb, HXTOptimizationTaoCallbackH *hcb, void *data, const char *petscOptions);
#endif

#endif

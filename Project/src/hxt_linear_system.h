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

#ifndef HEXTREME_LINEAR_SYSTEM_H
#define HEXTREME_LINEAR_SYSTEM_H
#include <stdint.h>
#include "hxt_api.h"

typedef struct HXTLinearSystemStruct HXTLinearSystem;

HXTStatus hxtLinearSystemAddToMatrix(HXTLinearSystem *lsys, int el0, int el1, const double *localMatrix);
HXTStatus hxtLinearSystemAddMatrixEntry(HXTLinearSystem *lsys, int node0, int field0, int node1, int field1, double entry);
HXTStatus hxtLinearSystemAddToRhs(HXTLinearSystem *lsys, double *rhs, int el0, const double *localVector);
HXTStatus hxtLinearSystemSize(HXTLinearSystem *lsys, int *size);
HXTStatus hxtLinearSystemZeroMatrix(HXTLinearSystem *lsys);
HXTStatus hxtLinearSystemSolve(HXTLinearSystem *lsys, double *rhs, double *solution);
HXTStatus hxtLinearSystemSetMatrixRowIdentity(HXTLinearSystem *lsys, int node, int field);
HXTStatus hxtLinearSystemSetMatrixRowFieldCombinaison(HXTLinearSystem *system, int node, int field, double *coeff);
HXTStatus hxtLinearSystemSetRhsEntry(HXTLinearSystem *lsys, double *rhs, int node, int field, double v);
HXTStatus hxtLinearSystemAddRhsEntry(HXTLinearSystem *lsys, double *rhs, int node, int field, double v);
HXTStatus hxtLinearSystemDelete(HXTLinearSystem **pSystem);
HXTStatus hxtLinearSystemGetRhsNorm(HXTLinearSystem *lsys, double *rhs, double *norm);
HXTStatus hxtLinearSystemHasConverged(HXTLinearSystem *lsys, int *converged);

HXTStatus hxtLinearSystemCreateLU(HXTLinearSystem **sys, int nElement, int nNodesByElement, int nFields, uint32_t *elements);
typedef struct HXTLinearSystemLUStruct HXTLinearSystemLU;
HXTStatus hxtLinearSystemGetLinearSystemLU(HXTLinearSystem *sys, HXTLinearSystemLU **psys);

#ifdef HXT_HAVE_PETSC
typedef struct HXTLinearSystemPETScStruct HXTLinearSystemPETSc;
HXTStatus hxtLinearSystemCreatePETSc(HXTLinearSystem **sys, int nElement, int nNodesByElement, int nFields, uint32_t *elements, const char *petscOptions);
HXTStatus hxtLinearSystemGetLinearSystemPETSc(HXTLinearSystem *sys, HXTLinearSystemPETSc **psys);
#endif

HXTStatus hxtInitializeLinearSystems(int *argc, char ***argv);

#endif

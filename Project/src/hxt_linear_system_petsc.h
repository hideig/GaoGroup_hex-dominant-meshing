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

#ifndef HEXTREME_LINEAR_SYSTEM_PETSC_H
#define HEXTREME_LINEAR_SYSTEM_PETSC_H
#include <stdint.h>
#include "hxt_api.h"
#include <petscmat.h>

typedef struct HXTLinearSystemPETScStruct HXTLinearSystemPETSc;

HXTStatus hxtLinearSystemPETScAddToMatrix(HXTLinearSystemPETSc *lsys, int el0, int el1, const double *localMatrix);
HXTStatus hxtLinearSystemPETScAddMatrixEntry(HXTLinearSystemPETSc *lsys, int node0, int field0, int node1, int field1, double entry);
HXTStatus hxtLinearSystemPETScAddToRhs(HXTLinearSystemPETSc *lsys, double *rhs, int el0, const double *localVector);
HXTStatus hxtLinearSystemPETScZeroMatrix(HXTLinearSystemPETSc *lsys);
HXTStatus hxtLinearSystemPETScSolve(HXTLinearSystemPETSc *lsys, double *rhs, double *solution);
HXTStatus hxtLinearSystemPETScSetMatrixRowIdentity(HXTLinearSystemPETSc *lsys, int node, int field);
HXTStatus hxtLinearSystemPETScSetMatrixRowFieldCombinaison(HXTLinearSystemPETSc *system, int node, int field, double *coeff);
HXTStatus hxtLinearSystemPETScSetRhsEntry(HXTLinearSystemPETSc *lsys, double *rhs, int node, int field, double v);
HXTStatus hxtLinearSystemPETScAddRhsEntry(HXTLinearSystemPETSc *lsys, double *rhs, int node, int field, double v);
HXTStatus hxtLinearSystemPETScDelete(HXTLinearSystemPETSc **pSystem);
HXTStatus hxtLinearSystemPETScGetRhsNorm(HXTLinearSystemPETSc *lsys, double *rhs, double *norm);
HXTStatus hxtLinearSystemPETScCreate(HXTLinearSystemPETSc **pSystem, int nElements, int nNodesByElement, int nFields, uint32_t *elements, const char *petscOptions);
HXTStatus hxtLinearSystemPETScHasConverged(HXTLinearSystemPETSc *lsys, int *converged);
HXTStatus hxtLinearSystemPETScSize(HXTLinearSystemPETSc *lsys, int *size);

HXTStatus hxtLinearSystemPETScGetMat(HXTLinearSystemPETSc *lsys, Mat *mat);
HXTStatus hxtLinearSystemPETScMapToVec(HXTLinearSystemPETSc *lsys, double *v, Vec vec);
HXTStatus hxtLinearSystemPETScMapFromVec(HXTLinearSystemPETSc *lsys, Vec vec, double *v);

HXTStatus hxtInitializePETSc(int *argc, char ***argv);

#endif

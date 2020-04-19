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

#ifndef HEXTREME_LINEAR_SYSTEM_LU_H
#define HEXTREME_LINEAR_SYSTEM_LU_H
#include <stdint.h>
#include "hxt_api.h"

typedef struct HXTLinearSystemLUStruct HXTLinearSystemLU;

HXTStatus hxtLinearSystemLUAddToMatrix(HXTLinearSystemLU *lsys, int el0, int el1, const double *localMatrix);
HXTStatus hxtLinearSystemLUAddMatrixEntry(HXTLinearSystemLU *lsys, int node0, int field0, int node1, int field1, double entry);
HXTStatus hxtLinearSystemLUAddToRhs(HXTLinearSystemLU *lsys, double *rhs, int el0, const double *localVector);
HXTStatus hxtLinearSystemLUZeroMatrix(HXTLinearSystemLU *lsys);
HXTStatus hxtLinearSystemLUSolve(HXTLinearSystemLU *lsys, double *rhs, double *solution);
HXTStatus hxtLinearSystemLUSetMatrixRowIdentity(HXTLinearSystemLU *lsys, int node, int field);
HXTStatus hxtLinearSystemLUSetMatrixRowFieldCombinaison(HXTLinearSystemLU *system, int node, int field, double *coeff);
HXTStatus hxtLinearSystemLUSetRhsEntry(HXTLinearSystemLU *lsys, double *rhs, int node, int field, double v);
HXTStatus hxtLinearSystemLUAddRhsEntry(HXTLinearSystemLU *lsys, double *rhs, int node, int field, double v);
HXTStatus hxtLinearSystemLUDelete(HXTLinearSystemLU **pSystem);
HXTStatus hxtLinearSystemLUHasConverged(HXTLinearSystemLU *lsys, int *converged); 
HXTStatus hxtLinearSystemLUGetRhsNorm(HXTLinearSystemLU *lsys, double *rhs, double *norm);
HXTStatus hxtLinearSystemLUCreate(HXTLinearSystemLU **pSystem, int nElements, int nNodesByElement, int nFields, uint32_t *elements);
HXTStatus hxtLinearSystemLUSize(HXTLinearSystemLU *lsys, int *size);

#endif

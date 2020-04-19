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

#ifndef        _HXT_BBOX_H_
#define _HXT_BBOX_H_

#include <float.h>
#include "hxt_api.h"
#include "hxt_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct hxtBboxStruct{
    double hxtDeclareAligned min[3];
    double hxtDeclareAligned max[3];
} hxtBbox;

// /* create a new bounding box */
// HXTStatus hxtBboxCreate(hxtBbox** bboxP);

// /* delete a bounding box */
// HXTStatus hxtBboxDelete(hxtBbox** bboxP);

static inline void hxtBboxInit(hxtBbox* bbox){
  bbox->min[0] = DBL_MAX;
  bbox->min[1] = DBL_MAX;
  bbox->min[2] = DBL_MAX;
  bbox->max[0] = -DBL_MAX;
  bbox->max[1] = -DBL_MAX;
  bbox->max[2] = -DBL_MAX;
}

/* update the bounding box with one new vertex */
HXTStatus hxtBboxAddOne(hxtBbox* bbox, double* coord);

/* update the bounding box with an array of n vertices at once (far quicker) */
HXTStatus hxtBboxAdd(hxtBbox* bbox, double* coord, uint32_t n);

#ifdef __cplusplus
}
#endif

#endif

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

#ifndef HEXTREME_EDGE_H
#define HEXTREME_EDGE_H

#include "hxt_tools.h"
#include "hxt_mesh.h"
#include <math.h>

typedef struct hxtLineLoopStruct HXTLineLoop;

typedef struct hxtBoundariesStruct HXTBoundaries;

typedef struct hxtEdgesStruct{

  HXTMesh *edg2mesh;
  uint64_t *global; // initial number of elements
  
  uint32_t numEdges;// number of edges

  uint32_t* node;// edg2vertices

  uint16_t* color;// to be used ...

  uint64_t* edg2tri;
  uint32_t* tri2edg;
  uint32_t *bdryedges;// physical lines/boundary
  uint64_t nEdgesBdry;  
  int edgesBdryTotal;// to be removed

} HXTEdges; 


double hxtEdgesLength(const HXTEdges *edges,uint32_t ie);

HXTStatus hxtLineLoopGetEdges(HXTLineLoop *loop, uint32_t **edges);
HXTStatus hxtLineLoopGetNumberOfEdges(HXTLineLoop *loop, int *number);
HXTStatus hxtLineLoopGetLength(HXTLineLoop *loops, double *length);

HXTStatus hxtBoundariesGetNumberOfLineLoops(HXTBoundaries *boundaries, int *number);
HXTStatus hxtBoundariesGetNumberOfBorderEdges(HXTBoundaries *boundaries, int *number);
HXTStatus hxtBoundariesGetNumberOfEdgesOfLineLoop(HXTBoundaries *boundaries, int i, int *number);
HXTStatus hxtBoundariesGetEdgesOfLineLoop(HXTBoundaries *boundaries, int i, uint32_t **edges);
HXTStatus hxtBoundariesGetLengthOfLineLoop(HXTBoundaries *boundaries, int i, double *length);
HXTStatus hxtBoundariesGetLineLoop(HXTBoundaries *boundaries, int i,HXTLineLoop **lineLoop);
HXTStatus hxtBoundariesGetSeamPoint(HXTBoundaries *boundaries, int *seamPoint);
HXTStatus hxtEdgesSetBoundaries(HXTEdges *edges, HXTBoundaries **boundaries);


HXTStatus hxtEdgesCreate(HXTMesh *mesh, HXTEdges** edges);
HXTStatus hxtEdgesDelete(HXTEdges **edges);

int hxtEdgesIsBoundary (HXTEdges *edges, uint32_t *e);

#endif

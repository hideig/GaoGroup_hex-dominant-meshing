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

#ifndef HEXTREME_CURVATURE_H
#define HEXTREME_CURVATURE_H

#include "hxt_mesh.h"
#include "hxt_edge.h"

/// Estimating Curvatures and Their Derivatives on Triangle Meshes, Szymon Rusinkiewicz
HXTStatus hxtCurvatureRusinkiewicz (HXTMesh *mesh, double **nodalCurvatures, double **crossField, HXTEdges* edges, int debug);

#endif

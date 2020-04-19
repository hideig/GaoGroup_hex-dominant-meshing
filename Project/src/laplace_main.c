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

#include "hxt_api.h"
#include "hxt_laplace.h"
#include "hxt_linear_system.h"
int main(int argc, char **argv) {
  if (argc != 2)
    return HXT_ERROR_MSG(HXT_STATUS_FAILED,"usage: laplace input.msh");
  HXT_CHECK(hxtInitializeLinearSystems(&argc, &argv));

  HXTContext *context;
  hxtContextCreate(&context);
  HXTMesh *mesh;
  HXT_CHECK(hxtMeshCreate(context, &mesh));
  HXT_CHECK(hxtMeshReadGmsh(mesh, argv[1]));
  HXT_CHECK(hxtLaplace(mesh));  
  return HXT_STATUS_OK;
}

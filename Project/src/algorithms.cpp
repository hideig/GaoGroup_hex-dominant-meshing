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

#include "algorithms.h"
#include "algorithms.h"

#include <stdio.h>  
#include <stdlib.h>  
#include <time.h>

#include <string>
#include <vector>

/**
* \authors Jeanne Pellerin, Kilian Verhetsel, Jonathan Lambrechts
*/

// Used to enumerate files in a directory.
#if defined WIN32
#include <io.h> 
#else
#include <glob.h>
#endif

void enumerateFiles(
  const std::string& directory,
  const std::string& extension,
  std::vector<std::string>& filenames)
{
  std::string filespecs = directory + "*." + extension;
#if defined WIN32
  _finddata_t data;
  intptr_t currentFile = _findfirst(filespecs.c_str(), &data);

  if (currentFile != -1) {
    int isTheEnd = 0;
    while (isTheEnd != -1) {
      std::string fullName = directory + data.name;
      filenames.push_back(fullName);
      isTheEnd = _findnext(currentFile, &data);
    }
    _findclose(currentFile);
  }
#else
  glob_t pglob;
  if (glob(filespecs.c_str(), 0, NULL, &pglob) == 0) {
    for (size_t i = 0; i < pglob.gl_pathc; ++i) {
      filenames.push_back(pglob.gl_pathv[i]);
    }
  }
  globfree(&pglob);
#endif
}

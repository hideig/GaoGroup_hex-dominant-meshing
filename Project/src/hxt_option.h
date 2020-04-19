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

#ifndef HXT_OPTION_H
#define HXT_OPTION_H

#include "hxt_api.h"
#include "hxt_message.h"


/* 
 * usage example:
 *
 
#include "hxt_option.h"

HXTStatus opencb(void *val, void *data) {
  //open_file(*(char**)val);
  return HXT_STATUS_OK;
}

int main(int argc, char **argv) {

  char *opts = NULL;
  char *opti = NULL;
  double optd = 3.;
  int optb = 0;

  HXTOptionList *list;
  HXT_CHECK(hxtOptionListCreate(&list, "Solve the laplace equation."));
  HXT_CHECK(hxtOptionListAddCallback(list, "filename", "file to open", HXT_OPTION_STRING, HXT_OPTION_REPEATED_ARGUMENT, &opti, opencb, NULL));
  HXT_CHECK(hxtOptionListAdd(list, "optd", "description of double option", HXT_OPTION_DOUBLE, 0, &optd));
  HXT_CHECK(hxtOptionListAdd(list, "opts", "description of string option", HXT_OPTION_STRING, 0, &opts));
  HXT_CHECK(hxtOptionListAdd(list, "optb", "description of flag option", HXT_OPTION_FLAG, 0, &optb));
  HXT_CHECK(hxtOptionListParseArgv(list, argc, argv, "laplace"));
  hxtOptionListPrintUsage(list, argv[0]);
  HXT_CHECK(hxtOptionListDelete(&list));
  ...
  }
*/

#ifdef __cplusplus
extern "C" {
#endif

typedef struct HXTOptionListStruct HXTOptionList;

typedef enum {HXT_OPTION_DEFAULT=0, HXT_OPTION_ARGUMENT=1, HXT_OPTION_REPEATED_ARGUMENT=2} HXTOptionFlag;
typedef enum {HXT_OPTION_DOUBLE, HXT_OPTION_INT, HXT_OPTION_STRING, HXT_OPTION_FLAG} hxtOptionType;

typedef HXTStatus HXTOptionCallback(void *valuep, void *data);

HXTStatus hxtOptionListCreate(HXTOptionList **list, const char *description);
HXTStatus hxtOptionListDelete(HXTOptionList **list);

HXTStatus hxtOptionListAdd(HXTOptionList *list, const char *name, const char *description, hxtOptionType otype, HXTOptionFlag oflag, void *data);
HXTStatus hxtOptionListAddCallback(HXTOptionList *list, const char *name, const char *description, hxtOptionType otype, HXTOptionFlag oflag, void *data, HXTOptionCallback *cb, void *cbdata);

HXTStatus hxtOptionListParseArgv(HXTOptionList *options, int argc, char **argv, const char *description);
HXTStatus hxtOptionListPrintUsage(const HXTOptionList *options, const char *exename);

#ifdef __cplusplus
}
#endif


#endif

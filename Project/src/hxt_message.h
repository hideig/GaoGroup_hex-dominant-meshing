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

#ifndef _MESSAGE_H_
#define _MESSAGE_H_

#include "hxt_api.h"

#include "hxt_omp.h"

#define STR(x) #x
#define STRINGIFY(x) STR(x)

#define HXT_INFO(...)                 hxtMessageInfo(__func__, __FILE__, STRINGIFY(__LINE__), ## __VA_ARGS__ )
#define HXT_INFO_COND(cond, ...)      ((cond)?HXT_INFO(__VA_ARGS__):HXT_STATUS_OK)
#define HXT_WARNING(...)              hxtMessageWarning(__func__, __FILE__, STRINGIFY(__LINE__), ## __VA_ARGS__ )

/* print an error message corresponding to a status */
#define HXT_ERROR_MSG(status, ...)    hxtMessageError(status, __func__, __FILE__, STRINGIFY(__LINE__), ## __VA_ARGS__ )
#define HXT_ERROR(status)             HXT_ERROR_MSG(status,NULL)

/* give trace message if status is not OK, but does not return */
#define HXT_TRACE_MSG(status, ...)    hxtMessageTraceError(status, __func__, __FILE__, STRINGIFY(__LINE__), ## __VA_ARGS__ )
#define HXT_TRACE(status)             HXT_TRACE_MSG(status, NULL)

/* give trace message and return the error status if not ok */
#define HXT_CHECK_MSG(status, ...)              \
  do {                                          \
    HXTStatus _tmp_ = status;                   \
    if(_tmp_<0) {                               \
      if (_tmp_>HXT_STATUS_INTERNAL) {          \
        HXT_TRACE_MSG(_tmp_, ## __VA_ARGS__);   \
      }                                         \
      return _tmp_;                             \
    }                                           \
  } while(0)
#define HXT_CHECK(status)             HXT_CHECK_MSG(status,NULL)

/* use to check some expression inside of function, throw error if exp is not true */
#ifdef DEBUG
  #define HXT_ASSERT_MSG(exp, ...)                                          \
    do {                                                                    \
      if (!(exp)) {                                                         \
        return HXT_ERROR_MSG(HXT_STATUS_ASSERTION_FAILED, ## __VA_ARGS__ ); \
    }                                                                       \
  } while(0)
  #define HXT_ASSERT(exp)             HXT_ASSERT_MSG(exp, "assertion (" #exp ") Failed")
#else
  #define HXT_ASSERT_MSG(exp, ...)
  #define HXT_ASSERT(exp)
#endif


#ifdef __cplusplus
extern "C" {
#endif

HXTStatus  hxtMessageInfo       ( const char* func, const char* file, const char* line, const char *fmt, ...);
HXTStatus  hxtMessageWarning    ( const char* func, const char* file, const char* line, const char *fmt, ...);
HXTStatus  hxtMessageError      ( HXTStatus status, const char* func, const char* file, const char* line, const char *fmt, ...);
HXTStatus  hxtMessageTraceError ( HXTStatus status, const char* func, const char* file, const char* line, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif

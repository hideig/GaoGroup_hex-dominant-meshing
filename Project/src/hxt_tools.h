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

#ifndef HXT_TOOLS_H
#define HXT_TOOLS_H

#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "hxt_message.h" // already include hxt_api (for HXT_status_t)


/* define SIMD ALIGNMENT */
#ifndef SIMD_ALIGN
#ifdef __AVX512F__
#define SIMD_ALIGN 64
#elif defined(__AVX2__)
#define SIMD_ALIGN 32
#else
#define SIMD_ALIGN 16 /* we align to 16 anyway even if no sse */
#endif
#endif

// declare alignement of pointer allocated on the stack or in a struct
#ifdef _MSC_VER
#define hxtDeclareAligned  __declspec(align(SIMD_ALIGN))
#ifndef __restrict__
#define __restrict__ __restrict
#endif
#else
#define hxtDeclareAligned  __attribute__((aligned(SIMD_ALIGN)))
#endif


#if !defined(__INTEL_COMPILER) || (__INTEL_COMPILER < 1300)
#define __assume(x)
#define __assume_aligned(x,y)
#endif


/*********************************************************
 * Hextreme malloc implementation
 *********************************************************/

/************************************************************************
  WARNING: the pointer passed (p) should be a double pointer !!

  use it like this:
___
|  int arrayLength = ...;
|  int *array;
|  HXT_CHECK( hxtMalloc(&array, sizeof(int)*arrayLength) );
|  
|  array[0] = ...;
|  [...]
|  array[arrayLength-1] = ...;
|
|  [...]
|
|  HXT_CHECK( hxtFree(&array) );
|___
 ************************************************************************/

static inline HXTStatus hxtMalloc(void* ptrToPtr, size_t size)
{
  void** p = (void**)ptrToPtr;
  *p = malloc(size);
  if (*p==NULL)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  return HXT_STATUS_OK;
}


// allocate num element of size size and zero the memory
static inline HXTStatus hxtCalloc(void* ptrToPtr, size_t num, size_t size)
{
  void** p = (void**)ptrToPtr;
  *p = calloc(num, size);
  if (*p==NULL)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  return HXT_STATUS_OK;
}


static inline HXTStatus hxtFree(void* ptrToPtr)
{
  void** p = (void**)ptrToPtr;
  free(*p);
  *p=NULL;
  return HXT_STATUS_OK;
}


static inline HXTStatus hxtRealloc(void* ptrToPtr, size_t size)
{
  void** p = (void**)ptrToPtr;
  void* newptr = realloc(*p, size);
  if (newptr==NULL && *p!=NULL && size!=0)
    return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
  *p = newptr;
  return HXT_STATUS_OK;
}


/*********************************************************
 * Hextreme aligned malloc implementation
 *********************************************************/
static inline HXTStatus hxtGetAlignedBlockSize(void* ptrToPtr, size_t* size)
{
  char** p2 = *(char***)(ptrToPtr);

  if(p2==NULL){
    *size = 0;
    return HXT_STATUS_OK;
  }

  *size = p2[-2]-p2[-1];

  if((*size^0xBF58476D1CE4E5B9ULL) != ((size_t)p2[-3]-(size_t)p2[-1])){
    return HXT_ERROR(HXT_STATUS_POINTER_ERROR);
  }

  return HXT_STATUS_OK;
}

static inline HXTStatus hxtAlignedMalloc(void* ptrToPtr, size_t size)
{
  char** p2;
  char* pstart;
  // malloc is at least aligned on sizeof(size_t), so the maximum shift is SIMD_ALIGN-sizeof(size_t)
  const size_t startOffset = SIMD_ALIGN - sizeof(size_t) + 3*sizeof(char*); // up to 80 bytes...

  const size_t safe_size = size + startOffset;
  if(safe_size<size) // account for possible overflow
    goto error;

  pstart = (char*) malloc(safe_size); // start of original block
  if(pstart==NULL)
    goto error;

  p2 = (char**)(((size_t)(pstart) + startOffset) & ~(SIMD_ALIGN - 1)); // only keep bits ge SIMD_ALIGN
  p2[-1] = pstart;
  p2[-2] = pstart + size;
  p2[-3] = pstart + (size ^ 0xBF58476D1CE4E5B9ULL); // makes a verification possible
  *(void**)ptrToPtr = p2;

  return HXT_STATUS_OK;

error:
  *(void**)ptrToPtr = NULL;
  return HXT_ERROR(HXT_STATUS_OUT_OF_MEMORY);
}


static inline HXTStatus hxtAlignedFree(void* ptrToPtr)
{
  void** p = (void**)ptrToPtr;

  if(*p!=NULL){
    char** p2 = (char**)(*p);
    if(((size_t)(p2[-2]-p2[-1])^0xBF58476D1CE4E5B9ULL) != (size_t)(p2[-3]-p2[-1]))
      return HXT_ERROR(HXT_STATUS_POINTER_ERROR);

    p2[-3] = NULL;
    free(p2[-1]);
    *p = NULL;
  }

  return HXT_STATUS_OK;
}


static inline HXTStatus hxtAlignedRealloc(void* ptrToPtr, size_t size)
{
  if (size == 0) {
    HXT_CHECK(hxtAlignedFree(ptrToPtr));
    return HXT_STATUS_OK;
  }
  
  size_t old_size;
  HXT_CHECK( hxtGetAlignedBlockSize(ptrToPtr, &old_size) );

  if(size>old_size || size+4096<old_size){ // we do not shrink block to gain less than 4096 bytes
    void* newptr = NULL;
    HXT_CHECK( hxtAlignedMalloc(&newptr, size));

    memcpy(newptr, *(void**)ptrToPtr, size>old_size?old_size:size);

    HXT_CHECK( hxtAlignedFree(ptrToPtr) );

    *(void**)ptrToPtr = newptr;
  }

  return HXT_STATUS_OK;
}


/*********************************************************
 * Matrix operations
 *********************************************************/
HXTStatus hxtDet3x3(double mat[3][3], double *det);
HXTStatus hxtInv3x3(double mat[3][3], double inv[3][3], double *det);

HXTStatus hxtInv4x4ColumnMajor(double mat[16], double inv[16], double *det);

#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327950
#endif // !M_PI

#ifdef __cplusplus
}
#endif


#endif

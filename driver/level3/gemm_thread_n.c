/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "common.h"
#include "pagemap.h"

//#define LOCK

int CNAME(int mode, blas_arg_t *arg, BLASLONG *range_m, BLASLONG *range_n, int (*function)(blas_arg_t*, BLASLONG*, BLASLONG*,FLOAT *, FLOAT *, BLASLONG), void *sa, void *sb, BLASLONG nthreads) {

  blas_queue_t queue[MAX_CPU_NUMBER];
  BLASLONG range[MAX_CPU_NUMBER + 1];

  BLASLONG width, i, num_cpu;

  if (!range_n) {
    range[0] = 0;
    i        = arg -> n;
  } else {
    range[0] = range_n[0];
    i        = range_n[1] - range_n[0];
  }

  num_cpu  = 0;

  while (i > 0){

    width  = blas_quickdivide(i + nthreads - num_cpu - 1, nthreads - num_cpu);

    i -= width;
    if (i < 0) width = width + i;

    range[num_cpu + 1] = range[num_cpu] + width;

    queue[num_cpu].mode    = mode;
    queue[num_cpu].routine = function;
    queue[num_cpu].args    = arg;
    queue[num_cpu].range_m = range_m;
    queue[num_cpu].range_n = &range[num_cpu];
#if 1   //defined(LOONGSON3A)
    queue[num_cpu].sa      = sa	+ GEMM_OFFSET_A * num_cpu;
    queue[num_cpu].sb      = sa + GEMM_OFFSET_A * (nthreads + 1) + GEMM_OFFSET_B * num_cpu;
    //printf("%lx %lx\n", queue[num_cpu].sa, queue[num_cpu].sb);
#else
    queue[num_cpu].sa      = NULL;
    queue[num_cpu].sb      = NULL;
#endif
    queue[num_cpu].next    = &queue[num_cpu + 1];
    num_cpu ++;
  }

  if (num_cpu) {
#if 1 //defined(LOONGSON3A)
    queue[0].sa = sa;
    queue[0].sb = sa + GEMM_OFFSET_A * (nthreads + 1);
#else
	queue[0].sa = sa;
	queue[0].sb = sb;
#endif
    queue[num_cpu - 1].next = NULL;

#ifdef LOCK
    pid_t pid = getpid();
    uintptr_t paddr = 0;
    int mem_fd;
    void* page = sa + GEMM_OFFSET_A * (nthreads + 1);
    //void* page = sa;
    if (lkmc_pagemap_virt_to_phys_user(&paddr, pid, ((uintptr_t)page))) {
        return -1;
    }
    //printf("virt address: %lx , phys address: %lx\n", (uintptr_t)page, paddr);
    unsigned long long * base;
    if((mem_fd = open("/dev/mem", O_RDWR|O_SYNC)) == -1) {
      printf("open /dev/mem failed!\n");
      return 0;
    }
    base = (unsigned long long *)mmap(0, 0x280, PROT_READ | PROT_WRITE,
          MAP_SHARED, mem_fd, ((WhereAmI() >> 4) << 44 ) | 0x1fe00000);
    if (base == MAP_FAILED) {
         printf("mmap for /dev/mem failed!\n");
         return 0;
    }
    /* cpu0 ~ cpu3 */
    *(base + 0x40 + 0) = ((unsigned long long)paddr | 0x8000000000000000UL);
    *(base + 0x48 + 0) = ~(unsigned long long)(1024 * 1024 * 7 - 1);
    /* cpu4 ~ cpu7 */
    *(base + 0x40 + 1) = ((unsigned long long)(paddr + 4 * GEMM_OFFSET_B) | 0x8000000000000000UL);
    *(base + 0x48 + 1) = ~(unsigned long long)(1024 * 1024 * 7 - 1);
    /* cpu8 ~ cpu11 */
    *(base + 0x40 + 2) = ((unsigned long long)(paddr + 8 * GEMM_OFFSET_B) | 0x8000000000000000UL);
    *(base + 0x48 + 2) = ~(unsigned long long)(1024 * 1024 * 7 - 1);
    /* cpu12 ~ cpu15 */
    *(base + 0x40 + 3) = ((unsigned long long)(paddr + 12 * GEMM_OFFSET_B) | 0x8000000000000000UL);
    *(base + 0x48 + 3) = ~(unsigned long long)(1024 * 1024 * 7 - 1);
    //printf("lock virt address: %lx, %lx, %lx\n", base, *(base + 0x40), *(base + 0x48));
#endif


    exec_blas(num_cpu,
	      queue);
#ifdef LOCK
    *(base + 0x40 + 0) = ((unsigned long long)paddr)&(0x0UL) ;   //base_reg
    *(base + 0x40 + 1) = ((unsigned long long)paddr)&(0x0UL) ;   //base_reg
    *(base + 0x40 + 2) = ((unsigned long long)paddr)&(0x0UL) ;   //base_reg
    *(base + 0x40 + 3) = ((unsigned long long)paddr)&(0x0UL) ;   //base_reg
    if(munmap(base, 0x280))
    {
        printf("munmap for /dev/mem failed!\n");
        return 0;
    }
    close(mem_fd);
#endif
  }

  return 0;
}

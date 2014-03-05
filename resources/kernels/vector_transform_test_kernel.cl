/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file vector_transform_test_kernel.cl
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
//#pragma OPENCL EXTENSION cl_khr_fp16 : enable

// disable dbg_printf for GPU
#define dbg_printf(format, ...)

// enable printf for CPU
//#pragma OPENCL EXTENSION cl_amd_printf : enable
//#define dbg_printf(format, ...) printf(format, ##__VA_ARGS__)

__kernel void testKernel(
    __global float* inputValuesX,
    __global float* inputValuesY,
    __global float* inputValuesZ,
    __global float* outputValuesX,
    __global float* outputValuesY,
    __global float* outputValuesZ)
{
    dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

    unsigned int i = get_global_id(0);
    //unsigned int global_size = get_global_size(0);

    // evaluate the function
    float4 theVector = (float4)(inputValuesX[i], inputValuesY[i], inputValuesZ[i], -1.2345f);

    evaluateVectorTransform(&theVector);

    if (theVector.w != -1.2345f) {
        outputValuesX[i] = 8888888.f;
        outputValuesY[i] = 8888888.f;
        outputValuesZ[i] = 8888888.f;
    } else {
        outputValuesX[i] = theVector.x;
        outputValuesY[i] = theVector.y;
        outputValuesZ[i] = theVector.z;
    }

    dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    dbg_printf("Kernel finished.\n");
}

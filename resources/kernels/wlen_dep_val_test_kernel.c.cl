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
 * @file wlen_dep_val_test_kernel.cl
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

__kernel void testKernel(__global float* xValues,
                         __global float* yValues,
                         uint mode)
{
    //dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

    unsigned int i = get_global_id(0);
    unsigned int global_size = get_global_size(0);

    // evaluate the function
    if (mode==0) {
        yValues[i] = evaluateFunction(xValues[i]);
    } else if (mode==1) {
        yValues[i] = evaluateDerivative(xValues[i]);
    } else {
        yValues[i] = 9999999.f;
    }

    //dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    //dbg_printf("Kernel finished.\n");
}

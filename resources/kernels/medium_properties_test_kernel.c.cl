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
 * @file medium_properties_test_kernel.cl
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */


__kernel void testKernel(__global ulong* MWC_RNG_x,
                         __global uint* MWC_RNG_a,
                         __global float* xValues,
                         __global float* yValues,
                         uint layer,
                         uint mode)
{
    //dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

    unsigned int i = get_global_id(0);
    unsigned int global_size = get_global_size(0);

    //download MWC RNG state
    ulong real_rnd_x = MWC_RNG_x[i];
    uint real_rnd_a = MWC_RNG_a[i];
    ulong *rnd_x = &real_rnd_x;
    uint *rnd_a = &real_rnd_a;

    // evaluate the function
    if (mode==0) {
        yValues[i] = getPhaseRefIndex(layer, xValues[i]);
#ifndef NO_DISPERSION
    } else if (mode==1) {
        yValues[i] = getDispersion(layer, xValues[i]);
#endif
    } else if (mode==2) {
        yValues[i] = getGroupVelocity(layer, xValues[i]);
    } else if (mode==3) {
        yValues[i] = getAbsorptionLength(layer, xValues[i]);
    } else if (mode==4) {
        yValues[i] = getScatteringLength(layer, xValues[i]);
    } else if (mode==5) {
        // this ignores the input data and just generates random numbers
        yValues[i] = makeScatteringCosAngle(RNG_ARGS_TO_CALL);
    } else {
        yValues[i] = 9999999.f;
    }

    //dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    //dbg_printf("Kernel finished.\n");

    //upload MWC RNG state
    MWC_RNG_x[i] = real_rnd_x;
    MWC_RNG_a[i] = real_rnd_a;
}

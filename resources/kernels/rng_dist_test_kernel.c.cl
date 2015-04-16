__kernel void testKernel(__global float* randomNumbers,
                         __global ulong* MWC_RNG_x,
                         __global uint* MWC_RNG_a)
{
    //dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

    unsigned int i = get_global_id(0);
    //unsigned int global_size = get_global_size(0);

    //download MWC RNG state
    ulong real_rnd_x = MWC_RNG_x[i];
    uint real_rnd_a = MWC_RNG_a[i];
    ulong *rnd_x = &real_rnd_x;
    uint *rnd_a = &real_rnd_a;

    // call the function and generate a random number according to the requested distribution
#ifdef DISTRIBUTION_ARGS
    randomNumbers[i] = generateRandomNumberAccordingToDistribution(RNG_ARGS_TO_CALL, DISTRIBUTION_ARGS);
#else
    randomNumbers[i] = generateRandomNumberAccordingToDistribution(RNG_ARGS_TO_CALL);
#endif
    
    //dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    //dbg_printf("Kernel finished.\n");

    //upload MWC RNG state
    MWC_RNG_x[i] = real_rnd_x;
    MWC_RNG_a[i] = real_rnd_a;
}

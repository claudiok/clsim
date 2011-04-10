#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
//#pragma OPENCL EXTENSION cl_khr_fp16 : enable

// disable dbg_printf for GPU
#define dbg_printf(format, ...)

// enable printf for CPU
//#pragma OPENCL EXTENSION cl_amd_printf : enable
//#define dbg_printf(format, ...) printf(format, ##__VA_ARGS__)

__kernel void testKernel(__read_only __global float* xValues,
                         __write_only __global float* yValues,
                         uint mode)
{
    dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

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

    dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    dbg_printf("Kernel finished.\n");
}

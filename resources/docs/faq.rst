..
.. Copyright (c) 2011, 2012
.. Claudio Kopper <claudio.kopper@icecube.wisc.edu>
.. and the IceCube Collaboration <http://www.icecube.wisc.edu>
..
.. Permission to use, copy, modify, and/or distribute this software for any
.. purpose with or without fee is hereby granted, provided that the above
.. copyright notice and this permission notice appear in all copies.
..
.. THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
.. WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
.. MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
.. SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
.. WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
.. OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
.. CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
..
..
.. $Id$
..
.. @file index.rst
.. @version $Revision$
.. @date $Date$
.. @author Claudio Kopper
..

.. highlight:: python

Frequently asked questions and common pitfalls
==============================================

1. Why do I get angry emails from my cluster adminstrator about using too many threads
or causing a really high load when submitting clsim jobs on a cluster (when running on
CPUs)?

    All existing OpenCL implementations for CPUs (either Intel's or AMD's version at the
    moment) will parallelize photon propagation, which means that they will run multiple
    threads for a single clsim process. On multi-core systems, CPU usages of 1000% or more
    have been observed. This is great for getting results quickly on an interactive system,
    but it is a bad idea when using it on cluster nodes where multiple processes are
    being run at the same time. (Almost all of the existing batch schedulers assume that
    a process will only ever use a single CPU core/thread.)
    
    In order to configure clsim to only use a single core for photon propagation, the
    standard tray segment provides the option ``DoNotParallelize`` which you should
    set to ``True`` when submitting to a cluster::
    
        tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
                        DoNotParallelize=True,
                        UseGPUs=False, UseCPUs=True)

2. How can I make sure that my clsim instance is using the correct GPU in a multi-GPU system?

    By default, clsim will take over all available GPUs whenever the ``UseGPUs`` option is
    set to ``True``. The easiest way to select a specific GPU on cluster nodes that have
    multiple GPUs available rely on the submitting users to restrict their jobs to a certain
    GPU is using the ``CUDA_VISIBLE_DEVICES`` environment variable. (This is specific to
    nVidida and will not work with AMD cards.) You just set 
    
    .. code-block:: bash
    
        export CUDA_VISIBLE_DEVICES="0,2,3"
    
    This example will make clsim use devices number 0,2,3 and will skip device number 1.
    (The numbering starts from 0.) The process will actually only see three of the four
    hypothetical devices.
    
    On a well-configured cluster, this variable should already be set for each running job.

3. So I'm running on a condor cluster that does not set CUDA_VISIBLE_DEVICES automatically
(npx3 in Madison), what's a good way to set it? I do get something called a "_CONDOR_SLOT"
which is a number starting at 1 that I am supposed to use to select the GPU.

    You can try adding the following block of python code to your script before any of the
    IceTray stuff. It will substract 1 from the slot id and set CUDA_VISIBLE_DEVICES
    to use the appropriate GPU.::
    
        import os
        if "_CONDOR_SLOT" in os.environ: # running in condor?
            if "CUDA_VISIBLE_DEVICES" in os.environ:
                print "running in CONDOR, but CUDA_VISIBLE_DEVICES is already set. no further configuration necessary."
            else:
                condorSlotNumber = int(os.environ["_CONDOR_SLOT"])
                print "script seems to be running in condor (slot %u). auto-configuring CUDA_VISIBLE_DEVICES!" % condorSlotNumber
                os.environ["CUDA_VISIBLE_DEVICES"] = str(condorSlotNumber-1)

4. I am running clsim on an nVidia GPU, but it seems to hang and is not generating photons.
When looking at the GPU utilization using ``nvidia-smi``, I do see 100%, but nothing seems
to be happening.
    
    This is a known bug in very old OpenCL libraries supplied by nVidia. clsim is known not
    to work with driver versions 260.19.21 and older. Versions 270.x.x and newer have been tested
    and are known to work, so just update to the most recent version available and clsim should work.

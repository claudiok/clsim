#
# Copyright (c) 2011, 2012
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 
# 
# $Id: AsyncTap.py 108199 2013-07-12 21:33:08Z nwhitehorn $
# 
# @file AsyncTap.py
# @version $Revision: 108199 $
# @date $Date: 2013-07-12 17:33:08 -0400 (Fri, 12 Jul 2013) $
# @author Claudio Kopper
#

from __future__ import print_function

# ZMQ is not functional yet
#try:
#    import zmq
#    has_zmq = True
#except ImportError:
#    has_zmq = False
has_zmq = False

import multiprocessing
import signal
import sys
import time

from icecube import icetray, dataclasses
from I3Tray import *


# run a writer tray, getting frames from a queue (this runs as a subprocess)
def RunAsyncTray(queue,segment,segmentArgs,childIsReady,debug):
    if has_zmq:
        zmq_context = zmq.Context()
    else:
        zmq_context = None
    
    # just pushes frame on a queue
    class AsyncReceiver(icetray.I3Module):
        def __init__(self, context):
            icetray.I3Module.__init__(self, context)
            self.AddParameter("Debug", "Output some status information for debugging", False)
            self.AddParameter("Queue", "", None)
            self.AddParameter("ZMQContext", "", None)
            self.AddParameter("ChildIsReady", "", None)
            self.AddOutBox("OutBox")
            self.nframes = 0
        def Configure(self):
            self.debug = self.GetParameter("Debug")
            self.queue = self.GetParameter("Queue")
            self.zmq_context = self.GetParameter("ZMQContext")
            self.childIsReady = self.GetParameter("ChildIsReady")
            
            if has_zmq:
                socketString = self.queue
                if self.debug:
                    print("setting up zmq to receive on", socketString)
                    sys.stdout.flush()
                self.queue = self.zmq_context.socket(zmq.PULL)
                if self.debug:
                    print("waiting for zmq socket connection")
                    sys.stdout.flush()
                self.queue.connect(socketString)
                if self.debug:
                    print("zmq ready")
                    sys.stdout.flush()
            
            # signal eneryone that we are ready to receive data
            if self.debug:
                print("setting child process state to 1 (ready)")
                sys.stdout.flush()
            self.childIsReady.value = 1
            if self.debug:
                print("child process is ready to receive data.")
                sys.stdout.flush()
            
        def Process(self):
            #print "waiting for frame..", self.nframes
            #sys.stdout.flush()
            if has_zmq:
                frame = self.queue.recv_pyobj()
            else:
                frame = self.queue.get()
            #if frame is not None:
            #    #print "frame received", self.nframes , frame.Stop
            #    #sys.stdout.flush()
        
            if frame is None:
                if self.debug: 
                    print("requesting suspension")
                    sys.stdout.flush()
                self.RequestSuspension()
            else:
                self.nframes += 1
                self.PushFrame(frame)
        def Finish(self):
            if self.debug:
                print("received", self.nframes, "frames")
                sys.stdout.flush()
    
    class RerouteSigIntToDefault(icetray.I3Module):
        def __init__(self, context):
            icetray.I3Module.__init__(self, context)
            self.AddOutBox("OutBox")
        def Configure(self):
            # Reroute ctrl+C (sigint) to the default handler.
            # (This needs to be done in a module, as IceTray re-routes
            # the signal in tray.Execute())
            signal.signal(signal.SIGINT, signal.SIG_IGN)
        def DAQ(self, frame):
            self.PushFrame(frame)
    
    writeTray = I3Tray()
    writeTray.AddModule(AsyncReceiver, "theAsyncReceiver",
        Queue=queue, ZMQContext=zmq_context, Debug=debug, ChildIsReady=childIsReady)
    writeTray.AddModule(RerouteSigIntToDefault, "rerouteSigIntToDefault")
    
    if hasattr(segment, '__i3traysegment__'):
        writeTray.AddSegment(segment,"theSegment", **segmentArgs)
    else:
        writeTray.AddModule(segment,"theModule", **segmentArgs)
    
    

    if debug:
        print("worker starting..")
        sys.stdout.flush()
    writeTray.Execute()
    writeTray.Finish()
    if debug:
        print("worker finished.")
        sys.stdout.flush()



class AsyncTap(icetray.I3ConditionalModule):
    """
    Starts a module or a segment on its own tray in its own process
    and pipes frames from the current tray to the child tray.
    Can be used for things like asynchronous file writing.
    The frames will never get back to the master module,
    so this is effecively a "fork".
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("Debug", "Output some status information for debugging", False)
        self.AddParameter("BufferSize", "How many frames should be buffered", 2000)
        self.AddParameter("Segment", "The tray segment to run asynchronously", None)
        self.AddParameter("Args", "A dictionary with keyword arguments for the asynchronous segment.", dict())
        self.AddOutBox("OutBox")
        self.nframes = 0
        self.suspensionRequested=False

    def CheckOnChildProcess(self):
        if self.suspensionRequested: return False
        
        # check to see if child process is still alive
        if not self.process.is_alive():
            print("ERROR: ****** child tray died unexpectedly. Terminating master tray. ******")
            self.RequestSuspension()
            self.suspensionRequested=True
            return False
        else:
            return True

    def Configure(self):
        self.debug = self.GetParameter("Debug")
        self.buffersize = self.GetParameter("BufferSize")
        self.segment = self.GetParameter("Segment")
        self.args = self.GetParameter("Args")

        if self.debug:
            print("starting child process..")
            sys.stdout.flush()
        self.childIsReady = multiprocessing.Value('i', 0)
        
        if has_zmq:
            #socketString = "ipc:///tmp/asynctap.0"
            socketString = "tcp://127.0.0.1:5557"
            self.process = multiprocessing.Process(target=RunAsyncTray, args=(socketString,self.segment,self.args,self.childIsReady,self.debug,))
        else:
            self.queueToProcess = multiprocessing.Queue(self.buffersize)
            self.process = multiprocessing.Process(target=RunAsyncTray, args=(self.queueToProcess,self.segment,self.args,self.childIsReady,self.debug,))


        if has_zmq:
            if self.debug:
                print("binding to zmq PUSH socket", socketString)
                sys.stdout.flush()
            self.zmq_context = zmq.Context()
            self.queueToProcess = self.zmq_context.socket(zmq.PUSH)
            self.queueToProcess.bind(socketString)
            if self.debug:
                print("zmq socket is bound.")
                sys.stdout.flush()

        self.process.start()
        if self.debug:
            print("child process running.")
            sys.stdout.flush()

        self.CheckOnChildProcess()
        
        if self.childIsReady.value == 0:
            while True:
                if not self.process.is_alive():
                    print("*** child process died unexpectedly")
                    raise RuntimeError("Child process died unexpectedly")
                    
                if self.debug:
                    print("child process not ready yet, waiting..")
                    sys.stdout.flush()
                time.sleep(1) # sleep 1 second
                if self.childIsReady.value != 0:
                    break
        if self.debug:
            print("child process ready to receive data.")
            sys.stdout.flush()

    def Process(self):
        frame = self.PopFrame()
        if self.debug: 
            print("sending frame..", self.nframes , frame.Stop)
            sys.stdout.flush()
        if self.CheckOnChildProcess():
            if has_zmq:
                self.queueToProcess.send_pyobj(frame)
            else:
                self.queueToProcess.put(frame)
        if self.debug: 
            print("frame sent", self.nframes , frame.Stop)
            sys.stdout.flush()
        self.nframes += 1
        self.PushFrame(frame)

    def Finish(self):
        if self.debug:
            print("sent", self.nframes, "frames")
            sys.stdout.flush()
        if self.CheckOnChildProcess():
            if has_zmq:
                self.queueToProcess.send_pyobj(None) # signal receiver to quit
            else:
                self.queueToProcess.put(None) # signal receiver to quit
        
        if self.debug:
            print("async module finished. waiting for child tray..")
            sys.stdout.flush()
        self.process.join()
        if self.debug:
            print("child tray finished.")
            sys.stdout.flush()


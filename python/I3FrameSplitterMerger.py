from icecube import icetray, dataclasses, simclasses
from copy import copy

class I3FrameMCPEMerger(icetray.I3ConditionalModule):
    """
    IceTray Python I3Module to merge frames and frameobjects
    that were broken apart because they cause clsim to eat
    all the RAM
    """
    def __init__(self, context):
        """
        Standard IceTray module __init__ defines 
        I3Module parameters
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('InputMCPESeriesMapName',
                          'Name of input I3MCPESeriesMapName',
                          'I3MCPESeriesMapName')
        self.AddParameter('OriginalFrameStream',
                          "Stream the original frame is now",
                          icetray.I3Frame.Stream('q'))
        self.AddParameter("SplitFrameStream",
                          "Stream of the split frames",
                          icetray.I3Frame.Stream('Q'))
        self.AddOutBox('OutBox')

    def Configure(self):
        self.input_mcpeseries_name = self.GetParameter("InputMCPESeriesMapName")
        self.original_stream = self.GetParameter("OriginalFrameStream")
        self.split_stream = self.GetParameter("SplitFrameStream")
        self.original_frame = None
        self.frame_count = 0

    def Process(self):
        """
        Process() to merge the frames. It will Flush the accumulated frame out
        if any frame that is not the split stream stop.
        """
        frame = self.PopFrame()
        if frame.Stop not in [self.original_stream, self.split_stream]:
            if self.frame_count > 0:
                self.Flush()
                self.frame_count = 0
            self.PushFrame(frame)
        if frame.Stop == self.original_stream:
            if self.frame_count > 0:
                self.Flush()
                self.frame_count = 0
            self.original_frame = I3FrameStreamChanger.ChangeFrame(frame, self.split_stream)
            self.frame_count += 1
            self.new_mcpeseries = simclasses.I3MCPESeriesMap()
        if frame.Stop == self.split_stream:
            self.new_mcpeseries.merge_nosort(frame[self.input_mcpeseries_name])
            del frame
    
    def Finish(self):
        self.Flush()
            
    def Flush(self):
        """
        Function that flushes the combined frame out. It will sort the
        new I3MCPESeriesMap before putting it in the original frame.
        """
        self.new_mcpeseries.sort()
        self.original_frame[self.input_mcpeseries_name] = self.new_mcpeseries
        self.PushFrame(self.original_frame)

# class I3FrameMCPEMerger(icetray.I3PacketModule):
#     def __init__(self, ctx):
#         icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.Stream('q'))
#         self.AddParameter('OriginalFrameStream',
#                           "Stream the original frame is now",
#                           icetray.I3Frame.Stream('q'))
#         self.AddParameter("SplitFrameStream",
#                           "Stream of the split frames",
#                           icetray.I3Frame.Stream('Q'))
#         self.AddParameter('InputMCPESeriesMapName',
#                           'Name of input I3MCPESeriesMapName',
#                           'I3MCPESeriesMapName')
#         self.AddOutBox("OutBox")
#
#     def Configure(self):
#         icetray.I3PacketModule.Configure(self)
#         self.original_stream = self.GetParameter("OriginalFrameStream")
#         self.split_stream = self.GetParameter("SplitFrameStream")
#         self.input_mcpeseries_name = self.GetParameter("InputMCPESeriesMapName")
#         self.packet_types = [self.original_stream, self.split_stream]
#
#     def FramePacket(self, frames):
#         original_frame = None
#         mcpe_series_map = None
#         for frame in frames:
#             if frame.Stop == self.original_stream:
#                 original_frame = I3FrameStreamChanger.ChangeFrame(frame, self.split_stream)
#             if frame.Has(self.input_mcpeseries_name):
#                 if not mcpe_series_map:
#                     mcpe_series_map = frame[self.input_mcpeseries_name]
#                 else:
#                     mcpe_series_map.merge_nosort(frame[self.input_mcpeseries_name])
#             # if frame.Stop == self.split_stream:
#             #     self.PushFrame(I3FrameStreamChanger.ChangeFrame(frame, self.original_stream))
#         mcpe_series_map.sort()
#         original_frame[self.input_mcpeseries_name] = mcpe_series_map
#         self.PushFrame(original_frame)

class I3FrameEnergySplitter(icetray.I3ConditionalModule):
    """
    IceTray Python I3Module that breaks apart frames into smaller frames. 
    Needed for ultra high energy events to decrease the memory usage
    """
    def __init__(self, context):
        """
        Standard IceTray module __init__ defines 
        I3Module parameters
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('InputMCTreeName',
                          'Name of input I3MCTree',
                          'I3MCTree')
        self.AddParameter('EnergyPerFrame',
                          'Energy per frame',
                          20.*icetray.I3Units.PeV)
        self.AddParameter('DummyStream',
                          "The stream this module operates on",
                          icetray.I3Frame.Stream('q'))
        
        self.AddOutBox('OutBox')
    
    def Configure(self):
        self.input_mctree_name = self.GetParameter("InputMCTreeName")
        self.energy_per_frame = self.GetParameter("EnergyPerFrame")
        self.dummy_stream = self.GetParameter("DummyStream")        
        # self.Register(self.stream, self.Splitter)
        
    def DAQ(self, frame):
        """
        DAQ() function that takes the frames and breaks them apart
        trying to keep the maximum energy per frame below a 
        set threshold
        """
        i = 0
        frame_list = []
        mctree = frame[self.input_mctree_name]
        # self.PushFrame()
        frame_list.append(I3FrameStreamChanger.ChangeFrame(frame, self.dummy_stream))
        energy_counter = 0.
        new_frame = icetray.I3Frame(icetray.I3Frame.DAQ)
        new_mctree = dataclasses.I3MCTree()
        dummy_primary = dataclasses.I3Particle(dataclasses.I3Particle.ParticleShape.Primary)
        new_mctree.add_primary(dummy_primary)
        for p in mctree:
            if p.shape == dataclasses.I3Particle.Dark or\
               p.location_type != dataclasses.I3Particle.InIce or\
               p.is_neutrino:
               continue
            new_mctree.append_child(dummy_primary, p)
            if p.type == dataclasses.I3Particle.MuMinus or \
               p.type == dataclasses.I3Particle.MuPlus:
                   energy_counter +=  p.length * 0.3*icetray.I3Units.GeV
            else:
                energy_counter += p.energy
            if energy_counter >= self.energy_per_frame:
                icetray.logging.log_info("Energy per frame %f" % (energy_counter / (1.*icetray.I3Units.PeV)))
                new_frame[self.input_mctree_name] = new_mctree
                i += 1
                # self.PushFrame(new_frame)
                self.PushFrame(new_frame)
                # frame_list.append(new_frame)
                new_frame = icetray.I3Frame(icetray.I3Frame.DAQ)
                new_mctree = dataclasses.I3MCTree()
                new_mctree.add_primary(dummy_primary)
                energy_counter = 0.
        icetray.logging.log_info("Energy per frame %f" % (energy_counter / (1.*icetray.I3Units.PeV)))
        # print (energy_counter / (1.*icetray.I3Units.PeV))
        new_frame[self.input_mctree_name] = new_mctree
        # self.PushFrame(new_frame)
        # frame_list.append(new_frame)
        self.PushFrame(new_frame)
        i += 1
        # for f in frame_list:
        #     if f.Stop == self.dummy_stream:
        #         f["SubFrameCount"] = icetray.I3Int(i)
        #     self.PushFrame(f)
        icetray.logging.log_info("--------------------------------------------------------")
        
        # self.PushFrame(icetray.I3Frame(icetray.I3Frame.Stream(",")))


class I3FrameStreamChanger(icetray.I3ConditionalModule):
    """
    Changes the stream identifier of P-frames
    """
    def __init__(self, ctx):
        super(I3FrameStreamChanger, self).__init__(ctx)
        self.AddParameter("NewStream",
                          "The new stream/stop for frames",
                          icetray.I3Frame.Stream('q'))
        self.AddParameter("OldStream",
                          "The old stream/stop for frames",
                          icetray.I3Frame.Stream('Q'))
        self.AddOutBox("OutBox")

    def Configure(self):
        self.new_stream = self.GetParameter("NewStream")
        self.old_stream = self.GetParameter("OldStream")
        self.Register(self.old_stream, self.Changer)

    def Changer(self, frame):
        if frame.Stop != self.old_stream:
            self.PushFrame(frame)
        self.PushFrame(I3FrameStreamChanger.ChangeFrame(frame, self.new_stream))
    
    @staticmethod
    def ChangeFrame(frame, new_stream):
        frame.purge() # deletes all non-native items
        
        for name in frame.keys():
            frame.change_stream(name, new_stream)
        
        new_frame = icetray.I3Frame(new_stream)
        new_frame.merge(frame)
        del frame
        return new_frame 
        
        # self.PushFrame(new_frame)
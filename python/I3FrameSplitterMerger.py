from icecube import icetray, dataclasses
from copy import copy

# class I3FrameMCPEMerger(icetray.I3ConditionalModule):
#     """
#
#     """
#     def __init__(self, context):
#         icetray.I3ConditionalModule.__init__(self, context)
#         self.AddParameter('InputMCPESeriesMapName',
#                           'Name of input I3MCPESeriesMapName',
#                           'I3MCPESeriesMapName')
#         self.AddParameter('OriginalFrameStream',
#                           "Stream the original frame is now",
#                           icetray.I3Frame.Stream('q'))
#         self.AddParameter("SplitFrameStream",
#                           "Stream of the split frames",
#                           icetray.I3Frame.Stream('Q'))
#         self.AddParameter("StopFrameStream",
#                           "Stream that indicates the end of split streams",
#                           icetray.I3Frame.Stream(','))
#         self.AddOutBox('OutBox')
#
#     def Configure(self):
#
#         self.input_mcpeseries_name = self.GetParameter("InputMCPESeriesMapName")
#         self.original_stream = self.GetParameter("OriginalFrameStream")
#         self.split_stream = self.GetParameter("SplitFrameStream")
#         self.stop_stream = self.GetParameter("StopFrameStream")
#         self.original_frame = None
#         self.split_frames = []
#
#     def Process(self):
#
#         frame = self.PopFrame()
#         if frame.Stop not in [self.original_stream, self.split_stream, self.stop_stream]:
#             self.PushFrame(frame)
#         if frame.Stop == self.original_stream:
#             self.original_frame = I3FrameStreamChanger.ChangeFrame(frame, self.split_stream)
#         if frame.Stop == self.split_stream:
#             self.split_frames.append(frame)
#         if frame.Stop == self.stop_stream:
#             for i, f in enumerate(self.split_frames):
#                 if i == 0:
#                     new_mcpeseriesmap = copy.copy(f[self.input_mcpeseries_name])
#                 else:
#                     new_mcpeseriesmap.merge_nosort(f[self.input_mcpeseries_name])
#             new_mcpeseriesmap.sort()
#             self.original_frame[self.input_mcpeseries_name] = new_mcpeseriesmap
#         self.PushFrame(self.original_frame)
        
        
        # print frame

class I3FrameMCPEMerger(icetray.I3PacketModule):
    def __init__(self, ctx):
        icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.Stream('q'))
        self.AddParameter('OriginalFrameStream',
                          "Stream the original frame is now",
                          icetray.I3Frame.Stream('q'))
        self.AddParameter("SplitFrameStream",
                          "Stream of the split frames",
                          icetray.I3Frame.Stream('Q'))
        self.AddParameter('InputMCPESeriesMapName',
                          'Name of input I3MCPESeriesMapName',
                          'I3MCPESeriesMapName')
        self.AddOutBox("OutBox")
    
    def Configure(self):
        icetray.I3PacketModule.Configure(self)
        original_stream = self.GetParameter("OriginalFrameStream")
        split_stream = self.GetParameter("SplitFrameStream")
        self.input_mcpeseries_name = self.GetParameter("InputMCPESeriesMapName")
        self.packet_types = [original_stream, split_stream]
        self.original_frame = None
        self.mcpe_series_map = None

    def FramePacket(self, frames):
        print len(frames)
        for frame in frames:
            if frame.Stop == icetray.I3Frame.Stream('q'):
                self.original_frame = I3FrameStreamChanger.ChangeFrame(frame, self.split_stream)
            if frame.Has("self.input_mcpeseries_name"):
                if not self.mcpe_series_map:
                    self.mcpe_series_map = frame[self.input_mcpeseries_name]
                else:
                    self.mcpe_series_map.merge_nosort(frame[self.input_mcpeseries_name])
        self.mcpe_series_map.sort()
        self.original_frame[self.input_mcpeseries_name] = self.mcpe_series_map
        self.PushFrame(self.original_frame)

# class I3FrameSplitterClean(icetray.I3Co):
    
class I3FrameEnergySplitter(icetray.I3ConditionalModule):
    """
    I3Module that breaks apart frames into smaller frames. 
    Needed for ultra high energy events to decrease the memory usage
    
    """
    def __init__(self, context):
        """
        
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
        
        mctree = frame[self.input_mctree_name]
        self.PushFrame(I3FrameStreamChanger.ChangeFrame(frame, self.dummy_stream))
        energy_counter = 0.
        new_frame = icetray.I3Frame(icetray.I3Frame.DAQ)
        new_mctree = dataclasses.I3MCTree()
        dummy_primary = dataclasses.I3Particle(dataclasses.I3Particle.ParticleShape.Primary)
        new_mctree.add_primary(dummy_primary)
        for p in mctree:
            if p.shape == dataclasses.I3Particle.Dark or\
               p.location_type != dataclasses.I3Particle.InIce:
               continue
            new_mctree.append_child(dummy_primary, p)
            if p.type == dataclasses.I3Particle.MuMinus or \
               p.type == dataclasses.I3Particle.MuPlus:
                   energy_counter += 45.*icetray.I3Units.GeV
            else:
                energy_counter += p.energy
            if energy_counter >= self.energy_per_frame:
                print energy_counter
                new_frame[self.input_mctree_name] = new_mctree
                self.PushFrame(new_frame)
                new_frame = icetray.I3Frame(icetray.I3Frame.DAQ)
                new_mctree = dataclasses.I3MCTree()
                new_mctree.add_primary(dummy_primary)
                energy_counter = 0.
        new_frame[self.input_mctree_name] = new_mctree
        self.PushFrame(new_frame)
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
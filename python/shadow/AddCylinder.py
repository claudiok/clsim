from icecube import icetray , dataclasses , clsim , simclasses
import os
from icecube.icetray import OMKey
from icecube.simclasses import I3CylinderMap
import numpy as np

source_directory = os.environ['I3_SRC']

file = open(source_directory + "/ppc/resources/ice/dx.dat")
contents = []
for line in file:
    contents += [line.split()]

DOMs = []

for key in contents:
    DOMs.append( [OMKey( int(key[0]) , int(key[1]) ) , float(key[2])] ) # i contains string , OM_number , orientation of cable , and error on orientation

class AddCylinder(icetray.I3Module):

    def __init__(self,context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter( "TopCylinder" , "Position of top of cylinder" , dataclasses.I3Position( 0 , 0 , 500 ))
        self.AddParameter( "BottomCylinder" , "Position of the bottom of cylinder" , dataclasses.I3Position( 0 , 0 , -500 ))
        self.AddParameter( "Radius" , "Radius of Cylinder" , 10)
        
    def Configure(self):
        self.top = self.GetParameter('TopCylinder')
        self.bottom = self.GetParameter('BottomCylinder')
        self.radius = self.GetParameter('Radius')
        

    def Geometry(self,frame):
        frame['Cable'] = simclasses.I3ExtraGeometryItemCylinder(self.top , self.bottom , self.radius)
        self.PushFrame(frame)

            
class AddCylinders(icetray.I3Module):
    def __init__(self,context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter( "Length_of_cylinder", "Length of the cable" , 1.0 )
        self.AddParameter( "Radius_of_cylinder" , "Radius of the cable" , 0.023 ) #The radius of the cables is 23 mm
        self.AddParameter( "Cable_map" , "Map of cables in geometry" , " " )
        self.AddParameter( "Radius_of_DOM" , "Radius of the DOM" , 0.5 )
        
    def Configure(self):
        self.height = self.GetParameter("Length_of_cylinder")
        self.radius = self.GetParameter("Radius_of_cylinder")
        self.cable_map = self.GetParameter("Cable_map")
        self.dom_radius = self.GetParameter("Radius_of_DOM")
       
    def Geometry(self,frame):
        geometry = frame["I3Geometry"]
        for i in DOMs:
            om_key = i[0]
            orientation = i[1] * (np.pi/180.0)
            position_x = geometry.omgeo[ om_key ].position.x + self.dom_radius + self.radius * np.cos( orientation )
            position_y = geometry.omgeo[ om_key ].position.y + 0.5 + self.dom_radius * np.sin( orientation )
            position_z = geometry.omgeo[ om_key ].position.z
            self.cable_map[om_key] = simclasses.I3ExtraGeometryItemCylinder(dataclasses.I3Position( position_x , position_y , position_z + self.height/2.0),
                                                                       dataclasses.I3Position( position_x , position_y , position_z - self.height/2.0),
                                                                       self.radius)

        frame['CableMap'] = self.cable_map
        self.PushFrame(frame)


        

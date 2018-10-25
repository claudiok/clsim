from icecube import icetray , dataclasses , clsim , simclasses
import os
from icecube.icetray import OMKey
from icecube.simclasses import I3CylinderMap

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
        frame['Cable'] = clsim.I3ExtraGeometryItemCylinder(self.top , self.bottom , self.radius)
        self.PushFrame(frame)

            
class AddCylinders(icetray.I3Module):
    def __init__(self,context):
        icetray.I3Module.__init__(self,context)
        self.AddParameter( "Length_of_cylinder", "Length of the cable" , 1.0 )
        self.AddParameter( "Radius_of_cylinder" , "Radius of the cable" , 0.023 ) #The radius of the cables is 23 mm
        
    def Configure(self):
        height = self.GetParameter("Length_of_cylinder")
        radius = self.GetParameter("Radius_of_cylinder")
       
    def Geometry(self,frame):
        geometry = frame["I3Geometry"]
        cylinderMap = icetray.I3CylinderMap
        for i in DOMs:
            om_key = i[0]
            orientation = i[1]
            position_x = geometry[ om_key ].position.GetX() + 0.5 + radius * cos( orientation )
            position_y = geometry[ om_key ].position.GetY() + 0.5 + radius * sin( orientation )
            position_z = geometry[ om_key ].position.GetZ()
            cylinderMap[om_key] = clsim.I3ExtraGeometryItemCylinder(dataclasses.I3Position( position_x , position_y , position_z + height/2.0 ),
                                                               dataclasses.I3Position( position_x , position_y , position_z - height/2.0 ),
                                                               radius)
            frame["cable_map"] = cylinderMap

        self.PushFrame(frame)


        

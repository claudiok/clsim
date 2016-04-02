from __future__ import print_function
from icecube import icetray, dataclasses
import math

scenario = window.gl.scenario
scenario.clear()

app.files.openFile('test_cascades_clsim.i3')
app.files.selectFrameByIdx( 4 )
pivot_pos = frame["I3MCTree"].most_energetic_cascade.pos

camera_x, camera_y, camera_z = (276.944, 542.8, 593.621)

window.gl.setCameraPivot(pivot_pos.x,pivot_pos.y,pivot_pos.z)
window.gl.setCameraLoc( camera_x, camera_y, camera_z )

dist = math.sqrt( (pivot_pos.x-camera_x)**2 + (pivot_pos.y-camera_y)**2 + (pivot_pos.z-camera_z)**2 )
print("dist camera-pivot:", dist)

window.timeline.setEventTimeWindow(0,2000)
window.timeline.minTime=0
window.timeline.maxTime=2000
window.timeline.time=0


geometry_artist = scenario.add( 'Detector Geometry', ['I3Geometry'] )
scenario.changeSetting( geometry_artist, 'linewidth', 4 )
paths_artist = scenario.add( 'Photon paths', ['PropagatedPhotons','I3MCTree_preMuonProp'] )
scenario.changeSetting( paths_artist, 'antares_ref_index', False )
scenario.changeSetting( paths_artist, 'linewidth', 4 )
scenario.changeSetting( paths_artist, 'min_residual', 0 )
scenario.changeSetting( paths_artist, 'max_residual', 200 )
scenario.add( 'Ice' )


# from icecube.steamshovel.util import movie

# oversample = 3
# width = 1920
# height = 1080

# M1 = movie.Movie(
# 	window,
# 	size=(width*oversample,height*oversample),
# 	outputsize=(width,height),
# 	rotation = -.1,
# 	zoom=0.,
# 	nframes=1000
# )
# M1.create('cascade_1TeV_in_ice.mp4', keepfiles=False)


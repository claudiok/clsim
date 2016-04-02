import math

scenario = window.gl.scenario
scenario.clear()


app.files.openFile('test_flashes_clsim.i3')
app.files.selectFrameByIdx( 4 )

flashing_om = frame["I3FlasherInfo_OMKeys"][0]
flasher_pos = frame["I3Geometry"].omgeo[flashing_om].position

dist_camera_pivot = 483.655638114
camera_x, camera_y, camera_z = ( 581.148, 558.933, 442.08 )

dx = flasher_pos.x-camera_x
dy = flasher_pos.y-camera_y
dz = flasher_pos.z-camera_z
dd = math.sqrt( dx**2 + dy**2 + dz**2 ) 

camera_x = flasher_pos.x-dx*dist_camera_pivot/dd
camera_y = flasher_pos.y-dy*dist_camera_pivot/dd
camera_z = flasher_pos.z-dz*dist_camera_pivot/dd

window.gl.setCameraPivot(flasher_pos.x,flasher_pos.y,flasher_pos.z)
window.gl.setCameraLoc( camera_x, camera_y, camera_z )

window.timeline.setEventTimeWindow(0,2000)
window.timeline.minTime=0
window.timeline.maxTime=2000
window.timeline.time=0


scenario.add( 'Detector Geometry', ['I3Geometry'] )

paths_artist = scenario.add( 'Photon paths for flashers', ['PropagatedPhotons','I3FlasherInfo_pulses'] )
scenario.changeSetting( paths_artist, 'antares_ref_index', False )
# scenario.changeSetting( paths_artist, 'linewidth', 4 )
scenario.changeSetting( paths_artist, 'min_residual', 0 )
scenario.changeSetting( paths_artist, 'max_residual', 200 )

scenario.add( 'Ice' )


from icecube.steamshovel.util import movie

# M1 = movie.Movie( window, size=(1920,1080), rotation = -.1, nframes=1000 )
# M1.create('track_100TeV_in_ice.mp4', keepfiles=False)


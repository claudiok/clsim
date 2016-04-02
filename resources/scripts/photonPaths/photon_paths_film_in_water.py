app.files.openFile('test_muons_clsim_ANTARES.i3')

scenario = window.gl.scenario
scenario.clear()
scenario.add( 'Detector Geometry', ['I3Geometry'] )
paths_artist = scenario.add( 'Photon paths', ['PropagatedPhotons','I3MCTree_preMuonProp'] )
scenario.changeSetting( paths_artist, 'antares_ref_index', True )
scenario.add( 'Ice' )

from icecube.steamshovel.util import movie

window.gl.setCameraPivot(0,0,0)

app.files.selectFrameByIdx( 4 )
window.timeline.maxTime=10000
window.timeline.minTime=0
window.gl.setCameraLoc( 872, 1490, 168 )

M1 = movie.Movie( window, size=(1920,1080), rotation = -.1, nframes=1000 )
M1.create('track_100TeV_in_water.mp4', keepfiles=False)


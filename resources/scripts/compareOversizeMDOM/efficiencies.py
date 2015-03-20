
from icecube.icetray import I3Units
from icecube import clsim
import numpy

glassThickness = 14.*I3Units.mm
glassAbslen = numpy.array(
               [ 0.00*I3Units.cm,  # 260nm
                 0.00*I3Units.cm,  # 270nm
                 0.00*I3Units.cm,  # 280nm
                 0.00*I3Units.cm,  # 290nm
                 
                 0.17*I3Units.cm,  # 300nm
                 0.39*I3Units.cm,  # 310nm
                 0.84*I3Units.cm,  # 320nm
                 1.82*I3Units.cm,  # 330nm
                 3.92*I3Units.cm,  # 340nm
                 8.41*I3Units.cm,  # 350nm
                18.09*I3Units.cm,  # 360nm
                27.21*I3Units.cm,  # 370nm
                19.23*I3Units.cm,  # 380nm
                61.84*I3Units.cm,  # 390nm
               128.04*I3Units.cm,  # 400nm
                81.25*I3Units.cm,  # 410nm
                73.02*I3Units.cm,  # 420nm
                77.30*I3Units.cm,  # 430nm
                65.66*I3Units.cm,  # 440nm
                81.63*I3Units.cm,  # 450nm
               109.23*I3Units.cm,  # 460nm
               116.08*I3Units.cm,  # 470nm
               113.90*I3Units.cm,  # 480nm
               118.86*I3Units.cm,  # 490nm
               126.55*I3Units.cm,  # 500nm
               139.70*I3Units.cm,  # 510nm
               145.68*I3Units.cm,  # 520nm
               150.88*I3Units.cm,  # 530nm
               151.80*I3Units.cm,  # 540nm
               147.16*I3Units.cm,  # 550nm
               142.40*I3Units.cm,  # 560nm
               138.27*I3Units.cm,  # 570nm
               134.58*I3Units.cm,  # 580nm
               135.64*I3Units.cm,  # 590nm
               142.87*I3Units.cm,  # 600nm
               148.37*I3Units.cm,  # 610nm
               
               148.37*I3Units.cm,  # 620nm
               148.37*I3Units.cm,  # 630nm
               148.37*I3Units.cm,  # 640nm
               148.37*I3Units.cm,  # 650nm
               148.37*I3Units.cm,  # 660nm
               148.37*I3Units.cm,  # 670nm
               148.37*I3Units.cm,  # 680nm
               148.37*I3Units.cm,  # 690nm
               
               ])
glassAbsLen = clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, glassAbslen)

gelAbslen = numpy.array(
            [  0.00*I3Units.cm, # 260nm
               0.00*I3Units.cm, # 270nm
               0.00*I3Units.cm, # 280nm
               0.00*I3Units.cm, # 290nm
               
               0.00*I3Units.cm, # 300nm
               8.00*I3Units.cm, # 310nm
              15.60*I3Units.cm, # 320nm
              23.08*I3Units.cm, # 330nm
              30.49*I3Units.cm, # 340nm
              37.14*I3Units.cm, # 350nm
              41.88*I3Units.cm, # 360nm
              45.71*I3Units.cm, # 370nm
              48.96*I3Units.cm, # 380nm
              53.29*I3Units.cm, # 390nm
              56.64*I3Units.cm, # 400nm
              59.38*I3Units.cm, # 410nm
              62.53*I3Units.cm, # 420nm
              64.48*I3Units.cm, # 430nm
              66.91*I3Units.cm, # 440nm
              68.05*I3Units.cm, # 450nm
              72.31*I3Units.cm, # 460nm
              74.55*I3Units.cm, # 470nm
              76.48*I3Units.cm, # 480nm
              78.18*I3Units.cm, # 490nm
              81.08*I3Units.cm, # 500nm
              84.49*I3Units.cm, # 510nm
              85.88*I3Units.cm, # 520nm
              86.95*I3Units.cm, # 530nm
              90.10*I3Units.cm, # 540nm
              89.09*I3Units.cm, # 550nm
              94.36*I3Units.cm, # 560nm
              96.42*I3Units.cm, # 570nm
              96.90*I3Units.cm, # 580nm
              99.89*I3Units.cm, # 590nm
              99.94*I3Units.cm, # 600nm
             100.81*I3Units.cm, # 610nm

             100.81*I3Units.cm, # 620nm
             100.81*I3Units.cm, # 630nm
             100.81*I3Units.cm, # 640nm
             100.81*I3Units.cm, # 650nm
             100.81*I3Units.cm, # 660nm
             100.81*I3Units.cm, # 670nm
             100.81*I3Units.cm, # 680nm
             100.81*I3Units.cm, # 690nm
             ])
gelAbsLen = clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, gelAbslen)

wavelengthacceptance = numpy.array(
                       [ 0.0 * 0.01, # 260nm
                         0.0 * 0.01, # 270nm
                         0.5 * 0.01, # 280nm
                         3.1 * 0.01, # 290nm
                         9.8 * 0.01, # 300nm
                        17.5 * 0.01, # 310nm
                        23.2 * 0.01, # 320nm
                        26.5 * 0.01, # 330nm
                        28.1 * 0.01, # 340nm
                        28.1 * 0.01, # 350nm
                        29.1 * 0.01, # 360nm
                        30.1 * 0.01, # 370nm
                        30.4 * 0.01, # 380nm
                        30.1 * 0.01, # 390nm
                        29.9 * 0.01, # 400nm
                        29.3 * 0.01, # 410nm
                        28.6 * 0.01, # 420nm
                        27.5 * 0.01, # 430nm
                        26.5 * 0.01, # 440nm
                        25.0 * 0.01, # 450nm
                        23.2 * 0.01, # 460nm
                        21.1 * 0.01, # 470nm
                        19.6 * 0.01, # 480nm
                        18.5 * 0.01, # 490nm
                        17.2 * 0.01, # 500nm
                        15.4 * 0.01, # 510nm
                        12.1 * 0.01, # 520nm
                         9.3 * 0.01, # 530nm
                         7.2 * 0.01, # 540nm
                         6.2 * 0.01, # 550nm
                         4.6 * 0.01, # 560nm
                         3.6 * 0.01, # 570nm
                         2.8 * 0.01, # 580nm
                         2.1 * 0.01, # 590nm
                         1.3 * 0.01, # 600nm
                         0.8 * 0.01, # 610nm
                         0.5 * 0.01, # 620nm
                         0.3 * 0.01, # 630nm
                         0.0 * 0.01, # 640nm
                         0.0 * 0.01, # 650nm
                         0.0 * 0.01, # 660nm
                         0.0 * 0.01, # 670nm
                         0.0 * 0.01, # 680nm
                         0.0 * 0.01, # 690nm
                         ])
wavelengthAcceptance = clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, wavelengthacceptance)

cosines = [-1.,   -0.95, -0.9,  -0.85, -0.8,  -0.75, -0.7,  -0.65, -0.6 , -0.55, -0.5 , -0.45,
           -0.4,  -0.35, -0.3,  -0.25, -0.2,  -0.15, -0.1,  -0.05,  0.,    0.05]
acceptanceRatios = [ 1.56901889,   1.23410288,   1.17370277,   1.13023417,   1.1384318,
                     1.12859394,   1.09788943,   1.10675797,   1.11971921,   1.11214495,
                     1.15454716,   1.15321001,   1.2245022,    1.27309651,   1.30277673,
                     1.40193728,   1.52872749,   1.92281438,   1.94452437,   1.91083794,
                     1.92,         1.92] # regularize the acceptance ratios here..
                     #2.91631343,  29.62337522]
winstonRatioFunc = clsim.util.interpolate.interp1d(cosines,acceptanceRatios)

cos_bins = numpy.linspace(-1.,1.,1001)
angularAcceptance = []
for cos_bin in cos_bins:
    if cos_bin <= 0.:
        angularAcceptance.append(0.)
    else:
        angularAcceptance.append(cos_bin * winstonRatioFunc(-cos_bin))
angularAcceptance = clsim.I3CLSimFunctionFromTable(cos_bins[0], cos_bins[1]-cos_bins[0], angularAcceptance)

dom2007a_eff_area = [
    0.0000064522,
    0.0000064522,
    0.0000064522,
    0.0000064522,
    0.0000021980,
    0.0001339040,
    0.0005556810,
    0.0016953000,
    0.0035997000,
    0.0061340900,
    0.0074592700,
    0.0090579800,
    0.0099246700,
    0.0105769000,
    0.0110961000,
    0.0114214000,
    0.0114425000,
    0.0111527000,
    0.0108086000,
    0.0104458000,
    0.0099763100,
    0.0093102500,
    0.0087516600,
    0.0083225800,
    0.0079767200,
    0.0075625100,
    0.0066377000,
    0.0053335800,
    0.0043789400,
    0.0037583500,
    0.0033279800,
    0.0029212500,
    0.0025334900,
    0.0021115400,
    0.0017363300,
    0.0013552700,
    0.0010546600,
    0.0007201020,
    0.0004843820,
    0.0002911110,
    0.0001782310,
    0.0001144300,
    0.0000509155,
    0]
dom2007a_eff_area = numpy.array(dom2007a_eff_area)*I3Units.meter2 # apply units (this is an effective area)
domRadius = 0.16510*I3Units.m
domArea = numpy.pi*domRadius**2.
dom2007a_efficiency = (dom2007a_eff_area/domArea)

mdom_efficiency = wavelengthacceptance*numpy.exp(-glassThickness/glassAbslen)*numpy.exp(-1*I3Units.cm/gelAbslen)
envelope_efficiency = numpy.maximum(dom2007a_efficiency, mdom_efficiency)

def GetIceCubeDOMAcceptance(oversize=1., efficiency=1.):
	return clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, efficiency*dom2007a_efficiency/oversize**2)
def GetMDOMAcceptance(oversize=1., efficiency=1.):
	return clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, efficiency*mdom_efficiency/oversize**2)
def GetAcceptanceEnvelope(oversize=1., efficiency=1.):
	"""
	Get the envelope of the maximal DOM and mDOM acceptances
	"""
	return clsim.I3CLSimFunctionFromTable(260.*I3Units.nanometer, 10.*I3Units.nanometer, efficiency*envelope_efficiency/oversize**2)

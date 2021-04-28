# Batch processing script.
# By default ParaView uses nearest neighbor interpolation.
#
# Hongyang Zhou, hyzhou@umich.edu 11/18/2019

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
testvtu = XMLUnstructuredGridReader(
    FileName=['../../data/test.vtu'])
testvtu.PointArrayStatus = ['B [nT]', 'U [km/s]', 'Rho [amu/cm^3]',
    'J [`mA/m^2]']

# create a new 'CSV Reader'
boundarytxt = CSVReader(
    FileName=['data/surface.txt'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=boundarytxt)
tableToPoints1.XColumn = 'X'
tableToPoints1.YColumn = 'Y'
tableToPoints1.ZColumn = 'Z'

# create a new 'Point Dataset Interpolator'
pointDatasetInterpolator1 = PointDatasetInterpolator(Input=testvtu,
    Source=tableToPoints1)
#pointDatasetInterpolator1.Kernel = 'VoronoiKernel'
pointDatasetInterpolator1.Kernel = 'LinearKernel'
pointDatasetInterpolator1.Kernel.Radius = 0.0625
pointDatasetInterpolator1.Locator = 'Static Point Locator'

# save data
SaveData('data/surface_value.txt',
    proxy=pointDatasetInterpolator1)

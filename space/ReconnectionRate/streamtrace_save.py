# Batch processing script.
#
# Hongyang Zhou, hyzhou@umich.edu 11/18/2019

#### import the simple module from the paraview
from paraview.simple import *

version = paraview.servermanager.vtkSMProxyManager.GetVersionMajor()*10 + \
          paraview.servermanager.vtkSMProxyManager.GetVersionMinor()

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
data1 = XMLUnstructuredGridReader(
	FileName=['../../data/test.vtu'])

# Properties modified on data1
data1.PointArrayStatus = ['B [nT]']

# create a new 'Clip'
clip1 = Clip(Input=data1)
clip1.ClipType = 'Box'
clip1.ClipType.Position = [-2.0, -3.0, 0.0]
if version >= 57:
    clip1.ClipType.Length = [6.5, 6.0, 2.1]
else: # On Paraview 5.6 or older
    clip1.ClipType.Scale = [6.5, 6.0, 2.1]

# create a new 'CSV Reader'
seeds = CSVReader(
	FileName=['data/seeds.txt'])

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=seeds)
tableToPoints1.XColumn = 'X'
tableToPoints1.YColumn = 'Y'
tableToPoints1.ZColumn = 'Z'

# create a new Stream Tracer
streamTracer = StreamTracerWithCustomSource(Input=clip1,
    SeedSource=tableToPoints1)
streamTracer.Vectors = ['POINTS', 'B [nT]']
streamTracer.MaximumStreamlineLength = 6.0

# Properties modified on streamTracer
streamTracer.IntegrationDirection = 'FORWARD'
streamTracer.ComputeVorticity = 0

# save data
SaveData('data/streamline.txt',
	proxy=streamTracer)

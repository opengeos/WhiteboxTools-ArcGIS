import arcpy
from WBT.whitebox_tools import WhiteboxTools
wbt = WhiteboxTools()

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        tools = []        
        tools.append(Help)
        tools.append(License)
        tools.append(Version)
        tools.append(ListTools)
        tools.append(ToolHelp)
        tools.append(ToolParameters)
        tools.append(ViewCode)
        tools.append(RunTool)
        tools.append(Update)

        tools.append(AddPointCoordinatesToTable)
        tools.append(ConvertNodataToZero)
        tools.append(ConvertRasterFormat)
        tools.append(ExportTableToCsv)
        tools.append(JoinTables)
        tools.append(LinesToPolygons)
        tools.append(MergeTableWithCsv)
        tools.append(MergeVectors)
        tools.append(MultiPartToSinglePart)
        tools.append(NewRasterFromBase)
        tools.append(PolygonsToLines)
        tools.append(PrintGeoTiffTags)
        tools.append(RasterToVectorLines)
        tools.append(RasterToVectorPoints)
        tools.append(ReinitializeAttributeTable)
        tools.append(RemovePolygonHoles)
        tools.append(SetNodataValue)
        tools.append(SinglePartToMultiPart)
        tools.append(VectorLinesToRaster)
        tools.append(VectorPointsToRaster)
        tools.append(VectorPolygonsToRaster)
        tools.append(AggregateRaster)
        tools.append(BlockMaximumGridding)
        tools.append(BlockMinimumGridding)
        tools.append(Centroid)
        tools.append(CentroidVector)
        tools.append(Clump)
        tools.append(ConstructVectorTin)
        tools.append(CreateHexagonalVectorGrid)
        tools.append(CreatePlane)
        tools.append(CreateRectangularVectorGrid)
        tools.append(Dissolve)
        tools.append(EliminateCoincidentPoints)
        tools.append(ExtendVectorLines)
        tools.append(ExtractNodes)
        tools.append(ExtractRasterValuesAtPoints)
        tools.append(FindLowestOrHighestPoints)
        tools.append(IdwInterpolation)
        tools.append(LayerFootprint)
        tools.append(Medoid)
        tools.append(MinimumBoundingBox)
        tools.append(MinimumBoundingCircle)
        tools.append(MinimumBoundingEnvelope)
        tools.append(MinimumConvexHull)
        tools.append(NearestNeighbourGridding)
        tools.append(PolygonArea)
        tools.append(PolygonLongAxis)
        tools.append(PolygonPerimeter)
        tools.append(PolygonShortAxis)
        tools.append(RasterCellAssignment)
        tools.append(Reclass)
        tools.append(ReclassEqualInterval)
        tools.append(ReclassFromFile)
        tools.append(SmoothVectors)
        tools.append(TinGridding)
        tools.append(VectorHexBinning)
        tools.append(VoronoiDiagram)
        tools.append(BufferRaster)
        tools.append(CostAllocation)
        tools.append(CostDistance)
        tools.append(CostPathway)
        tools.append(EuclideanAllocation)
        tools.append(EuclideanDistance)
        tools.append(AverageOverlay)
        tools.append(Clip)
        tools.append(ClipRasterToPolygon)
        tools.append(CountIf)
        tools.append(Difference)
        tools.append(Erase)
        tools.append(ErasePolygonFromRaster)
        tools.append(HighestPosition)
        tools.append(Intersect)
        tools.append(LineIntersections)
        tools.append(LowestPosition)
        tools.append(MaxAbsoluteOverlay)
        tools.append(MaxOverlay)
        tools.append(MinAbsoluteOverlay)
        tools.append(MinOverlay)
        tools.append(PercentEqualTo)
        tools.append(PercentGreaterThan)
        tools.append(PercentLessThan)
        tools.append(PickFromList)
        tools.append(Polygonize)
        tools.append(SplitWithLines)
        tools.append(SumOverlay)
        tools.append(SymmetricalDifference)
        tools.append(Union)
        tools.append(WeightedOverlay)
        tools.append(WeightedSum)
        tools.append(CompactnessRatio)
        tools.append(EdgeProportion)
        tools.append(ElongationRatio)
        tools.append(FindPatchOrClassEdgeCells)
        tools.append(HoleProportion)
        tools.append(LinearityIndex)
        tools.append(PatchOrientation)
        tools.append(PerimeterAreaRatio)
        tools.append(RadiusOfGyration)
        tools.append(RelatedCircumscribingCircle)
        tools.append(ShapeComplexityIndex)
        tools.append(Aspect)
        tools.append(CircularVarianceOfAspect)
        tools.append(DevFromMeanElev)
        tools.append(DiffFromMeanElev)
        tools.append(DirectionalRelief)
        tools.append(DownslopeIndex)
        tools.append(DrainagePreservingSmoothing)
        tools.append(EdgeDensity)
        tools.append(ElevAbovePit)
        tools.append(ElevPercentile)
        tools.append(ElevRelativeToMinMax)
        tools.append(ElevRelativeToWatershedMinMax)
        tools.append(FeaturePreservingDenoise)
        tools.append(FetchAnalysis)
        tools.append(FillMissingData)
        tools.append(FindRidges)
        tools.append(Hillshade)
        tools.append(HorizonAngle)
        tools.append(HypsometricAnalysis)
        tools.append(MaxAnisotropyDev)
        tools.append(MaxAnisotropyDevSignature)
        tools.append(MaxBranchLength)
        tools.append(MaxDifferenceFromMean)
        tools.append(MaxDownslopeElevChange)
        tools.append(MaxElevDevSignature)
        tools.append(MaxElevationDeviation)
        tools.append(MinDownslopeElevChange)
        tools.append(MultiscaleRoughness)
        tools.append(MultiscaleRoughnessSignature)
        tools.append(MultiscaleTopographicPositionImage)
        tools.append(NumDownslopeNeighbours)
        tools.append(NumUpslopeNeighbours)
        tools.append(PennockLandformClass)
        tools.append(PercentElevRange)
        tools.append(PlanCurvature)
        tools.append(Profile)
        tools.append(ProfileCurvature)
        tools.append(RelativeAspect)
        tools.append(RelativeStreamPowerIndex)
        tools.append(RelativeTopographicPosition)
        tools.append(RemoveOffTerrainObjects)
        tools.append(RuggednessIndex)
        tools.append(SedimentTransportIndex)
        tools.append(Slope)
        tools.append(SlopeVsElevationPlot)
        tools.append(StandardDeviationOfSlope)
        tools.append(SurfaceAreaRatio)
        tools.append(TangentialCurvature)
        tools.append(TotalCurvature)
        tools.append(Viewshed)
        tools.append(VisibilityIndex)
        tools.append(WetnessIndex)
        tools.append(AverageFlowpathSlope)
        tools.append(AverageUpslopeFlowpathLength)
        tools.append(Basins)
        tools.append(BreachDepressions)
        tools.append(BreachSingleCellPits)
        tools.append(D8FlowAccumulation)
        tools.append(D8MassFlux)
        tools.append(D8Pointer)
        tools.append(DInfFlowAccumulation)
        tools.append(DInfMassFlux)
        tools.append(DInfPointer)
        tools.append(DepthInSink)
        tools.append(DownslopeDistanceToStream)
        tools.append(DownslopeFlowpathLength)
        tools.append(ElevationAboveStream)
        tools.append(ElevationAboveStreamEuclidean)
        tools.append(Fd8FlowAccumulation)
        tools.append(Fd8Pointer)
        tools.append(FillBurn)
        tools.append(FillDepressions)
        tools.append(FillSingleCellPits)
        tools.append(FindNoFlowCells)
        tools.append(FindParallelFlow)
        tools.append(FlattenLakes)
        tools.append(FloodOrder)
        tools.append(FlowAccumulationFullWorkflow)
        tools.append(FlowLengthDiff)
        tools.append(Hillslopes)
        tools.append(ImpoundmentIndex)
        tools.append(Isobasins)
        tools.append(JensonSnapPourPoints)
        tools.append(LongestFlowpath)
        tools.append(MaxUpslopeFlowpathLength)
        tools.append(NumInflowingNeighbours)
        tools.append(RaiseWalls)
        tools.append(Rho8Pointer)
        tools.append(Sink)
        tools.append(SnapPourPoints)
        tools.append(StochasticDepressionAnalysis)
        tools.append(StrahlerOrderBasins)
        tools.append(Subbasins)
        tools.append(TraceDownslopeFlowpaths)
        tools.append(UnnestBasins)
        tools.append(Watershed)
        tools.append(ChangeVectorAnalysis)
        tools.append(Closing)
        tools.append(CreateColourComposite)
        tools.append(FlipImage)
        tools.append(IhsToRgb)
        tools.append(ImageStackProfile)
        tools.append(IntegralImage)
        tools.append(KMeansClustering)
        tools.append(LineThinning)
        tools.append(ModifiedKMeansClustering)
        tools.append(Mosaic)
        tools.append(MosaicWithFeathering)
        tools.append(NormalizedDifferenceVegetationIndex)
        tools.append(Opening)
        tools.append(RemoveSpurs)
        tools.append(Resample)
        tools.append(RgbToIhs)
        tools.append(SplitColourComposite)
        tools.append(ThickenRasterLine)
        tools.append(TophatTransform)
        tools.append(WriteFunctionMemoryInsertion)
        tools.append(AdaptiveFilter)
        tools.append(BilateralFilter)
        tools.append(ConservativeSmoothingFilter)
        tools.append(CornerDetection)
        tools.append(DiffOfGaussianFilter)
        tools.append(DiversityFilter)
        tools.append(EdgePreservingMeanFilter)
        tools.append(EmbossFilter)
        tools.append(FastAlmostGaussianFilter)
        tools.append(GaussianFilter)
        tools.append(HighPassFilter)
        tools.append(HighPassMedianFilter)
        tools.append(KNearestMeanFilter)
        tools.append(LaplacianFilter)
        tools.append(LaplacianOfGaussianFilter)
        tools.append(LeeFilter)
        tools.append(LineDetectionFilter)
        tools.append(MajorityFilter)
        tools.append(MaximumFilter)
        tools.append(MeanFilter)
        tools.append(MedianFilter)
        tools.append(MinimumFilter)
        tools.append(OlympicFilter)
        tools.append(PercentileFilter)
        tools.append(PrewittFilter)
        tools.append(RangeFilter)
        tools.append(RobertsCrossFilter)
        tools.append(ScharrFilter)
        tools.append(SobelFilter)
        tools.append(StandardDeviationFilter)
        tools.append(TotalFilter)
        tools.append(UnsharpMasking)
        tools.append(UserInedWeightsFilter)
        tools.append(BalanceContrastEnhancement)
        tools.append(CorrectVignetting)
        tools.append(DirectDecorrelationStretch)
        tools.append(GammaCorrection)
        tools.append(GaussianContrastStretch)
        tools.append(HistogramEqualization)
        tools.append(HistogramMatching)
        tools.append(HistogramMatchingTwoImages)
        tools.append(MinMaxContrastStretch)
        tools.append(PanchromaticSharpening)
        tools.append(PercentageContrastStretch)
        tools.append(SigmoidalContrastStretch)
        tools.append(StandardDeviationContrastStretch)
        tools.append(ClassifyOverlapPoints)
        tools.append(ClipLidarToPolygon)
        tools.append(ErasePolygonFromLidar)
        tools.append(FilterLidarScanAngles)
        tools.append(FindFlightlineEdgePoints)
        tools.append(FlightlineOverlap)
        tools.append(LasToAscii)
        tools.append(LasToMultipointShapefile)
        tools.append(LasToShapefile)
        tools.append(LidarBlockMaximum)
        tools.append(LidarBlockMinimum)
        tools.append(LidarClassifySubset)
        tools.append(LidarColourize)
        tools.append(LidarConstructVectorTin)
        tools.append(LidarElevationSlice)
        tools.append(LidarGroundPointFilter)
        tools.append(LidarHexBinning)
        tools.append(LidarHillshade)
        tools.append(LidarHistogram)
        tools.append(LidarIdwInterpolation)
        tools.append(LidarInfo)
        tools.append(LidarJoin)
        tools.append(LidarKappaIndex)
        tools.append(LidarNearestNeighbourGridding)
        tools.append(LidarPointDensity)
        tools.append(LidarPointStats)
        tools.append(LidarRemoveDuplicates)
        tools.append(LidarRemoveOutliers)
        tools.append(LidarSegmentation)
        tools.append(LidarSegmentationBasedFilter)
        tools.append(LidarThin)
        tools.append(LidarThinHighDensity)
        tools.append(LidarTile)
        tools.append(LidarTileFootprint)
        tools.append(LidarTinGridding)
        tools.append(LidarTophatTransform)
        tools.append(NormalVectors)
        tools.append(SelectTilesByPolygon)
        tools.append(And)
        tools.append(Not)
        tools.append(Or)
        tools.append(AbsoluteValue)
        tools.append(Add)
        tools.append(Anova)
        tools.append(ArcCos)
        tools.append(ArcSin)
        tools.append(ArcTan)
        tools.append(Atan2)
        tools.append(AttributeCorrelation)
        tools.append(AttributeHistogram)
        tools.append(AttributeScattergram)
        tools.append(Ceil)
        tools.append(Cos)
        tools.append(Cosh)
        tools.append(CrispnessIndex)
        tools.append(CrossTabulation)
        tools.append(CumulativeDistribution)
        tools.append(Decrement)
        tools.append(Divide)
        tools.append(EqualTo)
        tools.append(Exp)
        tools.append(Exp2)
        tools.append(ExtractRasterStatistics)
        tools.append(Floor)
        tools.append(GreaterThan)
        tools.append(ImageAutocorrelation)
        tools.append(ImageCorrelation)
        tools.append(ImageRegression)
        tools.append(InPlaceAdd)
        tools.append(InPlaceDivide)
        tools.append(InPlaceMultiply)
        tools.append(InPlaceSubtract)
        tools.append(Increment)
        tools.append(IntegerDivision)
        tools.append(IsNoData)
        tools.append(KappaIndex)
        tools.append(KsTestForNormality)
        tools.append(LessThan)
        tools.append(ListUniqueValues)
        tools.append(Ln)
        tools.append(Log10)
        tools.append(Log2)
        tools.append(Max)
        tools.append(Min)
        tools.append(Modulo)
        tools.append(Multiply)
        tools.append(Negate)
        tools.append(NotEqualTo)
        tools.append(Power)
        tools.append(PrincipalComponentAnalysis)
        tools.append(Quantiles)
        tools.append(RandomField)
        tools.append(RandomSample)
        tools.append(RasterHistogram)
        tools.append(RasterSummaryStats)
        tools.append(Reciprocal)
        tools.append(RescaleValueRange)
        tools.append(RootMeanSquareError)
        tools.append(Round)
        tools.append(Sin)
        tools.append(Sinh)
        tools.append(Square)
        tools.append(SquareRoot)
        tools.append(Subtract)
        tools.append(Tan)
        tools.append(Tanh)
        tools.append(ToDegrees)
        tools.append(ToRadians)
        tools.append(TrendSurface)
        tools.append(TrendSurfaceVectorPoints)
        tools.append(Truncate)
        tools.append(TurningBandsSimulation)
        tools.append(Xor)
        tools.append(ZScores)
        tools.append(DistanceToOutlet)
        tools.append(ExtractStreams)
        tools.append(ExtractValleys)
        tools.append(FarthestChannelHead)
        tools.append(FindMainStem)
        tools.append(HackStreamOrder)
        tools.append(HortonStreamOrder)
        tools.append(LengthOfUpstreamChannels)
        tools.append(LongProfile)
        tools.append(LongProfileFromPoints)
        tools.append(RasterStreamsToVector)
        tools.append(RasterizeStreams)
        tools.append(RemoveShortStreams)
        tools.append(ShreveStreamMagnitude)
        tools.append(StrahlerStreamOrder)
        tools.append(StreamLinkClass)
        tools.append(StreamLinkIdentifier)
        tools.append(StreamLinkLength)
        tools.append(StreamLinkSlope)
        tools.append(StreamSlopeContinuous)
        tools.append(TopologicalStreamOrder)
        tools.append(TributaryIdentifier)

        self.tools = tools

class Help(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Help"
        self.description = "Help description for WhiteboxTools"
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        messages.addMessage(wbt.help())
        return


class License(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "License"
        self.description = "License information for WhiteboxTools."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        messages.addMessage(wbt.license())
        return


class Version(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Version"
        self.description = "Version information for WhiteboxTools."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        messages.addMessage(wbt.version())
        return



class ListTools(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "List Tools"
        self.description = "All available tools in WhiteboxTools."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""

            # First parameter
        param0 = arcpy.Parameter(
            displayName="Keywords",
            name="keywords",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        params = [param0]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        param0 = parameters[0].valueAsText

        if param0 is None:
            tools = wbt.list_tools()
        else:
            tools = wbt.list_tools([param0])
            
        for index, tool in enumerate(tools):
            messages.addMessage("{}. {}: {}".format(index, tool, tools[tool]))
        return

## name conflict with ArcGIS Python Toolbox
# class Toolbox(object):
#     def __init__(self):
#         """Define the tool (tool name is the name of the class)."""
#         self.label = "Toolbox"
#         self.description = "The toolbox for a specific tool."
#         self.category = "About WhiteboxTools"

#     def getParameterInfo(self):
#         """Define parameter definitions"""
#         params = None
#         return params

#     def isLicensed(self):
#         """Set whether tool is licensed to execute."""
#         return True

#     def updateParameters(self, parameters):
#         """Modify the values and properties of parameters before internal
#         validation is performed.  This method is called whenever a parameter
#         has been changed."""
#         return

#     def updateMessages(self, parameters):
#         """Modify the messages created by internal validation for each tool
#         parameter.  This method is called after internal validation."""
#         return

#     def execute(self, parameters, messages):
#         """The source code of the tool."""
#         return


class ToolHelp(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Help"
        self.description = "Help description for a specific tool."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ToolParameters(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Parameters"
        self.description = "Tool parameter descriptions for a specific tool."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ViewCode(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "View Code"
        self.description = "Source code for a specific tool."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        # Set a value list of 1, 10 and 100
        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = ["Lidar Info", "And", "Or", "Not"]

        params = [tool_name]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        param0 = parameters[0].valueAsText
        tool_name = param0.replace(" ", "").strip()
        messages.addMessage(wbt.view_code(tool_name))
        return


class RunTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Run Tool"
        self.description = "Runs a tool and specifies tool arguments."
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Update(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Update"
        self.description = "Download the latest version of WhiteboxTools"
        self.category = "About WhiteboxTools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return
class AddPointCoordinatesToTable(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add Point Coordinates To Table"
        self.description = "Modifies the attribute table of a point vector by adding fields containing each point's X and Y coordinates."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ConvertNodataToZero(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert Nodata To Zero"
        self.description = "Converts nodata values in a raster to zero."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ConvertRasterFormat(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert Raster Format"
        self.description = "Converts raster data from one format to another."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExportTableToCsv(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Export Table To Csv"
        self.description = "Exports an attribute table to a CSV text file."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class JoinTables(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Join Tables"
        self.description = "Merge a vector's attribute table with another table based on a common field."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LinesToPolygons(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lines To Polygons"
        self.description = "Converts vector polylines to polygons."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MergeTableWithCsv(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Merge Table With Csv"
        self.description = "Merge a vector's attribute table with a table contained within a CSV text file."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MergeVectors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Merge Vectors"
        self.description = "Combines two or more input vectors of the same ShapeType creating a single, new output vector."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MultiPartToSinglePart(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multi Part To Single Part"
        self.description = "Converts a vector file containing multi-part features into a vector containing only single-part features."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NewRasterFromBase(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "New Raster From Base"
        self.description = "Creates a new raster using a base image."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PolygonsToLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygons To Lines"
        self.description = "Converts vector polygons to polylines."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PrintGeoTiffTags(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Print Geo Tiff Tags"
        self.description = "Prints the tags within a GeoTIFF."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterToVectorLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster To Vector Lines"
        self.description = "Converts a raster lines features into a vector of the POLYLINE shapetype."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterToVectorPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster To Vector Points"
        self.description = "Converts a raster dataset to a vector of the POINT shapetype."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ReinitializeAttributeTable(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reinitialize Attribute Table"
        self.description = "Reinitializes a vector's attribute table deleting all fields but the feature ID (FID)."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RemovePolygonHoles(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Polygon Holes"
        self.description = "Removes holes within the features of a vector polygon file."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SetNodataValue(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Set Nodata Value"
        self.description = "Assign a specified value in an input image to the NoData value."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SinglePartToMultiPart(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Single Part To Multi Part"
        self.description = "Converts a vector file containing multi-part features into a vector containing only single-part features."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class VectorLinesToRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Lines To Raster"
        self.description = "Converts a vector containing polylines into a raster."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class VectorPointsToRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Points To Raster"
        self.description = "Converts a vector containing points into a raster."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class VectorPolygonsToRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Polygons To Raster"
        self.description = "Converts a vector containing polygons into a raster."
        self.category = "Data Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AggregateRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Aggregate Raster"
        self.description = "Aggregates a raster to a lower resolution."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BlockMaximumGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Block Maximum Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using a block maximum scheme."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BlockMinimumGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Block Minimum Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using a block minimum scheme."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Centroid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Centroid"
        self.description = "Calculates the centroid, or average location, of raster polygon objects."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CentroidVector(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Centroid Vector"
        self.description = "Identifes the centroid point of a vector polyline or polygon feature or a group of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Clump(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clump"
        self.description = "Groups cells that form discrete areas, assigning them unique identifiers."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ConstructVectorTin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Construct Vector Tin"
        self.description = "Creates a vector triangular irregular network (TIN) for a set of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CreateHexagonalVectorGrid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Hexagonal Vector Grid"
        self.description = "Creates a hexagonal vector grid."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CreatePlane(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Plane"
        self.description = "Creates a raster image based on the equation for a simple plane."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CreateRectangularVectorGrid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Rectangular Vector Grid"
        self.description = "Creates a rectangular vector grid."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Dissolve(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Dissolve"
        self.description = "Removes the interior, or shared, boundaries within a vector polygon coverage."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EliminateCoincidentPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Eliminate Coincident Points"
        self.description = "Removes any coincident, or nearly coincident, points from a vector points file."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExtendVectorLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extend Vector Lines"
        self.description = "Extends vector lines by a specified distance."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExtractNodes(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Nodes"
        self.description = "Converts vector lines or polygons into vertex points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExtractRasterValuesAtPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Raster Values At Points"
        self.description = "Extracts the values of raster(s) at vector point locations."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindLowestOrHighestPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Lowest Or Highest Points"
        self.description = "Locates the lowest and/or highest valued cells in a raster."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class IdwInterpolation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Idw Interpolation"
        self.description = "Interpolates vector points into a raster surface using an inverse-distance weighted scheme."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LayerFootprint(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Layer Footprint"
        self.description = "Creates a vector polygon footprint of the area covered by a raster grid or vector layer."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Medoid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Medoid"
        self.description = "Calculates the medoid for a series of vector features contained in a shapefile."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinimumBoundingBox(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Bounding Box"
        self.description = "Creates a vector minimum bounding rectangle around vector features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinimumBoundingCircle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Bounding Circle"
        self.description = "Delineates the minimum bounding circle (i.e. smallest enclosing circle) for a group of vectors."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinimumBoundingEnvelope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Bounding Envelope"
        self.description = "Creates a vector axis-aligned minimum bounding rectangle (envelope) around vector features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinimumConvexHull(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Convex Hull"
        self.description = "Creates a vector convex polygon around vector features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NearestNeighbourGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Nearest Neighbour Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using the nearest neighbour."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PolygonArea(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Area"
        self.description = "Calculates the area of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PolygonLongAxis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Long Axis"
        self.description = "This tool can be used to map the long axis of polygon features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PolygonPerimeter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Perimeter"
        self.description = "Calculates the perimeter of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PolygonShortAxis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Short Axis"
        self.description = "This tool can be used to map the short axis of polygon features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterCellAssignment(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Cell Assignment"
        self.description = "Assign row or column number to cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Reclass(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reclass"
        self.description = "Reclassifies the values in a raster image."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ReclassEqualInterval(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reclass Equal Interval"
        self.description = "Reclassifies the values in a raster image based on equal-ranges."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ReclassFromFile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reclass From File"
        self.description = "Reclassifies the values in a raster image using reclass ranges in a text file."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SmoothVectors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Smooth Vectors"
        self.description = "Smooths a vector coverage of either a POLYLINE or POLYGON base ShapeType."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TinGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tin Gridding"
        self.description = "Creates a raster grid based on a triangular irregular network (TIN) fitted to vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class VectorHexBinning(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Hex Binning"
        self.description = "Hex-bins a set of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class VoronoiDiagram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Voronoi Diagram"
        self.description = "Creates a vector Voronoi diagram for a set of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BufferRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Buffer Raster"
        self.description = "Maps a distance-based buffer around each non-background (non-zero/non-nodata) grid cell in an input image."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CostAllocation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cost Allocation"
        self.description = "Identifies the source cell to which each grid cell is connected by a least-cost pathway in a cost-distance analysis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CostDistance(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cost Distance"
        self.description = "Performs cost-distance accumulation on a cost surface and a group of source cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CostPathway(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cost Pathway"
        self.description = "Performs cost-distance pathway analysis using a series of destination grid cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EuclideanAllocation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Euclidean Allocation"
        self.description = "Assigns grid cells in the output raster the value of the nearest target cell in the input image, measured by the Shih and Wu (2004) Euclidean distance transform."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EuclideanDistance(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Euclidean Distance"
        self.description = "Calculates the Shih and Wu (2004) Euclidean distance transform."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AverageOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average Overlay"
        self.description = "Calculates the average for each grid cell from a group of raster images."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Clip(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clip"
        self.description = "Extract all the features, or parts of features, that overlap with the features of the clip vector."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ClipRasterToPolygon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clip Raster To Polygon"
        self.description = "Clips a raster to a vector polygon."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CountIf(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Count If"
        self.description = "Counts the number of occurrences of a specified value in a cell-stack of rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Difference(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Difference"
        self.description = "Outputs the features that occur in one of the two vector inputs but not both, i.e. no overlapping features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Erase(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Erase"
        self.description = "Removes all the features, or parts of features, that overlap with the features of the erase vector polygon."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ErasePolygonFromRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Erase Polygon From Raster"
        self.description = "Erases (cuts out) a vector polygon from a raster."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HighestPosition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Highest Position"
        self.description = "Identifies the stack position of the maximum value within a raster stack on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Intersect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Intersect"
        self.description = "Identifies the parts of features in common between two input vector layers."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LineIntersections(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Line Intersections"
        self.description = "Identifies points where the features of two vector line layers intersect."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LowestPosition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lowest Position"
        self.description = "Identifies the stack position of the minimum value within a raster stack on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxAbsoluteOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Absolute Overlay"
        self.description = "Evaluates the maximum absolute value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Overlay"
        self.description = "Evaluates the maximum value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinAbsoluteOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Absolute Overlay"
        self.description = "Evaluates the minimum absolute value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Overlay"
        self.description = "Evaluates the minimum value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PercentEqualTo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Equal To"
        self.description = "Calculates the percentage of a raster stack that have cell values equal to an input on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PercentGreaterThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Greater Than"
        self.description = "Calculates the percentage of a raster stack that have cell values greather than an input on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PercentLessThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Less Than"
        self.description = "Calculates the percentage of a raster stack that have cell values less than an input on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PickFromList(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Pick From List"
        self.description = "Outputs the value from a raster stack specified by a position raster."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Polygonize(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygonize"
        self.description = "Creates a polygon layer from two or more intersecting line features contained in one or more input vector line files."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SplitWithLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Split With Lines"
        self.description = "Splits the lines or polygons in one layer using the lines in another layer."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SumOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sum Overlay"
        self.description = "Calculates the sum for each grid cell from a group of raster images."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SymmetricalDifference(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Symmetrical Difference"
        self.description = "Outputs the features that occur in one of the two vector inputs but not both, i.e. no overlapping features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Union(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Union"
        self.description = "Splits vector layers at their overlaps, creating a layer containing all the portions from both input and overlay layers."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class WeightedOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Weighted Overlay"
        self.description = "Performs a weighted sum on multiple input rasters after converting each image to a common scale. The tool performs a multi-criteria evaluation (MCE)."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class WeightedSum(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Weighted Sum"
        self.description = "Performs a weighted-sum overlay on multiple input raster images."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CompactnessRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Compactness Ratio"
        self.description = "Calculates the compactness ratio (A/P), a measure of shape complexity, for vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EdgeProportion(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Proportion"
        self.description = "Calculate the proportion of cells in a raster polygon that are edge cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElongationRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elongation Ratio"
        self.description = "Calculates the elongation ratio for vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindPatchOrClassEdgeCells(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Patch Or Class Edge Cells"
        self.description = "Finds all cells located on the edge of patch or class features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HoleProportion(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hole Proportion"
        self.description = "Calculates the proportion of the total area of a polygon's holes relative to the area of the polygon's hull."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LinearityIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Linearity Index"
        self.description = "Calculates the linearity index for vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PatchOrientation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Patch Orientation"
        self.description = "Calculates the orientation of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PerimeterAreaRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Perimeter Area Ratio"
        self.description = "Calculates the perimeter-area ratio of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RadiusOfGyration(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Radius Of Gyration"
        self.description = "Calculates the distance of cells from their polygon's centroid."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RelatedCircumscribingCircle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Related Circumscribing Circle"
        self.description = "Calculates the related circumscribing circle of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ShapeComplexityIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shape Complexity Index"
        self.description = "Calculates overall polygon shape complexity or irregularity."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Aspect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Aspect"
        self.description = "Calculates an aspect raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CircularVarianceOfAspect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Circular Variance Of Aspect"
        self.description = "Calculates the circular variance of aspect at a scale for a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DevFromMeanElev(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Dev From Mean Elev"
        self.description = "Calculates deviation from mean elevation."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DiffFromMeanElev(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Diff From Mean Elev"
        self.description = "Calculates difference from mean elevation (equivalent to a high-pass filter)."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DirectionalRelief(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Directional Relief"
        self.description = "Calculates relief for cells in an input DEM for a specified direction."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DownslopeIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Downslope Index"
        self.description = "Calculates the Hjerdt et al. (2004) downslope index."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DrainagePreservingSmoothing(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Drainage Preserving Smoothing"
        self.description = "Reduces short-scale variation in an input DEM while preserving breaks-in-slope and small drainage features using a modified Sun et al. (2007) algorithm."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EdgeDensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Density"
        self.description = "Calculates the density of edges, or breaks-in-slope within DEMs."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElevAbovePit(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Above Pit"
        self.description = "Calculate the elevation of each grid cell above the nearest downstream pit cell or grid edge cell."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElevPercentile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Percentile"
        self.description = "Calculates the elevation percentile raster from a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElevRelativeToMinMax(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Relative To Min Max"
        self.description = "Calculates the elevation of a location relative to the minimum and maximum elevations in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElevRelativeToWatershedMinMax(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Relative To Watershed Min Max"
        self.description = "Calculates the elevation of a location relative to the minimum and maximum elevations in a watershed."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FeaturePreservingDenoise(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Feature Preserving Denoise"
        self.description = "Reduces short-scale variation in an input DEM using a modified Sun et al. (2007) algorithm."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FetchAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fetch Analysis"
        self.description = "Performs an analysis of fetch or upwind distance to an obstacle."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FillMissingData(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Missing Data"
        self.description = "Fills NoData holes in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindRidges(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Ridges"
        self.description = "Identifies potential ridge and peak grid cells."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Hillshade(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hillshade"
        self.description = "Calculates a hillshade raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HorizonAngle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Horizon Angle"
        self.description = "Calculates horizon angle (maximum upwind slope) for each grid cell in an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HypsometricAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hypsometric Analysis"
        self.description = "Calculates a hypsometric curve for one or more DEMs."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxAnisotropyDev(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Anisotropy Dev"
        self.description = "Calculates the maximum anisotropy (directionality) in elevation deviation over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxAnisotropyDevSignature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Anisotropy Dev Signature"
        self.description = "Calculates the anisotropy in deviation from mean for points over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxBranchLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Branch Length"
        self.description = "Lindsay and Seibert's (2013) branch length index is used to map drainage divides or ridge lines."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxDifferenceFromMean(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Difference From Mean"
        self.description = "Calculates the maximum difference from mean elevation over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxDownslopeElevChange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Downslope Elev Change"
        self.description = "Calculates the maximum downslope change in elevation between a grid cell and its eight downslope neighbors."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxElevDevSignature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Elev Dev Signature"
        self.description = "Calculates the maximum elevation deviation over a range of spatial scales and for a set of points."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxElevationDeviation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Elevation Deviation"
        self.description = "Calculates the maximum elevation deviation over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinDownslopeElevChange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Downslope Elev Change"
        self.description = "Calculates the minimum downslope change in elevation between a grid cell and its eight downslope neighbors."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MultiscaleRoughness(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiscale Roughness"
        self.description = "Calculates surface roughness over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MultiscaleRoughnessSignature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiscale Roughness Signature"
        self.description = "Calculates the surface roughness for points over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MultiscaleTopographicPositionImage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiscale Topographic Position Image"
        self.description = "Creates a multiscale topographic position image from three DEVmax rasters of differing spatial scale ranges."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NumDownslopeNeighbours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Num Downslope Neighbours"
        self.description = "Calculates the number of downslope neighbours to each grid cell in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NumUpslopeNeighbours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Num Upslope Neighbours"
        self.description = "Calculates the number of upslope neighbours to each grid cell in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PennockLandformClass(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Pennock Landform Class"
        self.description = "Classifies hillslope zones based on slope, profile curvature, and plan curvature."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PercentElevRange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Elev Range"
        self.description = "Calculates percent of elevation range from a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PlanCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Plan Curvature"
        self.description = "Calculates a plan (contour) curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Profile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Profile"
        self.description = "Plots profiles from digital surface models."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ProfileCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Profile Curvature"
        self.description = "Calculates a profile curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RelativeAspect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Relative Aspect"
        self.description = "Calculates relative aspect (relative to a user-specified direction) from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RelativeStreamPowerIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Relative Stream Power Index"
        self.description = "Calculates the relative stream power index."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RelativeTopographicPosition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Relative Topographic Position"
        self.description = "Calculates the relative topographic position index from a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RemoveOffTerrainObjects(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Off Terrain Objects"
        self.description = "Removes off-terrain objects from a raster digital elevation model (DEM)."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RuggednessIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ruggedness Index"
        self.description = "Calculates the Riley et al.'s (1999) terrain ruggedness index from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SedimentTransportIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sediment Transport Index"
        self.description = "Calculates the sediment transport index."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Slope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Slope"
        self.description = "Calculates a slope raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SlopeVsElevationPlot(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Slope Vs Elevation Plot"
        self.description = "Creates a slope vs. elevation plot for one or more DEMs."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StandardDeviationOfSlope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Standard Deviation Of Slope"
        self.description = "Calculates the standard deviation of slope from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SurfaceAreaRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Surface Area Ratio"
        self.description = "Calculates a the surface area ratio of each grid cell in an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TangentialCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tangential Curvature"
        self.description = "Calculates a tangential curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TotalCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Total Curvature"
        self.description = "Calculates a total curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Viewshed(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Viewshed"
        self.description = "Identifies the viewshed for a point or set of points."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class VisibilityIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Visibility Index"
        self.description = "Estimates the relative visibility of sites in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class WetnessIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Wetness Index"
        self.description = "Calculates the topographic wetness index, Ln(A / tan(slope))."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AverageFlowpathSlope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average Flowpath Slope"
        self.description = "Measures the average slope gradient from each grid cell to all upslope divide cells."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AverageUpslopeFlowpathLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average Upslope Flowpath Length"
        self.description = "Measures the average length of all upslope flowpaths draining each grid cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Basins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Basins"
        self.description = "Identifies drainage basins that drain to the DEM edge."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BreachDepressions(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Breach Depressions"
        self.description = "Breaches all of the depressions in a DEM using Lindsay's (2016) algorithm. This should be preferred over depression filling in most cases."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BreachSingleCellPits(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Breach Single Cell Pits"
        self.description = "Removes single-cell pits from an input DEM by breaching."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class D8FlowAccumulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D8 Flow Accumulation"
        self.description = "Calculates a D8 flow accumulation raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class D8MassFlux(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D8 Mass Flux"
        self.description = "Performs a D8 mass flux calculation."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class D8Pointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D8 Pointer"
        self.description = "Calculates a D8 flow pointer raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DInfFlowAccumulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D Inf Flow Accumulation"
        self.description = "Calculates a D-infinity flow accumulation raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DInfMassFlux(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D Inf Mass Flux"
        self.description = "Performs a D-infinity mass flux calculation."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DInfPointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D Inf Pointer"
        self.description = "Calculates a D-infinity flow pointer (flow direction) raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DepthInSink(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Depth In Sink"
        self.description = "Measures the depth of sinks (depressions) in a DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DownslopeDistanceToStream(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Downslope Distance To Stream"
        self.description = "Measures distance to the nearest downslope stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DownslopeFlowpathLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Downslope Flowpath Length"
        self.description = "Calculates the downslope flowpath length from each cell to basin outlet."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElevationAboveStream(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elevation Above Stream"
        self.description = "Calculates the elevation of cells above the nearest downslope stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ElevationAboveStreamEuclidean(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elevation Above Stream Euclidean"
        self.description = "Calculates the elevation of cells above the nearest (Euclidean distance) stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Fd8FlowAccumulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fd8 Flow Accumulation"
        self.description = "Calculates an FD8 flow accumulation raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Fd8Pointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fd8 Pointer"
        self.description = "Calculates an FD8 flow pointer raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FillBurn(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Burn"
        self.description = "Burns streams into a DEM using the FillBurn (Saunders, 1999) method."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FillDepressions(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Depressions"
        self.description = "Fills all of the depressions in a DEM. Depression breaching should be preferred in most cases."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FillSingleCellPits(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Single Cell Pits"
        self.description = "Raises pit cells to the elevation of their lowest neighbour."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindNoFlowCells(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find No Flow Cells"
        self.description = "Finds grid cells with no downslope neighbours."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindParallelFlow(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Parallel Flow"
        self.description = "Finds areas of parallel flow in D8 flow direction rasters."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FlattenLakes(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flatten Lakes"
        self.description = "Flattens lake polygons in a raster DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FloodOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flood Order"
        self.description = "Assigns each DEM grid cell its order in the sequence of inundations that are encountered during a search starting from the edges, moving inward at increasing elevations."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FlowAccumulationFullWorkflow(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flow Accumulation Full Workflow"
        self.description = "Resolves all of the depressions in a DEM, outputting a breached DEM, an aspect-aligned non-divergent flow pointer, and a flow accumulation raster."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FlowLengthDiff(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flow Length Diff"
        self.description = "Calculates the local maximum absolute difference in downslope flowpath length, useful in mapping drainage divides and ridges."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Hillslopes(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hillslopes"
        self.description = "Identifies the individual hillslopes draining to each link in a stream network."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ImpoundmentIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Impoundment Index"
        self.description = "Calculates the impoundment size resulting from damming a DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Isobasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Isobasins"
        self.description = "Divides a landscape into nearly equal sized drainage basins (i.e. watersheds)."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class JensonSnapPourPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Jenson Snap Pour Points"
        self.description = "Moves outlet points used to specify points of interest in a watershedding operation to the nearest stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LongestFlowpath(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Longest Flowpath"
        self.description = "Delineates the longest flowpaths for a group of subbasins or watersheds."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaxUpslopeFlowpathLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Upslope Flowpath Length"
        self.description = "Measures the maximum length of all upslope flowpaths draining each grid cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NumInflowingNeighbours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Num Inflowing Neighbours"
        self.description = "Computes the number of inflowing neighbours to each cell in an input DEM based on the D8 algorithm."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RaiseWalls(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raise Walls"
        self.description = "Raises walls in a DEM along a line or around a polygon, e.g. a watershed."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Rho8Pointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rho8 Pointer"
        self.description = "Calculates a stochastic Rho8 flow pointer raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Sink(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sink"
        self.description = "Identifies the depressions in a DEM, giving each feature a unique identifier."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SnapPourPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Snap Pour Points"
        self.description = "Moves outlet points used to specify points of interest in a watershedding operation to the cell with the highest flow accumulation in its neighbourhood."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StochasticDepressionAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stochastic Depression Analysis"
        self.description = "Preforms a stochastic analysis of depressions within a DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StrahlerOrderBasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Strahler Order Basins"
        self.description = "Identifies Strahler-order basins from an input stream network."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Subbasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Subbasins"
        self.description = "Identifies the catchments, or sub-basin, draining to each link in a stream network."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TraceDownslopeFlowpaths(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trace Downslope Flowpaths"
        self.description = "Traces downslope flowpaths from one or more target sites (i.e. seed points)."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class UnnestBasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Unnest Basins"
        self.description = "Extract whole watersheds for a set of outlet points."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Watershed(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Watershed"
        self.description = "Identifies the watershed, or drainage basin, draining to a set of target cells."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ChangeVectorAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Change Vector Analysis"
        self.description = "Performs a change vector analysis on a two-date multi-spectral dataset."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Closing(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Closing"
        self.description = "A closing is a mathematical morphology operation involving an erosion (min filter) of a dilation (max filter) set."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CreateColourComposite(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Colour Composite"
        self.description = "Creates a colour-composite image from three bands of multispectral imagery."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FlipImage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flip Image"
        self.description = "Reflects an image in the vertical or horizontal axis."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class IhsToRgb(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ihs To Rgb"
        self.description = "Converts intensity, hue, and saturation (IHS) images into red, green, and blue (RGB) images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ImageStackProfile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Stack Profile"
        self.description = "Plots an image stack profile (i.e. signature) for a set of points and multispectral images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class IntegralImage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Integral Image"
        self.description = "Transforms an input image (summed area table) into its integral image equivalent."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class KMeansClustering(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "K Means Clustering"
        self.description = "Performs a k-means clustering operation on a multi-spectral dataset."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LineThinning(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Line Thinning"
        self.description = "Performs line thinning a on Boolean raster image; intended to be used with the RemoveSpurs tool."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ModifiedKMeansClustering(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Modified K Means Clustering"
        self.description = "Performs a modified k-means clustering operation on a multi-spectral dataset."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Mosaic(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mosaic"
        self.description = "Mosaics two or more images together."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MosaicWithFeathering(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mosaic With Feathering"
        self.description = "Mosaics two images together using a feathering technique in overlapping areas to reduce edge-effects."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NormalizedDifferenceVegetationIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Normalized Difference Vegetation Index"
        self.description = "Calculates the normalized difference vegetation index (NDVI) from near-infrared and red imagery."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Opening(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Opening"
        self.description = "An opening is a mathematical morphology operation involving a dilation (max filter) of an erosion (min filter) set."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RemoveSpurs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Spurs"
        self.description = "Removes the spurs (pruning operation) from a Boolean line image.; intended to be used on the output of the LineThinning tool."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Resample(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Resample"
        self.description = "Resamples one or more input images into a destination image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RgbToIhs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rgb To Ihs"
        self.description = "Converts red, green, and blue (RGB) images into intensity, hue, and saturation (IHS) images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SplitColourComposite(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Split Colour Composite"
        self.description = "This tool splits an RGB colour composite image into seperate multispectral images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ThickenRasterLine(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Thicken Raster Line"
        self.description = "Thickens single-cell wide lines within a raster image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TophatTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tophat Transform"
        self.description = "Performs either a white or black top-hat transform on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class WriteFunctionMemoryInsertion(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Write Function Memory Insertion"
        self.description = "Performs a write function memory insertion for single-band multi-date change detection."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AdaptiveFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Adaptive Filter"
        self.description = "Performs an adaptive filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BilateralFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bilateral Filter"
        self.description = "A bilateral filter is an edge-preserving smoothing filter introduced by Tomasi and Manduchi (1998)."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ConservativeSmoothingFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Conservative Smoothing Filter"
        self.description = "Performs a conservative-smoothing filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CornerDetection(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Corner Detection"
        self.description = "Identifies corner patterns in boolean images using hit-and-miss pattern mattching."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DiffOfGaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Diff Of Gaussian Filter"
        self.description = "Performs a Difference of Gaussian (DoG) filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DiversityFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Diversity Filter"
        self.description = "Assigns each cell in the output grid the number of different values in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EdgePreservingMeanFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Preserving Mean Filter"
        self.description = "Performs a simple edge-preserving mean filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EmbossFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Emboss Filter"
        self.description = "Performs an emboss filter on an image, similar to a hillshade operation."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FastAlmostGaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fast Almost Gaussian Filter"
        self.description = "Performs a fast approximate Gaussian filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class GaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gaussian Filter"
        self.description = "Performs a Gaussian filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HighPassFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "High Pass Filter"
        self.description = "Performs a high-pass filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HighPassMedianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "High Pass Median Filter"
        self.description = "Performs a high pass median filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class KNearestMeanFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "K Nearest Mean Filter"
        self.description = "A k-nearest mean filter is a type of edge-preserving smoothing filter."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LaplacianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Laplacian Filter"
        self.description = "Performs a Laplacian filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LaplacianOfGaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Laplacian Of Gaussian Filter"
        self.description = "Performs a Laplacian-of-Gaussian (LoG) filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LeeFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lee Filter"
        self.description = "Performs a Lee (Sigma) smoothing filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LineDetectionFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Line Detection Filter"
        self.description = "Performs a line-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MajorityFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Majority Filter"
        self.description = "Assigns each cell in the output grid the most frequently occurring value (mode) in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MaximumFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Maximum Filter"
        self.description = "Assigns each cell in the output grid the maximum value in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MeanFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mean Filter"
        self.description = "Performs a mean filter (low-pass filter) on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MedianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Median Filter"
        self.description = "Performs a median filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinimumFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Filter"
        self.description = "Assigns each cell in the output grid the minimum value in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class OlympicFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Olympic Filter"
        self.description = "Performs an olympic smoothing filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PercentileFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percentile Filter"
        self.description = "Performs a percentile filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PrewittFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Prewitt Filter"
        self.description = "Performs a Prewitt edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RangeFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Range Filter"
        self.description = "Assigns each cell in the output grid the range of values in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RobertsCrossFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Roberts Cross Filter"
        self.description = "Performs a Robert's cross edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ScharrFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Scharr Filter"
        self.description = "Performs a Scharr edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SobelFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sobel Filter"
        self.description = "Performs a Sobel edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StandardDeviationFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Standard Deviation Filter"
        self.description = "Assigns each cell in the output grid the standard deviation of values in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TotalFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Total Filter"
        self.description = "Performs a total filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class UnsharpMasking(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Unsharp Masking"
        self.description = "An image sharpening technique that enhances edges."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class UserInedWeightsFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "User Ined Weights Filter"
        self.description = "Performs a user-defined weights filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class BalanceContrastEnhancement(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Balance Contrast Enhancement"
        self.description = "Performs a balance contrast enhancement on a colour-composite image of multispectral data."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CorrectVignetting(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Correct Vignetting"
        self.description = "Corrects the darkening of images towards corners."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DirectDecorrelationStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Direct Decorrelation Stretch"
        self.description = "Performs a direct decorrelation stretch enhancement on a colour-composite image of multispectral data."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class GammaCorrection(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gamma Correction"
        self.description = "Performs a sigmoidal contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class GaussianContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gaussian Contrast Stretch"
        self.description = "Performs a Gaussian contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HistogramEqualization(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram Equalization"
        self.description = "Performs a histogram equalization contrast enhancment on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HistogramMatching(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram Matching"
        self.description = "Alters the statistical distribution of a raster image matching it to a specified PDF."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HistogramMatchingTwoImages(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram Matching Two Images"
        self.description = "This tool alters the cumulative distribution function of a raster image to that of another image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class MinMaxContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Max Contrast Stretch"
        self.description = "Performs a min-max contrast stretch on an input greytone image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PanchromaticSharpening(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Panchromatic Sharpening"
        self.description = "Increases the spatial resolution of image data by combining multispectral bands with panchromatic data."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PercentageContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percentage Contrast Stretch"
        self.description = "Performs a percentage linear contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SigmoidalContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sigmoidal Contrast Stretch"
        self.description = "Performs a sigmoidal contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StandardDeviationContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Standard Deviation Contrast Stretch"
        self.description = "Performs a standard-deviation contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ClassifyOverlapPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Classify Overlap Points"
        self.description = "Classifies or filters LAS points in regions of overlapping flight lines."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ClipLidarToPolygon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clip Lidar To Polygon"
        self.description = "Clips a LiDAR point cloud to a vector polygon or polygons."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ErasePolygonFromLidar(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Erase Polygon From Lidar"
        self.description = "Erases (cuts out) a vector polygon or polygons from a LiDAR point cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FilterLidarScanAngles(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Filter Lidar Scan Angles"
        self.description = "Removes points in a LAS file with scan angles greater than a threshold."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindFlightlineEdgePoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Flightline Edge Points"
        self.description = "Identifies points along a flightline's edge in a LAS file."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FlightlineOverlap(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flightline Overlap"
        self.description = "Reads a LiDAR (LAS) point file and outputs a raster containing the number of overlapping flight lines in each grid cell."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LasToAscii(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Las To Ascii"
        self.description = "Converts one or more LAS files into ASCII text files."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LasToMultipointShapefile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Las To Multipoint Shapefile"
        self.description = "Converts one or more LAS files into MultipointZ vector Shapefiles. When the input parameter is not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LasToShapefile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Las To Shapefile"
        self.description = "Converts one or more LAS files into a vector Shapefile of POINT ShapeType."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarBlockMaximum(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Block Maximum"
        self.description = "Creates a block-maximum raster from an input LAS file. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarBlockMinimum(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Block Minimum"
        self.description = "Creates a block-minimum raster from an input LAS file. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarClassifySubset(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Classify Subset"
        self.description = "Classifies the values in one LiDAR point cloud that correpond with points in a subset cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarColourize(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Colourize"
        self.description = "Adds the red-green-blue colour fields of a LiDAR (LAS) file based on an input image."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarConstructVectorTin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Construct Vector Tin"
        self.description = "Creates a vector triangular irregular network (TIN) fitted to LiDAR points."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarElevationSlice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Elevation Slice"
        self.description = "Outputs all of the points within a LiDAR (LAS) point file that lie between a specified elevation range."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarGroundPointFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Ground Point Filter"
        self.description = "Identifies ground points within LiDAR dataset using a slope-based method."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarHexBinning(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Hex Binning"
        self.description = "Hex-bins a set of LiDAR points."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarHillshade(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Hillshade"
        self.description = "Calculates a hillshade value for points within a LAS file and stores these data in the RGB field."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarHistogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Histogram"
        self.description = "Creates a histogram of LiDAR data."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarIdwInterpolation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Idw Interpolation"
        self.description = "Interpolates LAS files using an inverse-distance weighted (IDW) scheme. When the input/output parameters are not specified, the tool interpolates all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarInfo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Info"
        self.description = "Prints information about a LiDAR (LAS) dataset, including header, point return frequency, and classification data and information about the variable length records (VLRs) and geokeys."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarJoin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Join"
        self.description = "Joins multiple LiDAR (LAS) files into a single LAS file."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarKappaIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Kappa Index"
        self.description = "Performs a kappa index of agreement (KIA) analysis on the classifications of two LAS files."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarNearestNeighbourGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Nearest Neighbour Gridding"
        self.description = "Grids LAS files using nearest-neighbour scheme. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarPointDensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Point Density"
        self.description = "Calculates the spatial pattern of point density for a LiDAR data set. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarPointStats(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Point Stats"
        self.description = "Creates several rasters summarizing the distribution of LAS point data. When the input/output parameters are not specified, the tool works on all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarRemoveDuplicates(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Remove Duplicates"
        self.description = "Removes duplicate points from a LiDAR data set."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarRemoveOutliers(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Remove Outliers"
        self.description = "Removes outliers (high and low points) in a LiDAR point cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarSegmentation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Segmentation"
        self.description = "Segments a LiDAR point cloud based on normal vectors."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarSegmentationBasedFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Segmentation Based Filter"
        self.description = "Identifies ground points within LiDAR point clouds using a segmentation based approach."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarThin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Thin"
        self.description = "Thins a LiDAR point cloud, reducing point density."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarThinHighDensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Thin High Density"
        self.description = "Thins points from high density areas within a LiDAR point cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarTile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tile"
        self.description = "Tiles a LiDAR LAS file into multiple LAS files."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarTileFootprint(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tile Footprint"
        self.description = "Creates a vector polygon of the convex hull of a LiDAR point cloud. When the input/output parameters are not specified, the tool works with all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarTinGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tin Gridding"
        self.description = "Creates a raster grid based on a Delaunay triangular irregular network (TIN) fitted to LiDAR points."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LidarTophatTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tophat Transform"
        self.description = "Performs a white top-hat transform on a Lidar dataset; as an estimate of height above ground, this is useful for modelling the vegetation canopy."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NormalVectors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Normal Vectors"
        self.description = "Calculates normal vectors for points within a LAS file and stores these data (XYZ vector components) in the RGB field."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SelectTilesByPolygon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Select Tiles By Polygon"
        self.description = "Copies LiDAR tiles overlapping with a polygon into an output directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class And(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "And"
        self.description = "Performs a logical AND operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Not(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Not"
        self.description = "Performs a logical NOT operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Or(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Or"
        self.description = "Performs a logical OR operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AbsoluteValue(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Absolute Value"
        self.description = "Calculates the absolute value of every cell in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Add(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add"
        self.description = "Performs an addition operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Anova(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Anova"
        self.description = "Performs an analysis of variance (ANOVA) test on a raster dataset."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ArcCos(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Arc Cos"
        self.description = "Returns the inverse cosine (arccos) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ArcSin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Arc Sin"
        self.description = "Returns the inverse sine (arcsin) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ArcTan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Arc Tan"
        self.description = "Returns the inverse tangent (arctan) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Atan2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Atan2"
        self.description = "Returns the 2-argument inverse tangent (atan2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AttributeCorrelation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Attribute Correlation"
        self.description = "Performs a correlation analysis on attribute fields from a vector database."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AttributeHistogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Attribute Histogram"
        self.description = "Creates a histogram for the field values of a vector's attribute table."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class AttributeScattergram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Attribute Scattergram"
        self.description = "Creates a scattergram for two field values of a vector's attribute table."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Ceil(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ceil"
        self.description = "Returns the smallest (closest to negative infinity) value that is greater than or equal to the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Cos(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cos"
        self.description = "Returns the cosine (cos) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Cosh(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cosh"
        self.description = "Returns the hyperbolic cosine (cosh) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CrispnessIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Crispness Index"
        self.description = "Calculates the Crispness Index, which is used to quantify how crisp (or conversely how fuzzy) a probability image is."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CrossTabulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cross Tabulation"
        self.description = "Performs a cross-tabulation on two categorical images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class CumulativeDistribution(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cumulative Distribution"
        self.description = "Converts a raster image to its cumulative distribution function."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Decrement(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Decrement"
        self.description = "Decreases the values of each grid cell in an input raster by 1.0 (see also InPlaceSubtract)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Divide(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Divide"
        self.description = "Performs a division operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class EqualTo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Equal To"
        self.description = "Performs a equal-to comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Exp(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Exp"
        self.description = "Returns the exponential (base e) of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Exp2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Exp2"
        self.description = "Returns the exponential (base 2) of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExtractRasterStatistics(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Raster Statistics"
        self.description = "Extracts descriptive statistics for a group of patches in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Floor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Floor"
        self.description = "Returns the largest (closest to positive infinity) value that is less than or equal to the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class GreaterThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Greater Than"
        self.description = "Performs a greater-than comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ImageAutocorrelation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Autocorrelation"
        self.description = "Performs Moran's I analysis on two or more input images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ImageCorrelation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Correlation"
        self.description = "Performs image correlation on two or more input images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ImageRegression(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Regression"
        self.description = "Performs image regression analysis on two input images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class InPlaceAdd(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Add"
        self.description = "Performs an in-place addition operation (input1 += input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class InPlaceDivide(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Divide"
        self.description = "Performs an in-place division operation (input1 /= input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class InPlaceMultiply(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Multiply"
        self.description = "Performs an in-place multiplication operation (input1 *= input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class InPlaceSubtract(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Subtract"
        self.description = "Performs an in-place subtraction operation (input1 -= input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Increment(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Increment"
        self.description = "Increases the values of each grid cell in an input raster by 1.0. (see also InPlaceAdd)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class IntegerDivision(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Integer Division"
        self.description = "Performs an integer division operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class IsNoData(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Is No Data"
        self.description = "Identifies NoData valued pixels in an image."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class KappaIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Kappa Index"
        self.description = "Performs a kappa index of agreement (KIA) analysis on two categorical raster files."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class KsTestForNormality(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ks Test For Normality"
        self.description = "Evaluates whether the values in a raster are normally distributed."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LessThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Less Than"
        self.description = "Performs a less-than comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ListUniqueValues(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "List Unique Values"
        self.description = "Lists the unique values contained in a field witin a vector's attribute table."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Ln(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ln"
        self.description = "Returns the natural logarithm of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Log10(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Log10"
        self.description = "Returns the base-10 logarithm of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Log2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Log2"
        self.description = "Returns the base-2 logarithm of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Max(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max"
        self.description = "Performs a MAX operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Min(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min"
        self.description = "Performs a MIN operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Modulo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Modulo"
        self.description = "Performs a modulo operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Multiply(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiply"
        self.description = "Performs a multiplication operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Negate(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Negate"
        self.description = "Changes the sign of values in a raster or the 0-1 values of a Boolean raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class NotEqualTo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Not Equal To"
        self.description = "Performs a not-equal-to comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Power(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Power"
        self.description = "Raises the values in grid cells of one rasters, or a constant value, by values in another raster or constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class PrincipalComponentAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Principal Component Analysis"
        self.description = "Performs a principal component analysis (PCA) on a multi-spectral dataset."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Quantiles(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Quantiles"
        self.description = "Transforms raster values into quantiles."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RandomField(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Random Field"
        self.description = "Creates an image containing random values."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RandomSample(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Random Sample"
        self.description = "Creates an image containing randomly located sample grid cells with unique IDs."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterHistogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Histogram"
        self.description = "Creates a histogram from raster values."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterSummaryStats(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Summary Stats"
        self.description = "Measures a rasters min, max, average, standard deviation, num. non-nodata cells, and total."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Reciprocal(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reciprocal"
        self.description = "Returns the reciprocal (i.e. 1 / z) of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RescaleValueRange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rescale Value Range"
        self.description = "Performs a min-max contrast stretch on an input greytone image."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RootMeanSquareError(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Root Mean Square Error"
        self.description = "Calculates the RMSE and other accuracy statistics."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Round(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Round"
        self.description = "Rounds the values in an input raster to the nearest integer value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Sin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sin"
        self.description = "Returns the sine (sin) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Sinh(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sinh"
        self.description = "Returns the hyperbolic sine (sinh) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Square(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Square"
        self.description = "Squares the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class SquareRoot(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Square Root"
        self.description = "Returns the square root of the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Subtract(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Subtract"
        self.description = "Performs a differencing operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Tan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tan"
        self.description = "Returns the tangent (tan) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Tanh(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tanh"
        self.description = "Returns the hyperbolic tangent (tanh) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ToDegrees(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "To Degrees"
        self.description = "Converts a raster from radians to degrees."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ToRadians(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "To Radians"
        self.description = "Converts a raster from degrees to radians."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TrendSurface(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trend Surface"
        self.description = "Estimates the trend surface of an input raster file."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TrendSurfaceVectorPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trend Surface Vector Points"
        self.description = "Estimates a trend surface from vector points."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Truncate(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Truncate"
        self.description = "Truncates the values in a raster to the desired number of decimal places."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TurningBandsSimulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Turning Bands Simulation"
        self.description = "Creates an image containing random values based on a turning-bands simulation."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class Xor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Xor"
        self.description = "Performs a logical XOR operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ZScores(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Z Scores"
        self.description = "Standardizes the values in an input raster by converting to z-scores."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class DistanceToOutlet(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Distance To Outlet"
        self.description = "Calculates the distance of stream grid cells to the channel network outlet cell."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExtractStreams(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Streams"
        self.description = "Extracts stream grid cells from a flow accumulation raster."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ExtractValleys(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Valleys"
        self.description = "Identifies potential valley bottom grid cells based on local topolography alone."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FarthestChannelHead(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Farthest Channel Head"
        self.description = "Calculates the distance to the furthest upstream channel head for each stream cell."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class FindMainStem(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Main Stem"
        self.description = "Finds the main stem, based on stream lengths, of each stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HackStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hack Stream Order"
        self.description = "Assigns the Hack stream order to each tributary in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class HortonStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Horton Stream Order"
        self.description = "Assigns the Horton stream order to each tributary in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LengthOfUpstreamChannels(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Length Of Upstream Channels"
        self.description = "Calculates the total length of channels upstream."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LongProfile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Long Profile"
        self.description = "Plots the stream longitudinal profiles for one or more rivers."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class LongProfileFromPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Long Profile From Points"
        self.description = "Plots the longitudinal profiles from flow-paths initiating from a set of vector points."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterStreamsToVector(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Streams To Vector"
        self.description = "Converts a raster stream file into a vector file."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RasterizeStreams(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rasterize Streams"
        self.description = "Rasterizes vector streams based on Lindsay (2016) method."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class RemoveShortStreams(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Short Streams"
        self.description = "Removes short first-order streams from a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class ShreveStreamMagnitude(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shreve Stream Magnitude"
        self.description = "Assigns the Shreve stream magnitude to each link in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StrahlerStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Strahler Stream Order"
        self.description = "Assigns the Strahler stream order to each link in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StreamLinkClass(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Class"
        self.description = "Identifies the exterior/interior links and nodes in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StreamLinkIdentifier(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Identifier"
        self.description = "Assigns a unique identifier to each link in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StreamLinkLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Length"
        self.description = "Estimates the length of each link (or tributary) in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StreamLinkSlope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Slope"
        self.description = "Estimates the average slope of each link (or tributary) in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class StreamSlopeContinuous(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Slope Continuous"
        self.description = "Estimates the slope of each grid cell in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TopologicalStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Topological Stream Order"
        self.description = "Assigns each link in a stream network its topological order."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return


class TributaryIdentifier(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tributary Identifier"
        self.description = "Assigns a unique identifier to each tributary in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return



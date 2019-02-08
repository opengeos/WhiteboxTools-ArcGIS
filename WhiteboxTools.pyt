import arcpy


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        # self.tools = [Tool]
        tools = []

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

class AddPointCoordinatesToTable(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add Point Coordinates To Table"
        self.description = ""
        self.category = "Data Tools"

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

class ConvertNodataToZero(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert Nodata To Zero"
        self.description = ""
        self.category = "Data Tools"

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

class ConvertRasterFormat(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert Raster Format"
        self.description = ""
        self.category = "Data Tools"

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

class ExportTableToCsv(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Export Table To Csv"
        self.description = ""
        self.category = "Data Tools"

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

class JoinTables(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Join Tables"
        self.description = ""
        self.category = "Data Tools"

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

class LinesToPolygons(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lines To Polygons"
        self.description = ""
        self.category = "Data Tools"

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

class MergeTableWithCsv(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Merge Table With Csv"
        self.description = ""
        self.category = "Data Tools"

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

class MergeVectors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Merge Vectors"
        self.description = ""
        self.category = "Data Tools"

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

class MultiPartToSinglePart(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multi Part To Single Part"
        self.description = ""
        self.category = "Data Tools"

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

class NewRasterFromBase(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "New Raster From Base"
        self.description = ""
        self.category = "Data Tools"

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

class PolygonsToLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygons To Lines"
        self.description = ""
        self.category = "Data Tools"

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

class PrintGeoTiffTags(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Print Geo Tiff Tags"
        self.description = ""
        self.category = "Data Tools"

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

class RasterToVectorLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster To Vector Lines"
        self.description = ""
        self.category = "Data Tools"

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

class RasterToVectorPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster To Vector Points"
        self.description = ""
        self.category = "Data Tools"

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

class ReinitializeAttributeTable(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reinitialize Attribute Table"
        self.description = ""
        self.category = "Data Tools"

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

class RemovePolygonHoles(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Polygon Holes"
        self.description = ""
        self.category = "Data Tools"

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

class SetNodataValue(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Set Nodata Value"
        self.description = ""
        self.category = "Data Tools"

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

class SinglePartToMultiPart(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Single Part To Multi Part"
        self.description = ""
        self.category = "Data Tools"

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

class VectorLinesToRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Lines To Raster"
        self.description = ""
        self.category = "Data Tools"

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

class VectorPointsToRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Points To Raster"
        self.description = ""
        self.category = "Data Tools"

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

class VectorPolygonsToRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Polygons To Raster"
        self.description = ""
        self.category = "Data Tools"

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

class AggregateRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Aggregate Raster"
        self.description = ""
        self.category = "GIS Analysis"

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

class BlockMaximumGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Block Maximum Gridding"
        self.description = ""
        self.category = "GIS Analysis"

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

class BlockMinimumGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Block Minimum Gridding"
        self.description = ""
        self.category = "GIS Analysis"

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

class Centroid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Centroid"
        self.description = ""
        self.category = "GIS Analysis"

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

class CentroidVector(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Centroid Vector"
        self.description = ""
        self.category = "GIS Analysis"

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

class Clump(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clump"
        self.description = ""
        self.category = "GIS Analysis"

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

class ConstructVectorTin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Construct Vector Tin"
        self.description = ""
        self.category = "GIS Analysis"

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

class CreateHexagonalVectorGrid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Hexagonal Vector Grid"
        self.description = ""
        self.category = "GIS Analysis"

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

class CreatePlane(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Plane"
        self.description = ""
        self.category = "GIS Analysis"

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

class CreateRectangularVectorGrid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Rectangular Vector Grid"
        self.description = ""
        self.category = "GIS Analysis"

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

class Dissolve(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Dissolve"
        self.description = ""
        self.category = "GIS Analysis"

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

class EliminateCoincidentPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Eliminate Coincident Points"
        self.description = ""
        self.category = "GIS Analysis"

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

class ExtendVectorLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extend Vector Lines"
        self.description = ""
        self.category = "GIS Analysis"

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

class ExtractNodes(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Nodes"
        self.description = ""
        self.category = "GIS Analysis"

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

class ExtractRasterValuesAtPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Raster Values At Points"
        self.description = ""
        self.category = "GIS Analysis"

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

class FindLowestOrHighestPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Lowest Or Highest Points"
        self.description = ""
        self.category = "GIS Analysis"

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

class IdwInterpolation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Idw Interpolation"
        self.description = ""
        self.category = "GIS Analysis"

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

class LayerFootprint(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Layer Footprint"
        self.description = ""
        self.category = "GIS Analysis"

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

class Medoid(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Medoid"
        self.description = ""
        self.category = "GIS Analysis"

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

class MinimumBoundingBox(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Bounding Box"
        self.description = ""
        self.category = "GIS Analysis"

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

class MinimumBoundingCircle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Bounding Circle"
        self.description = ""
        self.category = "GIS Analysis"

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

class MinimumBoundingEnvelope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Bounding Envelope"
        self.description = ""
        self.category = "GIS Analysis"

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

class MinimumConvexHull(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Convex Hull"
        self.description = ""
        self.category = "GIS Analysis"

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

class NearestNeighbourGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Nearest Neighbour Gridding"
        self.description = ""
        self.category = "GIS Analysis"

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

class PolygonArea(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Area"
        self.description = ""
        self.category = "GIS Analysis"

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

class PolygonLongAxis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Long Axis"
        self.description = ""
        self.category = "GIS Analysis"

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

class PolygonPerimeter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Perimeter"
        self.description = ""
        self.category = "GIS Analysis"

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

class PolygonShortAxis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygon Short Axis"
        self.description = ""
        self.category = "GIS Analysis"

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

class RasterCellAssignment(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Cell Assignment"
        self.description = ""
        self.category = "GIS Analysis"

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

class Reclass(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reclass"
        self.description = ""
        self.category = "GIS Analysis"

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

class ReclassEqualInterval(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reclass Equal Interval"
        self.description = ""
        self.category = "GIS Analysis"

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

class ReclassFromFile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reclass From File"
        self.description = ""
        self.category = "GIS Analysis"

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

class SmoothVectors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Smooth Vectors"
        self.description = ""
        self.category = "GIS Analysis"

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

class TinGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tin Gridding"
        self.description = ""
        self.category = "GIS Analysis"

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

class VectorHexBinning(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Vector Hex Binning"
        self.description = ""
        self.category = "GIS Analysis"

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

class VoronoiDiagram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Voronoi Diagram"
        self.description = ""
        self.category = "GIS Analysis"

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

class BufferRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Buffer Raster"
        self.description = ""
        self.category = "GIS Analysis"

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

class CostAllocation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cost Allocation"
        self.description = ""
        self.category = "GIS Analysis"

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

class CostDistance(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cost Distance"
        self.description = ""
        self.category = "GIS Analysis"

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

class CostPathway(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cost Pathway"
        self.description = ""
        self.category = "GIS Analysis"

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

class EuclideanAllocation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Euclidean Allocation"
        self.description = ""
        self.category = "GIS Analysis"

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

class EuclideanDistance(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Euclidean Distance"
        self.description = ""
        self.category = "GIS Analysis"

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

class AverageOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class Clip(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clip"
        self.description = ""
        self.category = "GIS Analysis"

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

class ClipRasterToPolygon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clip Raster To Polygon"
        self.description = ""
        self.category = "GIS Analysis"

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

class CountIf(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Count If"
        self.description = ""
        self.category = "GIS Analysis"

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

class Difference(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Difference"
        self.description = ""
        self.category = "GIS Analysis"

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

class Erase(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Erase"
        self.description = ""
        self.category = "GIS Analysis"

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

class ErasePolygonFromRaster(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Erase Polygon From Raster"
        self.description = ""
        self.category = "GIS Analysis"

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

class HighestPosition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Highest Position"
        self.description = ""
        self.category = "GIS Analysis"

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

class Intersect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Intersect"
        self.description = ""
        self.category = "GIS Analysis"

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

class LineIntersections(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Line Intersections"
        self.description = ""
        self.category = "GIS Analysis"

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

class LowestPosition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lowest Position"
        self.description = ""
        self.category = "GIS Analysis"

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

class MaxAbsoluteOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Absolute Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class MaxOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class MinAbsoluteOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Absolute Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class MinOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class PercentEqualTo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Equal To"
        self.description = ""
        self.category = "GIS Analysis"

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

class PercentGreaterThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Greater Than"
        self.description = ""
        self.category = "GIS Analysis"

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

class PercentLessThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Less Than"
        self.description = ""
        self.category = "GIS Analysis"

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

class PickFromList(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Pick From List"
        self.description = ""
        self.category = "GIS Analysis"

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

class Polygonize(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Polygonize"
        self.description = ""
        self.category = "GIS Analysis"

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

class SplitWithLines(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Split With Lines"
        self.description = ""
        self.category = "GIS Analysis"

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

class SumOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sum Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class SymmetricalDifference(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Symmetrical Difference"
        self.description = ""
        self.category = "GIS Analysis"

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

class Union(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Union"
        self.description = ""
        self.category = "GIS Analysis"

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

class WeightedOverlay(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Weighted Overlay"
        self.description = ""
        self.category = "GIS Analysis"

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

class WeightedSum(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Weighted Sum"
        self.description = ""
        self.category = "GIS Analysis"

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

class CompactnessRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Compactness Ratio"
        self.description = ""
        self.category = "GIS Analysis"

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

class EdgeProportion(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Proportion"
        self.description = ""
        self.category = "GIS Analysis"

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

class ElongationRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elongation Ratio"
        self.description = ""
        self.category = "GIS Analysis"

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

class FindPatchOrClassEdgeCells(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Patch Or Class Edge Cells"
        self.description = ""
        self.category = "GIS Analysis"

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

class HoleProportion(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hole Proportion"
        self.description = ""
        self.category = "GIS Analysis"

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

class LinearityIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Linearity Index"
        self.description = ""
        self.category = "GIS Analysis"

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

class PatchOrientation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Patch Orientation"
        self.description = ""
        self.category = "GIS Analysis"

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

class PerimeterAreaRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Perimeter Area Ratio"
        self.description = ""
        self.category = "GIS Analysis"

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

class RadiusOfGyration(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Radius Of Gyration"
        self.description = ""
        self.category = "GIS Analysis"

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

class RelatedCircumscribingCircle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Related Circumscribing Circle"
        self.description = ""
        self.category = "GIS Analysis"

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

class ShapeComplexityIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shape Complexity Index"
        self.description = ""
        self.category = "GIS Analysis"

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

class Aspect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Aspect"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class CircularVarianceOfAspect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Circular Variance Of Aspect"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class DevFromMeanElev(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Dev From Mean Elev"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class DiffFromMeanElev(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Diff From Mean Elev"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class DirectionalRelief(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Directional Relief"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class DownslopeIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Downslope Index"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class DrainagePreservingSmoothing(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Drainage Preserving Smoothing"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class EdgeDensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Density"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class ElevAbovePit(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Above Pit"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class ElevPercentile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Percentile"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class ElevRelativeToMinMax(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Relative To Min Max"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class ElevRelativeToWatershedMinMax(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elev Relative To Watershed Min Max"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class FeaturePreservingDenoise(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Feature Preserving Denoise"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class FetchAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fetch Analysis"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class FillMissingData(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Missing Data"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class FindRidges(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Ridges"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class Hillshade(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hillshade"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class HorizonAngle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Horizon Angle"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class HypsometricAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hypsometric Analysis"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxAnisotropyDev(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Anisotropy Dev"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxAnisotropyDevSignature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Anisotropy Dev Signature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxBranchLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Branch Length"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxDifferenceFromMean(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Difference From Mean"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxDownslopeElevChange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Downslope Elev Change"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxElevDevSignature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Elev Dev Signature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MaxElevationDeviation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Elevation Deviation"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MinDownslopeElevChange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Downslope Elev Change"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MultiscaleRoughness(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiscale Roughness"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MultiscaleRoughnessSignature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiscale Roughness Signature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class MultiscaleTopographicPositionImage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiscale Topographic Position Image"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class NumDownslopeNeighbours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Num Downslope Neighbours"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class NumUpslopeNeighbours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Num Upslope Neighbours"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class PennockLandformClass(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Pennock Landform Class"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class PercentElevRange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percent Elev Range"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class PlanCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Plan Curvature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class Profile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Profile"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class ProfileCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Profile Curvature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class RelativeAspect(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Relative Aspect"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class RelativeStreamPowerIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Relative Stream Power Index"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class RelativeTopographicPosition(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Relative Topographic Position"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class RemoveOffTerrainObjects(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Off Terrain Objects"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class RuggednessIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ruggedness Index"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class SedimentTransportIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sediment Transport Index"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class Slope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Slope"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class SlopeVsElevationPlot(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Slope Vs Elevation Plot"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class StandardDeviationOfSlope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Standard Deviation Of Slope"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class SurfaceAreaRatio(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Surface Area Ratio"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class TangentialCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tangential Curvature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class TotalCurvature(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Total Curvature"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class Viewshed(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Viewshed"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class VisibilityIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Visibility Index"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class WetnessIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Wetness Index"
        self.description = ""
        self.category = "Geomorphometric Analysis"

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

class AverageFlowpathSlope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average Flowpath Slope"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class AverageUpslopeFlowpathLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Average Upslope Flowpath Length"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Basins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Basins"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class BreachDepressions(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Breach Depressions"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class BreachSingleCellPits(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Breach Single Cell Pits"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class D8FlowAccumulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D8 Flow Accumulation"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class D8MassFlux(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D8 Mass Flux"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class D8Pointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D8 Pointer"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class DInfFlowAccumulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D Inf Flow Accumulation"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class DInfMassFlux(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D Inf Mass Flux"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class DInfPointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "D Inf Pointer"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class DepthInSink(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Depth In Sink"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class DownslopeDistanceToStream(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Downslope Distance To Stream"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class DownslopeFlowpathLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Downslope Flowpath Length"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class ElevationAboveStream(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elevation Above Stream"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class ElevationAboveStreamEuclidean(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Elevation Above Stream Euclidean"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Fd8FlowAccumulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fd8 Flow Accumulation"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Fd8Pointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fd8 Pointer"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FillBurn(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Burn"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FillDepressions(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Depressions"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FillSingleCellPits(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fill Single Cell Pits"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FindNoFlowCells(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find No Flow Cells"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FindParallelFlow(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Parallel Flow"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FlattenLakes(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flatten Lakes"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FloodOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flood Order"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FlowAccumulationFullWorkflow(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flow Accumulation Full Workflow"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class FlowLengthDiff(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flow Length Diff"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Hillslopes(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hillslopes"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class ImpoundmentIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Impoundment Index"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Isobasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Isobasins"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class JensonSnapPourPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Jenson Snap Pour Points"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class LongestFlowpath(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Longest Flowpath"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class MaxUpslopeFlowpathLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max Upslope Flowpath Length"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class NumInflowingNeighbours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Num Inflowing Neighbours"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class RaiseWalls(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raise Walls"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Rho8Pointer(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rho8 Pointer"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Sink(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sink"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class SnapPourPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Snap Pour Points"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class StochasticDepressionAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stochastic Depression Analysis"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class StrahlerOrderBasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Strahler Order Basins"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Subbasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Subbasins"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class TraceDownslopeFlowpaths(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trace Downslope Flowpaths"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class UnnestBasins(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Unnest Basins"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class Watershed(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Watershed"
        self.description = ""
        self.category = "Hydrological Analysis"

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

class ChangeVectorAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Change Vector Analysis"
        self.description = ""
        self.category = "Image Processing Tools"

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

class Closing(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Closing"
        self.description = ""
        self.category = "Image Processing Tools"

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

class CreateColourComposite(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Colour Composite"
        self.description = ""
        self.category = "Image Processing Tools"

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

class FlipImage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flip Image"
        self.description = ""
        self.category = "Image Processing Tools"

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

class IhsToRgb(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ihs To Rgb"
        self.description = ""
        self.category = "Image Processing Tools"

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

class ImageStackProfile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Stack Profile"
        self.description = ""
        self.category = "Image Processing Tools"

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

class IntegralImage(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Integral Image"
        self.description = ""
        self.category = "Image Processing Tools"

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

class KMeansClustering(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "K Means Clustering"
        self.description = ""
        self.category = "Image Processing Tools"

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

class LineThinning(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Line Thinning"
        self.description = ""
        self.category = "Image Processing Tools"

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

class ModifiedKMeansClustering(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Modified K Means Clustering"
        self.description = ""
        self.category = "Image Processing Tools"

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

class Mosaic(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mosaic"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MosaicWithFeathering(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mosaic With Feathering"
        self.description = ""
        self.category = "Image Processing Tools"

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

class NormalizedDifferenceVegetationIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Normalized Difference Vegetation Index"
        self.description = ""
        self.category = "Image Processing Tools"

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

class Opening(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Opening"
        self.description = ""
        self.category = "Image Processing Tools"

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

class RemoveSpurs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Spurs"
        self.description = ""
        self.category = "Image Processing Tools"

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

class Resample(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Resample"
        self.description = ""
        self.category = "Image Processing Tools"

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

class RgbToIhs(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rgb To Ihs"
        self.description = ""
        self.category = "Image Processing Tools"

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

class SplitColourComposite(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Split Colour Composite"
        self.description = ""
        self.category = "Image Processing Tools"

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

class ThickenRasterLine(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Thicken Raster Line"
        self.description = ""
        self.category = "Image Processing Tools"

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

class TophatTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tophat Transform"
        self.description = ""
        self.category = "Image Processing Tools"

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

class WriteFunctionMemoryInsertion(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Write Function Memory Insertion"
        self.description = ""
        self.category = "Image Processing Tools"

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

class AdaptiveFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Adaptive Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class BilateralFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Bilateral Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class ConservativeSmoothingFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Conservative Smoothing Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class CornerDetection(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Corner Detection"
        self.description = ""
        self.category = "Image Processing Tools"

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

class DiffOfGaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Diff Of Gaussian Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class DiversityFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Diversity Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class EdgePreservingMeanFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Preserving Mean Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class EmbossFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Emboss Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class FastAlmostGaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fast Almost Gaussian Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class GaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gaussian Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class HighPassFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "High Pass Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class HighPassMedianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "High Pass Median Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class KNearestMeanFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "K Nearest Mean Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class LaplacianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Laplacian Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class LaplacianOfGaussianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Laplacian Of Gaussian Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class LeeFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lee Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class LineDetectionFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Line Detection Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MajorityFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Majority Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MaximumFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Maximum Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MeanFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mean Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MedianFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Median Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MinimumFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Minimum Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class OlympicFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Olympic Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class PercentileFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percentile Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class PrewittFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Prewitt Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class RangeFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Range Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class RobertsCrossFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Roberts Cross Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class ScharrFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Scharr Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class SobelFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sobel Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class StandardDeviationFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Standard Deviation Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class TotalFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Total Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class UnsharpMasking(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Unsharp Masking"
        self.description = ""
        self.category = "Image Processing Tools"

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

class UserInedWeightsFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "User Ined Weights Filter"
        self.description = ""
        self.category = "Image Processing Tools"

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

class BalanceContrastEnhancement(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Balance Contrast Enhancement"
        self.description = ""
        self.category = "Image Processing Tools"

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

class CorrectVignetting(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Correct Vignetting"
        self.description = ""
        self.category = "Image Processing Tools"

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

class DirectDecorrelationStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Direct Decorrelation Stretch"
        self.description = ""
        self.category = "Image Processing Tools"

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

class GammaCorrection(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gamma Correction"
        self.description = ""
        self.category = "Image Processing Tools"

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

class GaussianContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Gaussian Contrast Stretch"
        self.description = ""
        self.category = "Image Processing Tools"

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

class HistogramEqualization(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram Equalization"
        self.description = ""
        self.category = "Image Processing Tools"

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

class HistogramMatching(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram Matching"
        self.description = ""
        self.category = "Image Processing Tools"

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

class HistogramMatchingTwoImages(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Histogram Matching Two Images"
        self.description = ""
        self.category = "Image Processing Tools"

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

class MinMaxContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min Max Contrast Stretch"
        self.description = ""
        self.category = "Image Processing Tools"

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

class PanchromaticSharpening(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Panchromatic Sharpening"
        self.description = ""
        self.category = "Image Processing Tools"

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

class PercentageContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Percentage Contrast Stretch"
        self.description = ""
        self.category = "Image Processing Tools"

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

class SigmoidalContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sigmoidal Contrast Stretch"
        self.description = ""
        self.category = "Image Processing Tools"

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

class StandardDeviationContrastStretch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Standard Deviation Contrast Stretch"
        self.description = ""
        self.category = "Image Processing Tools"

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

class ClassifyOverlapPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Classify Overlap Points"
        self.description = ""
        self.category = "LiDAR Tools"

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

class ClipLidarToPolygon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Clip Lidar To Polygon"
        self.description = ""
        self.category = "LiDAR Tools"

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

class ErasePolygonFromLidar(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Erase Polygon From Lidar"
        self.description = ""
        self.category = "LiDAR Tools"

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

class FilterLidarScanAngles(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Filter Lidar Scan Angles"
        self.description = ""
        self.category = "LiDAR Tools"

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

class FindFlightlineEdgePoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Flightline Edge Points"
        self.description = ""
        self.category = "LiDAR Tools"

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

class FlightlineOverlap(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flightline Overlap"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LasToAscii(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Las To Ascii"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LasToMultipointShapefile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Las To Multipoint Shapefile"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LasToShapefile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Las To Shapefile"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarBlockMaximum(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Block Maximum"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarBlockMinimum(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Block Minimum"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarClassifySubset(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Classify Subset"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarColourize(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Colourize"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarConstructVectorTin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Construct Vector Tin"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarElevationSlice(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Elevation Slice"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarGroundPointFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Ground Point Filter"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarHexBinning(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Hex Binning"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarHillshade(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Hillshade"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarHistogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Histogram"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarIdwInterpolation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Idw Interpolation"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarInfo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Info"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarJoin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Join"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarKappaIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Kappa Index"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarNearestNeighbourGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Nearest Neighbour Gridding"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarPointDensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Point Density"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarPointStats(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Point Stats"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarRemoveDuplicates(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Remove Duplicates"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarRemoveOutliers(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Remove Outliers"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarSegmentation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Segmentation"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarSegmentationBasedFilter(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Segmentation Based Filter"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarThin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Thin"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarThinHighDensity(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Thin High Density"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarTile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tile"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarTileFootprint(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tile Footprint"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarTinGridding(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tin Gridding"
        self.description = ""
        self.category = "LiDAR Tools"

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

class LidarTophatTransform(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Lidar Tophat Transform"
        self.description = ""
        self.category = "LiDAR Tools"

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

class NormalVectors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Normal Vectors"
        self.description = ""
        self.category = "LiDAR Tools"

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

class SelectTilesByPolygon(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Select Tiles By Polygon"
        self.description = ""
        self.category = "LiDAR Tools"

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

class And(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "And"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Not(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Not"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Or(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Or"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class AbsoluteValue(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Absolute Value"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Add(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Add"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Anova(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Anova"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ArcCos(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Arc Cos"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ArcSin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Arc Sin"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ArcTan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Arc Tan"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Atan2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Atan2"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class AttributeCorrelation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Attribute Correlation"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class AttributeHistogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Attribute Histogram"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class AttributeScattergram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Attribute Scattergram"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Ceil(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ceil"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Cos(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cos"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Cosh(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cosh"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class CrispnessIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Crispness Index"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class CrossTabulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cross Tabulation"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class CumulativeDistribution(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Cumulative Distribution"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Decrement(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Decrement"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Divide(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Divide"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class EqualTo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Equal To"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Exp(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Exp"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Exp2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Exp2"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ExtractRasterStatistics(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Raster Statistics"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Floor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Floor"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class GreaterThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Greater Than"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ImageAutocorrelation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Autocorrelation"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ImageCorrelation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Correlation"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ImageRegression(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Image Regression"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class InPlaceAdd(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Add"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class InPlaceDivide(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Divide"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class InPlaceMultiply(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Multiply"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class InPlaceSubtract(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "In Place Subtract"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Increment(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Increment"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class IntegerDivision(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Integer Division"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class IsNoData(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Is No Data"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class KappaIndex(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Kappa Index"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class KsTestForNormality(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ks Test For Normality"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class LessThan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Less Than"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ListUniqueValues(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "List Unique Values"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Ln(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Ln"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Log10(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Log10"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Log2(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Log2"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Max(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Max"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Min(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Min"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Modulo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Modulo"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Multiply(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Multiply"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Negate(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Negate"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class NotEqualTo(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Not Equal To"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Power(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Power"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class PrincipalComponentAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Principal Component Analysis"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Quantiles(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Quantiles"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class RandomField(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Random Field"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class RandomSample(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Random Sample"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class RasterHistogram(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Histogram"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class RasterSummaryStats(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Summary Stats"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Reciprocal(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Reciprocal"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class RescaleValueRange(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rescale Value Range"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class RootMeanSquareError(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Root Mean Square Error"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Round(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Round"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Sin(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sin"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Sinh(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Sinh"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Square(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Square"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class SquareRoot(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Square Root"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Subtract(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Subtract"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Tan(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tan"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Tanh(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tanh"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ToDegrees(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "To Degrees"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ToRadians(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "To Radians"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class TrendSurface(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trend Surface"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class TrendSurfaceVectorPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trend Surface Vector Points"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Truncate(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Truncate"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class TurningBandsSimulation(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Turning Bands Simulation"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class Xor(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Xor"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class ZScores(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Z Scores"
        self.description = ""
        self.category = "Math and Stats Tools"

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

class DistanceToOutlet(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Distance To Outlet"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class ExtractStreams(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Streams"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class ExtractValleys(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Extract Valleys"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class FarthestChannelHead(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Farthest Channel Head"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class FindMainStem(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Find Main Stem"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class HackStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Hack Stream Order"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class HortonStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Horton Stream Order"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class LengthOfUpstreamChannels(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Length Of Upstream Channels"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class LongProfile(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Long Profile"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class LongProfileFromPoints(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Long Profile From Points"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class RasterStreamsToVector(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Streams To Vector"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class RasterizeStreams(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rasterize Streams"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class RemoveShortStreams(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Remove Short Streams"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class ShreveStreamMagnitude(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shreve Stream Magnitude"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class StrahlerStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Strahler Stream Order"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class StreamLinkClass(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Class"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class StreamLinkIdentifier(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Identifier"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class StreamLinkLength(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Length"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class StreamLinkSlope(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Link Slope"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class StreamSlopeContinuous(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Stream Slope Continuous"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class TopologicalStreamOrder(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Topological Stream Order"
        self.description = ""
        self.category = "Stream Network Analysis"

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

class TributaryIdentifier(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tributary Identifier"
        self.description = ""
        self.category = "Stream Network Analysis"

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


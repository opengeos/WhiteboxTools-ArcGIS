import arcpy
from WBT.whitebox_tools import WhiteboxTools
wbt = WhiteboxTools()
wbt.set_verbose_mode(True)

tool_labels = []
tool_labels.append("Absolute Value")
tool_labels.append("Adaptive Filter")
tool_labels.append("Add")
tool_labels.append("Add Point Coordinates To Table")
tool_labels.append("Aggregate Raster")
tool_labels.append("And")
tool_labels.append("Anova")
tool_labels.append("Arc Cos")
tool_labels.append("Arc Sin")
tool_labels.append("Arc Tan")
tool_labels.append("Aspect")
tool_labels.append("Atan2")
tool_labels.append("Attribute Correlation")
tool_labels.append("Attribute Histogram")
tool_labels.append("Attribute Scattergram")
tool_labels.append("Average Flowpath Slope")
tool_labels.append("Average Overlay")
tool_labels.append("Average Upslope Flowpath Length")
tool_labels.append("Balance Contrast Enhancement")
tool_labels.append("Basins")
tool_labels.append("Bilateral Filter")
tool_labels.append("Block Maximum Gridding")
tool_labels.append("Block Minimum Gridding")
tool_labels.append("Breach Depressions")
tool_labels.append("Breach Single Cell Pits")
tool_labels.append("Buffer Raster")
tool_labels.append("Ceil")
tool_labels.append("Centroid")
tool_labels.append("Centroid Vector")
tool_labels.append("Change Vector Analysis")
tool_labels.append("Circular Variance Of Aspect")
tool_labels.append("Classify Overlap Points")
tool_labels.append("Clip")
tool_labels.append("Clip Lidar To Polygon")
tool_labels.append("Clip Raster To Polygon")
tool_labels.append("Closing")
tool_labels.append("Clump")
tool_labels.append("Compactness Ratio")
tool_labels.append("Conservative Smoothing Filter")
tool_labels.append("Construct Vector Tin")
tool_labels.append("Convert Nodata To Zero")
tool_labels.append("Convert Raster Format")
tool_labels.append("Corner Detection")
tool_labels.append("Correct Vignetting")
tool_labels.append("Cos")
tool_labels.append("Cosh")
tool_labels.append("Cost Allocation")
tool_labels.append("Cost Distance")
tool_labels.append("Cost Pathway")
tool_labels.append("Count If")
tool_labels.append("Create Colour Composite")
tool_labels.append("Create Hexagonal Vector Grid")
tool_labels.append("Create Plane")
tool_labels.append("Create Rectangular Vector Grid")
tool_labels.append("Crispness Index")
tool_labels.append("Cross Tabulation")
tool_labels.append("Cumulative Distribution")
tool_labels.append("D Inf Flow Accumulation")
tool_labels.append("D Inf Mass Flux")
tool_labels.append("D Inf Pointer")
tool_labels.append("D8 Flow Accumulation")
tool_labels.append("D8 Mass Flux")
tool_labels.append("D8 Pointer")
tool_labels.append("Decrement")
tool_labels.append("Depth In Sink")
tool_labels.append("Dev From Mean Elev")
tool_labels.append("Diff From Mean Elev")
tool_labels.append("Diff Of Gaussian Filter")
tool_labels.append("Difference")
tool_labels.append("Direct Decorrelation Stretch")
tool_labels.append("Directional Relief")
tool_labels.append("Dissolve")
tool_labels.append("Distance To Outlet")
tool_labels.append("Diversity Filter")
tool_labels.append("Divide")
tool_labels.append("Downslope Distance To Stream")
tool_labels.append("Downslope Flowpath Length")
tool_labels.append("Downslope Index")
tool_labels.append("Drainage Preserving Smoothing")
tool_labels.append("Edge Density")
tool_labels.append("Edge Preserving Mean Filter")
tool_labels.append("Edge Proportion")
tool_labels.append("Elev Above Pit")
tool_labels.append("Elev Percentile")
tool_labels.append("Elev Relative To Min Max")
tool_labels.append("Elev Relative To Watershed Min Max")
tool_labels.append("Elevation Above Stream")
tool_labels.append("Elevation Above Stream Euclidean")
tool_labels.append("Eliminate Coincident Points")
tool_labels.append("Elongation Ratio")
tool_labels.append("Emboss Filter")
tool_labels.append("Equal To")
tool_labels.append("Erase")
tool_labels.append("Erase Polygon From Lidar")
tool_labels.append("Erase Polygon From Raster")
tool_labels.append("Euclidean Allocation")
tool_labels.append("Euclidean Distance")
tool_labels.append("Exp")
tool_labels.append("Exp2")
tool_labels.append("Export Table To Csv")
tool_labels.append("Extend Vector Lines")
tool_labels.append("Extract Nodes")
tool_labels.append("Extract Raster Statistics")
tool_labels.append("Extract Raster Values At Points")
tool_labels.append("Extract Streams")
tool_labels.append("Extract Valleys")
tool_labels.append("Farthest Channel Head")
tool_labels.append("Fast Almost Gaussian Filter")
tool_labels.append("Fd8 Flow Accumulation")
tool_labels.append("Fd8 Pointer")
tool_labels.append("Feature Preserving Denoise")
tool_labels.append("Fetch Analysis")
tool_labels.append("Fill Burn")
tool_labels.append("Fill Depressions")
tool_labels.append("Fill Missing Data")
tool_labels.append("Fill Single Cell Pits")
tool_labels.append("Filter Lidar Scan Angles")
tool_labels.append("Find Flightline Edge Points")
tool_labels.append("Find Lowest Or Highest Points")
tool_labels.append("Find Main Stem")
tool_labels.append("Find No Flow Cells")
tool_labels.append("Find Parallel Flow")
tool_labels.append("Find Patch Or Class Edge Cells")
tool_labels.append("Find Ridges")
tool_labels.append("Flatten Lakes")
tool_labels.append("Flightline Overlap")
tool_labels.append("Flip Image")
tool_labels.append("Flood Order")
tool_labels.append("Floor")
tool_labels.append("Flow Accumulation Full Workflow")
tool_labels.append("Flow Length Diff")
tool_labels.append("Gamma Correction")
tool_labels.append("Gaussian Contrast Stretch")
tool_labels.append("Gaussian Filter")
tool_labels.append("Greater Than")
tool_labels.append("Hack Stream Order")
tool_labels.append("High Pass Filter")
tool_labels.append("High Pass Median Filter")
tool_labels.append("Highest Position")
tool_labels.append("Hillshade")
tool_labels.append("Hillslopes")
tool_labels.append("Histogram Equalization")
tool_labels.append("Histogram Matching")
tool_labels.append("Histogram Matching Two Images")
tool_labels.append("Hole Proportion")
tool_labels.append("Horizon Angle")
tool_labels.append("Horton Stream Order")
tool_labels.append("Hypsometric Analysis")
tool_labels.append("Idw Interpolation")
tool_labels.append("Ihs To Rgb")
tool_labels.append("Image Autocorrelation")
tool_labels.append("Image Correlation")
tool_labels.append("Image Regression")
tool_labels.append("Image Stack Profile")
tool_labels.append("Impoundment Size Index")
tool_labels.append("In Place Add")
tool_labels.append("In Place Divide")
tool_labels.append("In Place Multiply")
tool_labels.append("In Place Subtract")
tool_labels.append("Increment")
tool_labels.append("Integer Division")
tool_labels.append("Integral Image")
tool_labels.append("Intersect")
tool_labels.append("Is No Data")
tool_labels.append("Isobasins")
tool_labels.append("Jenson Snap Pour Points")
tool_labels.append("Join Tables")
tool_labels.append("K Means Clustering")
tool_labels.append("K Nearest Mean Filter")
tool_labels.append("Kappa Index")
tool_labels.append("Ks Test For Normality")
tool_labels.append("Laplacian Filter")
tool_labels.append("Laplacian Of Gaussian Filter")
tool_labels.append("Las To Ascii")
tool_labels.append("Las To Multipoint Shapefile")
tool_labels.append("Las To Shapefile")
tool_labels.append("Layer Footprint")
tool_labels.append("Lee Filter")
tool_labels.append("Length Of Upstream Channels")
tool_labels.append("Less Than")
tool_labels.append("Lidar Block Maximum")
tool_labels.append("Lidar Block Minimum")
tool_labels.append("Lidar Classify Subset")
tool_labels.append("Lidar Colourize")
tool_labels.append("Lidar Construct Vector Tin")
tool_labels.append("Lidar Elevation Slice")
tool_labels.append("Lidar Ground Point Filter")
tool_labels.append("Lidar Hex Binning")
tool_labels.append("Lidar Hillshade")
tool_labels.append("Lidar Histogram")
tool_labels.append("Lidar Idw Interpolation")
tool_labels.append("Lidar Info")
tool_labels.append("Lidar Join")
tool_labels.append("Lidar Kappa Index")
tool_labels.append("Lidar Nearest Neighbour Gridding")
tool_labels.append("Lidar Point Density")
tool_labels.append("Lidar Point Stats")
tool_labels.append("Lidar Remove Duplicates")
tool_labels.append("Lidar Remove Outliers")
tool_labels.append("Lidar Segmentation")
tool_labels.append("Lidar Segmentation Based Filter")
tool_labels.append("Lidar Thin")
tool_labels.append("Lidar Thin High Density")
tool_labels.append("Lidar Tile")
tool_labels.append("Lidar Tile Footprint")
tool_labels.append("Lidar Tin Gridding")
tool_labels.append("Lidar Tophat Transform")
tool_labels.append("Line Detection Filter")
tool_labels.append("Line Intersections")
tool_labels.append("Line Thinning")
tool_labels.append("Linearity Index")
tool_labels.append("Lines To Polygons")
tool_labels.append("List Unique Values")
tool_labels.append("Ln")
tool_labels.append("Log10")
tool_labels.append("Log2")
tool_labels.append("Long Profile")
tool_labels.append("Long Profile From Points")
tool_labels.append("Longest Flowpath")
tool_labels.append("Lowest Position")
tool_labels.append("Majority Filter")
tool_labels.append("Max")
tool_labels.append("Max Absolute Overlay")
tool_labels.append("Max Anisotropy Dev")
tool_labels.append("Max Anisotropy Dev Signature")
tool_labels.append("Max Branch Length")
tool_labels.append("Max Difference From Mean")
tool_labels.append("Max Downslope Elev Change")
tool_labels.append("Max Elev Dev Signature")
tool_labels.append("Max Elevation Deviation")
tool_labels.append("Max Overlay")
tool_labels.append("Max Upslope Flowpath Length")
tool_labels.append("Maximum Filter")
tool_labels.append("Mean Filter")
tool_labels.append("Median Filter")
tool_labels.append("Medoid")
tool_labels.append("Merge Table With Csv")
tool_labels.append("Merge Vectors")
tool_labels.append("Min")
tool_labels.append("Min Absolute Overlay")
tool_labels.append("Min Downslope Elev Change")
tool_labels.append("Min Max Contrast Stretch")
tool_labels.append("Min Overlay")
tool_labels.append("Minimum Bounding Box")
tool_labels.append("Minimum Bounding Circle")
tool_labels.append("Minimum Bounding Envelope")
tool_labels.append("Minimum Convex Hull")
tool_labels.append("Minimum Filter")
tool_labels.append("Modified K Means Clustering")
tool_labels.append("Modulo")
tool_labels.append("Mosaic")
tool_labels.append("Mosaic With Feathering")
tool_labels.append("Multi Part To Single Part")
tool_labels.append("Multiply")
tool_labels.append("Multiscale Roughness")
tool_labels.append("Multiscale Roughness Signature")
tool_labels.append("Multiscale Topographic Position Image")
tool_labels.append("Nearest Neighbour Gridding")
tool_labels.append("Negate")
tool_labels.append("New Raster From Base")
tool_labels.append("Normal Vectors")
tool_labels.append("Normalized Difference Vegetation Index")
tool_labels.append("Not")
tool_labels.append("Not Equal To")
tool_labels.append("Num Downslope Neighbours")
tool_labels.append("Num Inflowing Neighbours")
tool_labels.append("Num Upslope Neighbours")
tool_labels.append("Olympic Filter")
tool_labels.append("Opening")
tool_labels.append("Or")
tool_labels.append("Panchromatic Sharpening")
tool_labels.append("Patch Orientation")
tool_labels.append("Pennock Landform Class")
tool_labels.append("Percent Elev Range")
tool_labels.append("Percent Equal To")
tool_labels.append("Percent Greater Than")
tool_labels.append("Percent Less Than")
tool_labels.append("Percentage Contrast Stretch")
tool_labels.append("Percentile Filter")
tool_labels.append("Perimeter Area Ratio")
tool_labels.append("Pick From List")
tool_labels.append("Plan Curvature")
tool_labels.append("Polygon Area")
tool_labels.append("Polygon Long Axis")
tool_labels.append("Polygon Perimeter")
tool_labels.append("Polygon Short Axis")
tool_labels.append("Polygonize")
tool_labels.append("Polygons To Lines")
tool_labels.append("Power")
tool_labels.append("Prewitt Filter")
tool_labels.append("Principal Component Analysis")
tool_labels.append("Print Geo Tiff Tags")
tool_labels.append("Profile")
tool_labels.append("Profile Curvature")
tool_labels.append("Quantiles")
tool_labels.append("Radius Of Gyration")
tool_labels.append("Raise Walls")
tool_labels.append("Random Field")
tool_labels.append("Random Sample")
tool_labels.append("Range Filter")
tool_labels.append("Raster Area")
tool_labels.append("Raster Cell Assignment")
tool_labels.append("Raster Histogram")
tool_labels.append("Raster Streams To Vector")
tool_labels.append("Raster Summary Stats")
tool_labels.append("Raster To Vector Lines")
tool_labels.append("Raster To Vector Points")
tool_labels.append("Rasterize Streams")
tool_labels.append("Reciprocal")
tool_labels.append("Reclass")
tool_labels.append("Reclass Equal Interval")
tool_labels.append("Reclass From File")
tool_labels.append("Reinitialize Attribute Table")
tool_labels.append("Related Circumscribing Circle")
tool_labels.append("Relative Aspect")
tool_labels.append("Relative Stream Power Index")
tool_labels.append("Relative Topographic Position")
tool_labels.append("Remove Off Terrain Objects")
tool_labels.append("Remove Polygon Holes")
tool_labels.append("Remove Short Streams")
tool_labels.append("Remove Spurs")
tool_labels.append("Resample")
tool_labels.append("Rescale Value Range")
tool_labels.append("Rgb To Ihs")
tool_labels.append("Rho8 Pointer")
tool_labels.append("Roberts Cross Filter")
tool_labels.append("Root Mean Square Error")
tool_labels.append("Round")
tool_labels.append("Ruggedness Index")
tool_labels.append("Scharr Filter")
tool_labels.append("Sediment Transport Index")
tool_labels.append("Select Tiles By Polygon")
tool_labels.append("Set Nodata Value")
tool_labels.append("Shape Complexity Index")
tool_labels.append("Shreve Stream Magnitude")
tool_labels.append("Sigmoidal Contrast Stretch")
tool_labels.append("Sin")
tool_labels.append("Single Part To Multi Part")
tool_labels.append("Sinh")
tool_labels.append("Sink")
tool_labels.append("Slope")
tool_labels.append("Slope Vs Elevation Plot")
tool_labels.append("Smooth Vectors")
tool_labels.append("Snap Pour Points")
tool_labels.append("Sobel Filter")
tool_labels.append("Split Colour Composite")
tool_labels.append("Split With Lines")
tool_labels.append("Square")
tool_labels.append("Square Root")
tool_labels.append("Standard Deviation Contrast Stretch")
tool_labels.append("Standard Deviation Filter")
tool_labels.append("Standard Deviation Of Slope")
tool_labels.append("Stochastic Depression Analysis")
tool_labels.append("Strahler Order Basins")
tool_labels.append("Strahler Stream Order")
tool_labels.append("Stream Link Class")
tool_labels.append("Stream Link Identifier")
tool_labels.append("Stream Link Length")
tool_labels.append("Stream Link Slope")
tool_labels.append("Stream Slope Continuous")
tool_labels.append("Subbasins")
tool_labels.append("Subtract")
tool_labels.append("Sum Overlay")
tool_labels.append("Surface Area Ratio")
tool_labels.append("Symmetrical Difference")
tool_labels.append("Tan")
tool_labels.append("Tangential Curvature")
tool_labels.append("Tanh")
tool_labels.append("Thicken Raster Line")
tool_labels.append("Tin Gridding")
tool_labels.append("To Degrees")
tool_labels.append("To Radians")
tool_labels.append("Tophat Transform")
tool_labels.append("Topological Stream Order")
tool_labels.append("Total Curvature")
tool_labels.append("Total Filter")
tool_labels.append("Trace Downslope Flowpaths")
tool_labels.append("Trend Surface")
tool_labels.append("Trend Surface Vector Points")
tool_labels.append("Tributary Identifier")
tool_labels.append("Truncate")
tool_labels.append("Turning Bands Simulation")
tool_labels.append("Union")
tool_labels.append("Unnest Basins")
tool_labels.append("Unsharp Masking")
tool_labels.append("User Defined Weights Filter")
tool_labels.append("Vector Hex Binning")
tool_labels.append("Vector Lines To Raster")
tool_labels.append("Vector Points To Raster")
tool_labels.append("Vector Polygons To Raster")
tool_labels.append("Viewshed")
tool_labels.append("Visibility Index")
tool_labels.append("Voronoi Diagram")
tool_labels.append("Watershed")
tool_labels.append("Weighted Overlay")
tool_labels.append("Weighted Sum")
tool_labels.append("Wetness Index")
tool_labels.append("Write Function Memory Insertion")
tool_labels.append("Xor")
tool_labels.append("Z Scores")


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
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
        tools.append(RasterArea)
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
        tools.append(ImpoundmentSizeIndex)
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
        tools.append(UserDefinedWeightsFilter)
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
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
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
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
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
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
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
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        # First parameter
        param0 = arcpy.Parameter(
            displayName="Keywords",
            name="keywords",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        # param0.multiValue = True
        param0.value = "lidar"
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
            messages.addMessage("{}. {}: {}".format(index + 1, tool, tools[tool]))
        return


class ToolHelp(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Help"
        self.description = "Help description for a specific tool."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

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
        messages.addMessage(wbt.tool_help(tool_name))
        return


class ToolParameters(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Tool Parameters"
        self.description = "Tool parameter descriptions for a specific tool."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

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
        messages.addMessage(wbt.tool_parameters(tool_name))
        return


class ViewCode(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "View Code"
        self.description = "Source code for a specific tool."
        self.category = "About WhiteboxTools"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

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
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        tool_name = arcpy.Parameter(
            displayName="Select a tool",
            name="tool_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        tool_name.value = "Lidar Info"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

        args = arcpy.Parameter(
            displayName="Arguments",
            name="agrs",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        params = [tool_name, args]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        param0 = parameters[0].valueAsText
        args = parameters[1].valueAsText
        tool_name = param0.replace(" ", "").strip()
        messages.addMessage(wbt.run_tool(tool_name, args))
        return


class AddPointCoordinatesToTable(object):
    def __init__(self):
        self.label = "Add Point Coordinates To Table"
        self.description = "Modifies the attribute table of a point vector by adding fields containing each point's X and Y coordinates."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.add_point_coordinates_to_table(input))
        return


class ConvertNodataToZero(object):
    def __init__(self):
        self.label = "Convert Nodata To Zero"
        self.description = "Converts nodata values in a raster to zero."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.convert_nodata_to_zero(input, output))
        return


class ConvertRasterFormat(object):
    def __init__(self):
        self.label = "Convert Raster Format"
        self.description = "Converts raster data from one format to another."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.convert_raster_format(input, output))
        return


class ExportTableToCsv(object):
    def __init__(self):
        self.label = "Export Table To Csv"
        self.description = "Exports an attribute table to a CSV text file."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["csv"]

        headers = arcpy.Parameter(
            displayName="Export field names as file header?",
            name="headers",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        headers.value = "true"

        params = [input, output, headers]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        headers = parameters[2].valueAsText
        messages.addMessage(wbt.export_table_to_csv(input, output, headers))
        return


class JoinTables(object):
    def __init__(self):
        self.label = "Join Tables"
        self.description = "Merge a vector's attribute table with another table based on a common field."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Primary Vector File",
            name="input1",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        pkey = arcpy.Parameter(
            displayName="Primary Key Field",
            name="pkey",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        pkey.parameterDependencies = [input1.name]

        input2 = arcpy.Parameter(
            displayName="Input Foreign Vector File",
            name="input2",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        fkey = arcpy.Parameter(
            displayName="Foreign Key Field",
            name="fkey",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        fkey.parameterDependencies = [input2.name]

        import_field = arcpy.Parameter(
            displayName="Imported Field",
            name="import_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        import_field.parameterDependencies = [input2.name]

        params = [input1, pkey, input2, fkey, import_field]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        pkey = parameters[1].valueAsText
        input2 = parameters[2].valueAsText
        fkey = parameters[3].valueAsText
        import_field = parameters[4].valueAsText
        messages.addMessage(wbt.join_tables(input1, pkey, input2, fkey, import_field))
        return


class LinesToPolygons(object):
    def __init__(self):
        self.label = "Lines To Polygons"
        self.description = "Converts vector polylines to polygons."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Line File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polyline"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.lines_to_polygons(input, output))
        return


class MergeTableWithCsv(object):
    def __init__(self):
        self.label = "Merge Table With Csv"
        self.description = "Merge a vector's attribute table with a table contained within a CSV text file."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Primary Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        pkey = arcpy.Parameter(
            displayName="Primary Key Field",
            name="pkey",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        pkey.parameterDependencies = [input.name]

        csv = arcpy.Parameter(
            displayName="Input CSV File",
            name="csv",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        csv.filter.list = ["csv"]

        fkey = arcpy.Parameter(
            displayName="Foreign Key Field",
            name="fkey",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        fkey.parameterDependencies = [csv.name]

        import_field = arcpy.Parameter(
            displayName="Imported Field",
            name="import_field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        import_field.parameterDependencies = [csv.name]

        params = [input, pkey, csv, fkey, import_field]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        pkey = parameters[1].valueAsText
        csv = parameters[2].valueAsText
        fkey = parameters[3].valueAsText
        import_field = parameters[4].valueAsText
        messages.addMessage(wbt.merge_table_with_csv(input, pkey, csv, fkey, import_field))
        return


class MergeVectors(object):
    def __init__(self):
        self.label = "Merge Vectors"
        self.description = "Combines two or more input vectors of the same ShapeType creating a single, new output vector."
        self.category = "Data Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Vector Files",
            name="inputs",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.merge_vectors(inputs, output))
        return


class MultiPartToSinglePart(object):
    def __init__(self):
        self.label = "Multi Part To Single Part"
        self.description = "Converts a vector file containing multi-part features into a vector containing only single-part features."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Line or Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Line or Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        exclude_holes = arcpy.Parameter(
            displayName="Exclude hole parts?",
            name="exclude_holes",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        exclude_holes.value = "true"

        params = [input, output, exclude_holes]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        exclude_holes = parameters[2].valueAsText
        messages.addMessage(wbt.multi_part_to_single_part(input, output, exclude_holes))
        return


class NewRasterFromBase(object):
    def __init__(self):
        self.label = "New Raster From Base"
        self.description = "Creates a new raster using a base image."
        self.category = "Data Tools"

    def getParameterInfo(self):
        base = arcpy.Parameter(
            displayName="Input Base File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        value = arcpy.Parameter(
            displayName="Constant Value",
            name="value",
            datatype=["GPString", "GPDouble"],
            parameterType="Optional",
            direction="Input")

        value.value = "nodata"

        data_type = arcpy.Parameter(
            displayName="Data Type",
            name="data_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        data_type.filter.type = "ValueList"
        data_type.filter.list = ['double', 'float', 'integer']

        data_type.value = "float"

        params = [base, output, value, data_type]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        base = parameters[0].valueAsText
        output = parameters[1].valueAsText
        value = parameters[2].valueAsText
        data_type = parameters[3].valueAsText
        messages.addMessage(wbt.new_raster_from_base(base, output, value, data_type))
        return


class PolygonsToLines(object):
    def __init__(self):
        self.label = "Polygons To Lines"
        self.description = "Converts vector polygons to polylines."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output Line File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polyline"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.polygons_to_lines(input, output))
        return


class PrintGeoTiffTags(object):
    def __init__(self):
        self.label = "Print Geo Tiff Tags"
        self.description = "Prints the tags within a GeoTIFF."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input GeoTIFF Raster File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.print_geo_tiff_tags(input))
        return


class RasterToVectorLines(object):
    def __init__(self):
        self.label = "Raster To Vector Lines"
        self.description = "Converts a raster lines features into a vector of the POLYLINE shapetype."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Raster Lines File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polyline"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.raster_to_vector_lines(input, output))
        return


class RasterToVectorPoints(object):
    def __init__(self):
        self.label = "Raster To Vector Points"
        self.description = "Converts a raster dataset to a vector of the POINT shapetype."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Raster File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Points File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.raster_to_vector_points(input, output))
        return


class ReinitializeAttributeTable(object):
    def __init__(self):
        self.label = "Reinitialize Attribute Table"
        self.description = "Reinitializes a vector's attribute table deleting all fields but the feature ID (FID)."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.reinitialize_attribute_table(input))
        return


class RemovePolygonHoles(object):
    def __init__(self):
        self.label = "Remove Polygon Holes"
        self.description = "Removes holes within the features of a vector polygon file."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.remove_polygon_holes(input, output))
        return


class SetNodataValue(object):
    def __init__(self):
        self.label = "Set Nodata Value"
        self.description = "Assign a specified value in an input image to the NoData value."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        back_value = arcpy.Parameter(
            displayName="Background Value",
            name="back_value",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        back_value.value = "0.0"

        params = [input, output, back_value]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        back_value = parameters[2].valueAsText
        messages.addMessage(wbt.set_nodata_value(input, output, back_value))
        return


class SinglePartToMultiPart(object):
    def __init__(self):
        self.label = "Single Part To Multi Part"
        self.description = "Converts a vector file containing multi-part features into a vector containing only single-part features."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Line or Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        field = arcpy.Parameter(
            displayName="Grouping ID Field Name",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        output = arcpy.Parameter(
            displayName="Output Line or Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, field, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.single_part_to_multi_part(input, field, output))
        return


class VectorLinesToRaster(object):
    def __init__(self):
        self.label = "Vector Lines To Raster"
        self.description = "Converts a vector containing polylines into a raster."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polyline"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        field.value = "FID"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        nodata = arcpy.Parameter(
            displayName="Background value is NoData?",
            name="nodata",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        nodata.value = "true"

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        params = [input, field, output, nodata, cell_size, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        nodata = parameters[3].valueAsText
        cell_size = parameters[4].valueAsText
        base = parameters[5].valueAsText
        messages.addMessage(wbt.vector_lines_to_raster(input, field, output, nodata, cell_size, base))
        return


class VectorPointsToRaster(object):
    def __init__(self):
        self.label = "Vector Points To Raster"
        self.description = "Converts a vector containing points into a raster."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        field.value = "FID"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        assign = arcpy.Parameter(
            displayName="Assignment Operation",
            name="assign",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        assign.filter.type = "ValueList"
        assign.filter.list = ['first', 'last', 'min', 'max', 'sum']

        assign.value = "last"

        nodata = arcpy.Parameter(
            displayName="Background value is NoData?",
            name="nodata",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        nodata.value = "true"

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        params = [input, field, output, assign, nodata, cell_size, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        assign = parameters[3].valueAsText
        nodata = parameters[4].valueAsText
        cell_size = parameters[5].valueAsText
        base = parameters[6].valueAsText
        messages.addMessage(wbt.vector_points_to_raster(input, field, output, assign, nodata, cell_size, base))
        return


class VectorPolygonsToRaster(object):
    def __init__(self):
        self.label = "Vector Polygons To Raster"
        self.description = "Converts a vector containing polygons into a raster."
        self.category = "Data Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        field.value = "FID"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        nodata = arcpy.Parameter(
            displayName="Background value is NoData?",
            name="nodata",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        nodata.value = "true"

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        params = [input, field, output, nodata, cell_size, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        nodata = parameters[3].valueAsText
        cell_size = parameters[4].valueAsText
        base = parameters[5].valueAsText
        messages.addMessage(wbt.vector_polygons_to_raster(input, field, output, nodata, cell_size, base))
        return


class AggregateRaster(object):
    def __init__(self):
        self.label = "Aggregate Raster"
        self.description = "Aggregates a raster to a lower resolution."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        agg_factor = arcpy.Parameter(
            displayName="Aggregation Factor (pixels)",
            name="agg_factor",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        agg_factor.value = "2"

        type = arcpy.Parameter(
            displayName="Aggregation Type",
            name="type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        type.filter.type = "ValueList"
        type.filter.list = ['mean', 'sum', 'maximum', 'minimum', 'range']

        type.value = "mean"

        params = [input, output, agg_factor, type]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        agg_factor = parameters[2].valueAsText
        type = parameters[3].valueAsText
        messages.addMessage(wbt.aggregate_raster(input, output, agg_factor, type))
        return


class BlockMaximumGridding(object):
    def __init__(self):
        self.label = "Block Maximum Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using a block maximum scheme."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        use_z = arcpy.Parameter(
            displayName="Use z-coordinate instead of field?",
            name="use_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        use_z.value = "false"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        params = [input, field, use_z, output, cell_size, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        use_z = parameters[2].valueAsText
        output = parameters[3].valueAsText
        cell_size = parameters[4].valueAsText
        base = parameters[5].valueAsText
        messages.addMessage(wbt.block_maximum_gridding(input, field, use_z, output, cell_size, base))
        return


class BlockMinimumGridding(object):
    def __init__(self):
        self.label = "Block Minimum Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using a block minimum scheme."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        use_z = arcpy.Parameter(
            displayName="Use z-coordinate instead of field?",
            name="use_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        use_z.value = "false"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        params = [input, field, use_z, output, cell_size, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        use_z = parameters[2].valueAsText
        output = parameters[3].valueAsText
        cell_size = parameters[4].valueAsText
        base = parameters[5].valueAsText
        messages.addMessage(wbt.block_minimum_gridding(input, field, use_z, output, cell_size, base))
        return


class Centroid(object):
    def __init__(self):
        self.label = "Centroid"
        self.description = "Calculates the centroid, or average location, of raster polygon objects."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        text_output = arcpy.Parameter(
            displayName="Output text?",
            name="text_output",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        params = [input, output, text_output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        text_output = parameters[2].valueAsText
        messages.addMessage(wbt.centroid(input, output, text_output))
        return


class CentroidVector(object):
    def __init__(self):
        self.label = "Centroid Vector"
        self.description = "Identifes the centroid point of a vector polyline or polygon feature or a group of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.centroid_vector(input, output))
        return


class Clump(object):
    def __init__(self):
        self.label = "Clump"
        self.description = "Groups cells that form discrete areas, assigning them unique identifiers."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        diag = arcpy.Parameter(
            displayName="Include diagonal connections?",
            name="diag",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        diag.value = "true"

        zero_back = arcpy.Parameter(
            displayName="Treat zero values as background?",
            name="zero_back",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        params = [input, output, diag, zero_back]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        diag = parameters[2].valueAsText
        zero_back = parameters[3].valueAsText
        messages.addMessage(wbt.clump(input, output, diag, zero_back))
        return


class ConstructVectorTin(object):
    def __init__(self):
        self.label = "Construct Vector Tin"
        self.description = "Creates a vector triangular irregular network (TIN) for a set of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        use_z = arcpy.Parameter(
            displayName="Use Shapefile 'z' values?",
            name="use_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        use_z.value = "false"

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, field, use_z, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        use_z = parameters[2].valueAsText
        output = parameters[3].valueAsText
        messages.addMessage(wbt.construct_vector_tin(input, field, use_z, output))
        return


class CreateHexagonalVectorGrid(object):
    def __init__(self):
        self.label = "Create Hexagonal Vector Grid"
        self.description = "Creates a hexagonal vector grid."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Base File",
            name="input",
            datatype=["DERasterDataset", "DEShapefile"],
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        width = arcpy.Parameter(
            displayName="Hexagon Width",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        orientation = arcpy.Parameter(
            displayName="Grid Orientation",
            name="orientation",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        orientation.filter.type = "ValueList"
        orientation.filter.list = ['horizontal', 'vertical']

        orientation.value = "horizontal"

        params = [input, output, width, orientation]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        width = parameters[2].valueAsText
        orientation = parameters[3].valueAsText
        messages.addMessage(wbt.create_hexagonal_vector_grid(input, output, width, orientation))
        return


class CreatePlane(object):
    def __init__(self):
        self.label = "Create Plane"
        self.description = "Creates a raster image based on the equation for a simple plane."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        base = arcpy.Parameter(
            displayName="Input Base File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        gradient = arcpy.Parameter(
            displayName="Gradient",
            name="gradient",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        gradient.value = "15.0"

        aspect = arcpy.Parameter(
            displayName="Aspect",
            name="aspect",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        aspect.value = "90.0"

        constant = arcpy.Parameter(
            displayName="Constant",
            name="constant",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        constant.value = "0.0"

        params = [base, output, gradient, aspect, constant]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        base = parameters[0].valueAsText
        output = parameters[1].valueAsText
        gradient = parameters[2].valueAsText
        aspect = parameters[3].valueAsText
        constant = parameters[4].valueAsText
        messages.addMessage(wbt.create_plane(base, output, gradient, aspect, constant))
        return


class CreateRectangularVectorGrid(object):
    def __init__(self):
        self.label = "Create Rectangular Vector Grid"
        self.description = "Creates a rectangular vector grid."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Base File",
            name="input",
            datatype=["DERasterDataset", "DEShapefile"],
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        width = arcpy.Parameter(
            displayName="Grid Cell Width",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        height = arcpy.Parameter(
            displayName="Grid Cell Height",
            name="height",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        xorig = arcpy.Parameter(
            displayName="Grid origin x-coordinate",
            name="xorig",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        xorig.value = "0"

        yorig = arcpy.Parameter(
            displayName="Grid origin y-coordinate",
            name="yorig",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        yorig.value = "0"

        params = [input, output, width, height, xorig, yorig]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        width = parameters[2].valueAsText
        height = parameters[3].valueAsText
        xorig = parameters[4].valueAsText
        yorig = parameters[5].valueAsText
        messages.addMessage(wbt.create_rectangular_vector_grid(input, output, width, height, xorig, yorig))
        return


class Dissolve(object):
    def __init__(self):
        self.label = "Dissolve"
        self.description = "Removes the interior, or shared, boundaries within a vector polygon coverage."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        field = arcpy.Parameter(
            displayName="Dissolve Field Attribute",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        snap = arcpy.Parameter(
            displayName="Snap Tolerance",
            name="snap",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        snap.value = "0.0"

        params = [input, field, output, snap]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        snap = parameters[3].valueAsText
        messages.addMessage(wbt.dissolve(input, field, output, snap))
        return


class EliminateCoincidentPoints(object):
    def __init__(self):
        self.label = "Eliminate Coincident Points"
        self.description = "Removes any coincident, or nearly coincident, points from a vector points file."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        tolerance = arcpy.Parameter(
            displayName="Distance Tolerance",
            name="tolerance",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [input, output, tolerance]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        tolerance = parameters[2].valueAsText
        messages.addMessage(wbt.eliminate_coincident_points(input, output, tolerance))
        return


class ExtendVectorLines(object):
    def __init__(self):
        self.label = "Extend Vector Lines"
        self.description = "Extends vector lines by a specified distance."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polyline"]

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polyline"]

        dist = arcpy.Parameter(
            displayName="Extend Distance",
            name="dist",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        extend = arcpy.Parameter(
            displayName="Extend Direction",
            name="extend",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        extend.filter.type = "ValueList"
        extend.filter.list = ['both ends', 'line start', 'line end']

        extend.value = "both ends"

        params = [input, output, dist, extend]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        dist = parameters[2].valueAsText
        extend = parameters[3].valueAsText
        messages.addMessage(wbt.extend_vector_lines(input, output, dist, extend))
        return


class ExtractNodes(object):
    def __init__(self):
        self.label = "Extract Nodes"
        self.description = "Converts vector lines or polygons into vertex points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Points File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.extract_nodes(input, output))
        return


class ExtractRasterValuesAtPoints(object):
    def __init__(self):
        self.label = "Extract Raster Values At Points"
        self.description = "Extracts the values of raster(s) at vector point locations."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        points = arcpy.Parameter(
            displayName="Input Points File",
            name="points",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        points.filter.list = ["Point"]

        out_text = arcpy.Parameter(
            displayName="Output text?",
            name="out_text",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        out_text.value = "false"

        params = [inputs, points, out_text]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        points = parameters[1].valueAsText
        out_text = parameters[2].valueAsText
        messages.addMessage(wbt.extract_raster_values_at_points(inputs, points, out_text))
        return


class FindLowestOrHighestPoints(object):
    def __init__(self):
        self.label = "Find Lowest Or Highest Points"
        self.description = "Locates the lowest and/or highest valued cells in a raster."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Raster File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Points File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['lowest', 'highest', 'both']

        out_type.value = "lowest"

        params = [input, output, out_type]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_type = parameters[2].valueAsText
        messages.addMessage(wbt.find_lowest_or_highest_points(input, output, out_type))
        return


class IdwInterpolation(object):
    def __init__(self):
        self.label = "Idw Interpolation"
        self.description = "Interpolates vector points into a raster surface using an inverse-distance weighted scheme."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        use_z = arcpy.Parameter(
            displayName="Use z-coordinate instead of field?",
            name="use_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        use_z.value = "false"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        weight = arcpy.Parameter(
            displayName="IDW Weight (Exponent) Value",
            name="weight",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        weight.value = "2.0"

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        min_points = arcpy.Parameter(
            displayName="Min. Number of Points",
            name="min_points",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        params = [input, field, use_z, output, weight, radius, min_points, cell_size, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        use_z = parameters[2].valueAsText
        output = parameters[3].valueAsText
        weight = parameters[4].valueAsText
        radius = parameters[5].valueAsText
        min_points = parameters[6].valueAsText
        cell_size = parameters[7].valueAsText
        base = parameters[8].valueAsText
        messages.addMessage(wbt.idw_interpolation(input, field, use_z, output, weight, radius, min_points, cell_size, base))
        return


class LayerFootprint(object):
    def __init__(self):
        self.label = "Layer Footprint"
        self.description = "Creates a vector polygon footprint of the area covered by a raster grid or vector layer."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Raster or Vector File",
            name="input",
            datatype=["DERasterDataset", "DEShapefile"],
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.layer_footprint(input, output))
        return


class Medoid(object):
    def __init__(self):
        self.label = "Medoid"
        self.description = "Calculates the medoid for a series of vector features contained in a shapefile."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.medoid(input, output))
        return


class MinimumBoundingBox(object):
    def __init__(self):
        self.label = "Minimum Bounding Box"
        self.description = "Creates a vector minimum bounding rectangle around vector features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        criterion = arcpy.Parameter(
            displayName="Minimization Criterion",
            name="criterion",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        criterion.filter.type = "ValueList"
        criterion.filter.list = ['area', 'length', 'width', 'perimeter']

        criterion.value = "area"

        features = arcpy.Parameter(
            displayName="Find bounding rectangles around each individual feature.",
            name="features",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        features.value = "true"

        params = [input, output, criterion, features]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        criterion = parameters[2].valueAsText
        features = parameters[3].valueAsText
        messages.addMessage(wbt.minimum_bounding_box(input, output, criterion, features))
        return


class MinimumBoundingCircle(object):
    def __init__(self):
        self.label = "Minimum Bounding Circle"
        self.description = "Delineates the minimum bounding circle (i.e. smallest enclosing circle) for a group of vectors."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        features = arcpy.Parameter(
            displayName="Find bounding circle around each individual feature.",
            name="features",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        features.value = "true"

        params = [input, output, features]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        features = parameters[2].valueAsText
        messages.addMessage(wbt.minimum_bounding_circle(input, output, features))
        return


class MinimumBoundingEnvelope(object):
    def __init__(self):
        self.label = "Minimum Bounding Envelope"
        self.description = "Creates a vector axis-aligned minimum bounding rectangle (envelope) around vector features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        features = arcpy.Parameter(
            displayName="Find bounding envelop around each individual feature.",
            name="features",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        features.value = "true"

        params = [input, output, features]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        features = parameters[2].valueAsText
        messages.addMessage(wbt.minimum_bounding_envelope(input, output, features))
        return


class MinimumConvexHull(object):
    def __init__(self):
        self.label = "Minimum Convex Hull"
        self.description = "Creates a vector convex polygon around vector features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        features = arcpy.Parameter(
            displayName="Find hulls around each individual feature.",
            name="features",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        features.value = "true"

        params = [input, output, features]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        features = parameters[2].valueAsText
        messages.addMessage(wbt.minimum_convex_hull(input, output, features))
        return


class NearestNeighbourGridding(object):
    def __init__(self):
        self.label = "Nearest Neighbour Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using the nearest neighbour."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        use_z = arcpy.Parameter(
            displayName="Use z-coordinate instead of field?",
            name="use_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        use_z.value = "false"

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File (optional)",
            name="base",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        max_dist = arcpy.Parameter(
            displayName="Maximum Search Distance",
            name="max_dist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, field, use_z, output, cell_size, base, max_dist]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        use_z = parameters[2].valueAsText
        output = parameters[3].valueAsText
        cell_size = parameters[4].valueAsText
        base = parameters[5].valueAsText
        max_dist = parameters[6].valueAsText
        messages.addMessage(wbt.nearest_neighbour_gridding(input, field, use_z, output, cell_size, base, max_dist))
        return


class PolygonArea(object):
    def __init__(self):
        self.label = "Polygon Area"
        self.description = "Calculates the area of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.polygon_area(input))
        return


class PolygonLongAxis(object):
    def __init__(self):
        self.label = "Polygon Long Axis"
        self.description = "This tool can be used to map the long axis of polygon features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.polygon_long_axis(input, output))
        return


class PolygonPerimeter(object):
    def __init__(self):
        self.label = "Polygon Perimeter"
        self.description = "Calculates the perimeter of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.polygon_perimeter(input))
        return


class PolygonShortAxis(object):
    def __init__(self):
        self.label = "Polygon Short Axis"
        self.description = "This tool can be used to map the short axis of polygon features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.polygon_short_axis(input, output))
        return


class RasterArea(object):
    def __init__(self):
        self.label = "Raster Area"
        self.description = "Calculates the area of polygons or classes within a raster image."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_text = arcpy.Parameter(
            displayName="Output text?",
            name="out_text",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        units = arcpy.Parameter(
            displayName="Units",
            name="units",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        units.filter.type = "ValueList"
        units.filter.list = ['grid cells', 'map units']

        units.value = "grid cells"

        zero_back = arcpy.Parameter(
            displayName="Treat zero values as background?",
            name="zero_back",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        params = [input, output, out_text, units, zero_back]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_text = parameters[2].valueAsText
        units = parameters[3].valueAsText
        zero_back = parameters[4].valueAsText
        messages.addMessage(wbt.raster_area(input, output, out_text, units, zero_back))
        return


class RasterCellAssignment(object):
    def __init__(self):
        self.label = "Raster Cell Assignment"
        self.description = "Assign row or column number to cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        assign = arcpy.Parameter(
            displayName="Which spatial variable should be assigned?",
            name="assign",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        assign.filter.type = "ValueList"
        assign.filter.list = ['column', 'row', 'x', 'y']

        assign.value = "column"

        params = [input, output, assign]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        assign = parameters[2].valueAsText
        messages.addMessage(wbt.raster_cell_assignment(input, output, assign))
        return


class Reclass(object):
    def __init__(self):
        self.label = "Reclass"
        self.description = "Reclassifies the values in a raster image."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        reclass_vals = arcpy.Parameter(
            displayName="Class Interval Size",
            name="reclass_vals",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        assign_mode = arcpy.Parameter(
            displayName="Operate in assign mode? (i.e. Reclass data are pair values rather than triplets)",
            name="assign_mode",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input, output, reclass_vals, assign_mode]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        reclass_vals = parameters[2].valueAsText
        assign_mode = parameters[3].valueAsText
        messages.addMessage(wbt.reclass(input, output, reclass_vals, assign_mode))
        return


class ReclassEqualInterval(object):
    def __init__(self):
        self.label = "Reclass Equal Interval"
        self.description = "Reclassifies the values in a raster image based on equal-ranges."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        interval = arcpy.Parameter(
            displayName="Class Interval Size",
            name="interval",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        interval.value = "10.0"

        start_val = arcpy.Parameter(
            displayName="Starting Value",
            name="start_val",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        end_val = arcpy.Parameter(
            displayName="Ending Value",
            name="end_val",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, interval, start_val, end_val]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        interval = parameters[2].valueAsText
        start_val = parameters[3].valueAsText
        end_val = parameters[4].valueAsText
        messages.addMessage(wbt.reclass_equal_interval(input, output, interval, start_val, end_val))
        return


class ReclassFromFile(object):
    def __init__(self):
        self.label = "Reclass From File"
        self.description = "Reclassifies the values in a raster image using reclass ranges in a text file."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        reclass_file = arcpy.Parameter(
            displayName="Input Reclass Text File",
            name="reclass_file",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, reclass_file, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        reclass_file = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.reclass_from_file(input, reclass_file, output))
        return


class SmoothVectors(object):
    def __init__(self):
        self.label = "Smooth Vectors"
        self.description = "Smooths a vector coverage of either a POLYLINE or POLYGON base ShapeType."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        filter = arcpy.Parameter(
            displayName="Filter Size",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "3"

        params = [input, output, filter]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        messages.addMessage(wbt.smooth_vectors(input, output, filter))
        return


class TinGridding(object):
    def __init__(self):
        self.label = "Tin Gridding"
        self.description = "Creates a raster grid based on a triangular irregular network (TIN) fitted to vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Optional",
            direction="Input")
        field.parameterDependencies = [input.name]

        use_z = arcpy.Parameter(
            displayName="Use Shapefile 'z' values?",
            name="use_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        use_z.value = "false"

        output = arcpy.Parameter(
            displayName="Output Raster File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [input, field, use_z, output, resolution]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        use_z = parameters[2].valueAsText
        output = parameters[3].valueAsText
        resolution = parameters[4].valueAsText
        messages.addMessage(wbt.tin_gridding(input, field, use_z, output, resolution))
        return


class VectorHexBinning(object):
    def __init__(self):
        self.label = "Vector Hex Binning"
        self.description = "Hex-bins a set of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Base File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        width = arcpy.Parameter(
            displayName="Hexagon Width",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        orientation = arcpy.Parameter(
            displayName="Grid Orientation",
            name="orientation",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        orientation.filter.type = "ValueList"
        orientation.filter.list = ['horizontal', 'vertical']

        orientation.value = "horizontal"

        params = [input, output, width, orientation]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        width = parameters[2].valueAsText
        orientation = parameters[3].valueAsText
        messages.addMessage(wbt.vector_hex_binning(input, output, width, orientation))
        return


class VoronoiDiagram(object):
    def __init__(self):
        self.label = "Voronoi Diagram"
        self.description = "Creates a vector Voronoi diagram for a set of vector points."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.voronoi_diagram(input, output))
        return


class BufferRaster(object):
    def __init__(self):
        self.label = "Buffer Raster"
        self.description = "Maps a distance-based buffer around each non-background (non-zero/non-nodata) grid cell in an input image."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        size = arcpy.Parameter(
            displayName="Buffer Size",
            name="size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        gridcells = arcpy.Parameter(
            displayName="Buffer size measured in grid cells?",
            name="gridcells",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input, output, size, gridcells]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        size = parameters[2].valueAsText
        gridcells = parameters[3].valueAsText
        messages.addMessage(wbt.buffer_raster(input, output, size, gridcells))
        return


class CostAllocation(object):
    def __init__(self):
        self.label = "Cost Allocation"
        self.description = "Identifies the source cell to which each grid cell is connected by a least-cost pathway in a cost-distance analysis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        source = arcpy.Parameter(
            displayName="Input Source File",
            name="source",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        backlink = arcpy.Parameter(
            displayName="Input Backlink File",
            name="backlink",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [source, backlink, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        source = parameters[0].valueAsText
        backlink = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.cost_allocation(source, backlink, output))
        return


class CostDistance(object):
    def __init__(self):
        self.label = "Cost Distance"
        self.description = "Performs cost-distance accumulation on a cost surface and a group of source cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        source = arcpy.Parameter(
            displayName="Input Source File",
            name="source",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        cost = arcpy.Parameter(
            displayName="Input Cost (Friction) File",
            name="cost",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        out_accum = arcpy.Parameter(
            displayName="Output Cost Accumulation File",
            name="out_accum",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_accum.filter.list = ["tif"]

        out_backlink = arcpy.Parameter(
            displayName="Output Backlink File",
            name="out_backlink",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_backlink.filter.list = ["tif"]

        params = [source, cost, out_accum, out_backlink]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        source = parameters[0].valueAsText
        cost = parameters[1].valueAsText
        out_accum = parameters[2].valueAsText
        out_backlink = parameters[3].valueAsText
        messages.addMessage(wbt.cost_distance(source, cost, out_accum, out_backlink))
        return


class CostPathway(object):
    def __init__(self):
        self.label = "Cost Pathway"
        self.description = "Performs cost-distance pathway analysis using a series of destination grid cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        destination = arcpy.Parameter(
            displayName="Input Destination File",
            name="destination",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        backlink = arcpy.Parameter(
            displayName="Input Backlink File",
            name="backlink",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zero_background = arcpy.Parameter(
            displayName="Treat zero values as background?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        params = [destination, backlink, output, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        destination = parameters[0].valueAsText
        backlink = parameters[1].valueAsText
        output = parameters[2].valueAsText
        zero_background = parameters[3].valueAsText
        messages.addMessage(wbt.cost_pathway(destination, backlink, output, zero_background))
        return


class EuclideanAllocation(object):
    def __init__(self):
        self.label = "Euclidean Allocation"
        self.description = "Assigns grid cells in the output raster the value of the nearest target cell in the input image, measured by the Shih and Wu (2004) Euclidean distance transform."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.euclidean_allocation(input, output))
        return


class EuclideanDistance(object):
    def __init__(self):
        self.label = "Euclidean Distance"
        self.description = "Calculates the Shih and Wu (2004) Euclidean distance transform."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.euclidean_distance(input, output))
        return


class AverageOverlay(object):
    def __init__(self):
        self.label = "Average Overlay"
        self.description = "Calculates the average for each grid cell from a group of raster images."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.average_overlay(inputs, output))
        return


class Clip(object):
    def __init__(self):
        self.label = "Clip"
        self.description = "Extract all the features, or parts of features, that overlap with the features of the clip vector."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        clip = arcpy.Parameter(
            displayName="Input Clip Polygon Vector File",
            name="clip",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        clip.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, clip, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        clip = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.clip(input, clip, output))
        return


class ClipRasterToPolygon(object):
    def __init__(self):
        self.label = "Clip Raster To Polygon"
        self.description = "Clips a raster to a vector polygon."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        polygons = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="polygons",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        polygons.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        maintain_dimensions = arcpy.Parameter(
            displayName="Maintain input raster dimensions?",
            name="maintain_dimensions",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        maintain_dimensions.value = "false"

        params = [input, polygons, output, maintain_dimensions]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        polygons = parameters[1].valueAsText
        output = parameters[2].valueAsText
        maintain_dimensions = parameters[3].valueAsText
        messages.addMessage(wbt.clip_raster_to_polygon(input, polygons, output, maintain_dimensions))
        return


class CountIf(object):
    def __init__(self):
        self.label = "Count If"
        self.description = "Counts the number of occurrences of a specified value in a cell-stack of rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        value = arcpy.Parameter(
            displayName="Value",
            name="value",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [inputs, output, value]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        value = parameters[2].valueAsText
        messages.addMessage(wbt.count_if(inputs, output, value))
        return


class Difference(object):
    def __init__(self):
        self.label = "Difference"
        self.description = "Outputs the features that occur in one of the two vector inputs but not both, i.e. no overlapping features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        overlay = arcpy.Parameter(
            displayName="Input Overlay Vector File",
            name="overlay",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, overlay, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        overlay = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.difference(input, overlay, output))
        return


class Erase(object):
    def __init__(self):
        self.label = "Erase"
        self.description = "Removes all the features, or parts of features, that overlap with the features of the erase vector polygon."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        erase = arcpy.Parameter(
            displayName="Input Erase Polygon Vector File",
            name="erase",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        erase.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, erase, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        erase = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.erase(input, erase, output))
        return


class ErasePolygonFromRaster(object):
    def __init__(self):
        self.label = "Erase Polygon From Raster"
        self.description = "Erases (cuts out) a vector polygon from a raster."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        polygons = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="polygons",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        polygons.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, polygons, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        polygons = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.erase_polygon_from_raster(input, polygons, output))
        return


class HighestPosition(object):
    def __init__(self):
        self.label = "Highest Position"
        self.description = "Identifies the stack position of the maximum value within a raster stack on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.highest_position(inputs, output))
        return


class Intersect(object):
    def __init__(self):
        self.label = "Intersect"
        self.description = "Identifies the parts of features in common between two input vector layers."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        overlay = arcpy.Parameter(
            displayName="Input Overlay Vector File",
            name="overlay",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        snap = arcpy.Parameter(
            displayName="Snap Tolerance",
            name="snap",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        snap.value = "0.0"

        params = [input, overlay, output, snap]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        overlay = parameters[1].valueAsText
        output = parameters[2].valueAsText
        snap = parameters[3].valueAsText
        messages.addMessage(wbt.intersect(input, overlay, output, snap))
        return


class LineIntersections(object):
    def __init__(self):
        self.label = "Line Intersections"
        self.description = "Identifies points where the features of two vector line layers intersect."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="input1",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input1.filter.list = ["Polyline", "Polygon"]

        input2 = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="input2",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input2.filter.list = ["Polyline", "Polygon"]

        output = arcpy.Parameter(
            displayName="Output Vector Point File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.line_intersections(input1, input2, output))
        return


class LowestPosition(object):
    def __init__(self):
        self.label = "Lowest Position"
        self.description = "Identifies the stack position of the minimum value within a raster stack on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.lowest_position(inputs, output))
        return


class MaxAbsoluteOverlay(object):
    def __init__(self):
        self.label = "Max Absolute Overlay"
        self.description = "Evaluates the maximum absolute value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.max_absolute_overlay(inputs, output))
        return


class MaxOverlay(object):
    def __init__(self):
        self.label = "Max Overlay"
        self.description = "Evaluates the maximum value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.max_overlay(inputs, output))
        return


class MinAbsoluteOverlay(object):
    def __init__(self):
        self.label = "Min Absolute Overlay"
        self.description = "Evaluates the minimum absolute value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.min_absolute_overlay(inputs, output))
        return


class MinOverlay(object):
    def __init__(self):
        self.label = "Min Overlay"
        self.description = "Evaluates the minimum value for each grid cell from a stack of input rasters."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.min_overlay(inputs, output))
        return


class PercentEqualTo(object):
    def __init__(self):
        self.label = "Percent Equal To"
        self.description = "Calculates the percentage of a raster stack that have cell values equal to an input on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        comparison = arcpy.Parameter(
            displayName="Input Comparison File",
            name="comparison",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, comparison, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        comparison = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.percent_equal_to(inputs, comparison, output))
        return


class PercentGreaterThan(object):
    def __init__(self):
        self.label = "Percent Greater Than"
        self.description = "Calculates the percentage of a raster stack that have cell values greather than an input on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        comparison = arcpy.Parameter(
            displayName="Input Comparison File",
            name="comparison",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, comparison, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        comparison = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.percent_greater_than(inputs, comparison, output))
        return


class PercentLessThan(object):
    def __init__(self):
        self.label = "Percent Less Than"
        self.description = "Calculates the percentage of a raster stack that have cell values less than an input on a cell-by-cell basis."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        comparison = arcpy.Parameter(
            displayName="Input Comparison File",
            name="comparison",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, comparison, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        comparison = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.percent_less_than(inputs, comparison, output))
        return


class PickFromList(object):
    def __init__(self):
        self.label = "Pick From List"
        self.description = "Outputs the value from a raster stack specified by a position raster."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        pos_input = arcpy.Parameter(
            displayName="Input Position File",
            name="pos_input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, pos_input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        pos_input = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.pick_from_list(inputs, pos_input, output))
        return


class Polygonize(object):
    def __init__(self):
        self.label = "Polygonize"
        self.description = "Creates a polygon layer from two or more intersecting line features contained in one or more input vector line files."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="inputs",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        inputs.filter.list = ["Polyline"]

        output = arcpy.Parameter(
            displayName="Output Vector Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.polygonize(inputs, output))
        return


class SplitWithLines(object):
    def __init__(self):
        self.label = "Split With Lines"
        self.description = "Splits the lines or polygons in one layer using the lines in another layer."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Lines or Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        split = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="split",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        split.filter.list = ["Polyline"]

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        params = [input, split, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        split = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.split_with_lines(input, split, output))
        return


class SumOverlay(object):
    def __init__(self):
        self.label = "Sum Overlay"
        self.description = "Calculates the sum for each grid cell from a group of raster images."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.sum_overlay(inputs, output))
        return


class SymmetricalDifference(object):
    def __init__(self):
        self.label = "Symmetrical Difference"
        self.description = "Outputs the features that occur in one of the two vector inputs but not both, i.e. no overlapping features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        overlay = arcpy.Parameter(
            displayName="Input Overlay Vector File",
            name="overlay",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        snap = arcpy.Parameter(
            displayName="Snap Tolerance",
            name="snap",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        snap.value = "0.0"

        params = [input, overlay, output, snap]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        overlay = parameters[1].valueAsText
        output = parameters[2].valueAsText
        snap = parameters[3].valueAsText
        messages.addMessage(wbt.symmetrical_difference(input, overlay, output, snap))
        return


class Union(object):
    def __init__(self):
        self.label = "Union"
        self.description = "Splits vector layers at their overlaps, creating a layer containing all the portions from both input and overlay layers."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        overlay = arcpy.Parameter(
            displayName="Input Overlay Vector File",
            name="overlay",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Vector File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")

        snap = arcpy.Parameter(
            displayName="Snap Tolerance",
            name="snap",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        snap.value = "0.0"

        params = [input, overlay, output, snap]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        overlay = parameters[1].valueAsText
        output = parameters[2].valueAsText
        snap = parameters[3].valueAsText
        messages.addMessage(wbt.union(input, overlay, output, snap))
        return


class WeightedOverlay(object):
    def __init__(self):
        self.label = "Weighted Overlay"
        self.description = "Performs a weighted sum on multiple input rasters after converting each image to a common scale. The tool performs a multi-criteria evaluation (MCE)."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        factors = arcpy.Parameter(
            displayName="Input Factor Files",
            name="factors",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        weights = arcpy.Parameter(
            displayName="Weight Values (e.g. 1.7;3.5;1.2)",
            name="weights",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        cost = arcpy.Parameter(
            displayName="Cost Factor? (e.g. false;true;true)",
            name="cost",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        constraints = arcpy.Parameter(
            displayName="Input Constraints Files",
            name="constraints",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        scale_max = arcpy.Parameter(
            displayName="Suitability Scale Maximum",
            name="scale_max",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        scale_max.value = "1.0"

        params = [factors, weights, cost, constraints, output, scale_max]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        factors = parameters[0].valueAsText
        weights = parameters[1].valueAsText
        cost = parameters[2].valueAsText
        constraints = parameters[3].valueAsText
        output = parameters[4].valueAsText
        scale_max = parameters[5].valueAsText
        messages.addMessage(wbt.weighted_overlay(factors, weights, cost, constraints, output, scale_max))
        return


class WeightedSum(object):
    def __init__(self):
        self.label = "Weighted Sum"
        self.description = "Performs a weighted-sum overlay on multiple input raster images."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        weights = arcpy.Parameter(
            displayName="Weight Values (e.g. 1.7;3.5;1.2)",
            name="weights",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [inputs, weights, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        weights = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.weighted_sum(inputs, weights, output))
        return


class CompactnessRatio(object):
    def __init__(self):
        self.label = "Compactness Ratio"
        self.description = "Calculates the compactness ratio (A/P), a measure of shape complexity, for vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.compactness_ratio(input))
        return


class EdgeProportion(object):
    def __init__(self):
        self.label = "Edge Proportion"
        self.description = "Calculate the proportion of cells in a raster polygon that are edge cells."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        output_text = arcpy.Parameter(
            displayName="Output a text report?",
            name="output_text",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input, output, output_text]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        output_text = parameters[2].valueAsText
        messages.addMessage(wbt.edge_proportion(input, output, output_text))
        return


class ElongationRatio(object):
    def __init__(self):
        self.label = "Elongation Ratio"
        self.description = "Calculates the elongation ratio for vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.elongation_ratio(input))
        return


class FindPatchOrClassEdgeCells(object):
    def __init__(self):
        self.label = "Find Patch Or Class Edge Cells"
        self.description = "Finds all cells located on the edge of patch or class features."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.find_patch_or_class_edge_cells(input, output))
        return


class HoleProportion(object):
    def __init__(self):
        self.label = "Hole Proportion"
        self.description = "Calculates the proportion of the total area of a polygon's holes relative to the area of the polygon's hull."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.hole_proportion(input))
        return


class LinearityIndex(object):
    def __init__(self):
        self.label = "Linearity Index"
        self.description = "Calculates the linearity index for vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.linearity_index(input))
        return


class PatchOrientation(object):
    def __init__(self):
        self.label = "Patch Orientation"
        self.description = "Calculates the orientation of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.patch_orientation(input))
        return


class PerimeterAreaRatio(object):
    def __init__(self):
        self.label = "Perimeter Area Ratio"
        self.description = "Calculates the perimeter-area ratio of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.perimeter_area_ratio(input))
        return


class RadiusOfGyration(object):
    def __init__(self):
        self.label = "Radius Of Gyration"
        self.description = "Calculates the distance of cells from their polygon's centroid."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        text_output = arcpy.Parameter(
            displayName="Output text?",
            name="text_output",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")

        params = [input, output, text_output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        text_output = parameters[2].valueAsText
        messages.addMessage(wbt.radius_of_gyration(input, output, text_output))
        return


class RelatedCircumscribingCircle(object):
    def __init__(self):
        self.label = "Related Circumscribing Circle"
        self.description = "Calculates the related circumscribing circle of vector polygons."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.related_circumscribing_circle(input))
        return


class ShapeComplexityIndex(object):
    def __init__(self):
        self.label = "Shape Complexity Index"
        self.description = "Calculates overall polygon shape complexity or irregularity."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Polygon"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.shape_complexity_index(input))
        return


class Aspect(object):
    def __init__(self):
        self.label = "Aspect"
        self.description = "Calculates an aspect raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.aspect(dem, output, zfactor))
        return


class CircularVarianceOfAspect(object):
    def __init__(self):
        self.label = "Circular Variance Of Aspect"
        self.description = "Calculates the circular variance of aspect at a scale for a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Roughness Scale File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Dimension",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "11"

        params = [dem, output, filter]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        messages.addMessage(wbt.circular_variance_of_aspect(dem, output, filter))
        return


class DevFromMeanElev(object):
    def __init__(self):
        self.label = "Dev From Mean Elev"
        self.description = "Calculates deviation from mean elevation."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [dem, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.dev_from_mean_elev(dem, output, filterx, filtery))
        return


class DiffFromMeanElev(object):
    def __init__(self):
        self.label = "Diff From Mean Elev"
        self.description = "Calculates difference from mean elevation (equivalent to a high-pass filter)."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [dem, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.diff_from_mean_elev(dem, output, filterx, filtery))
        return


class DirectionalRelief(object):
    def __init__(self):
        self.label = "Directional Relief"
        self.description = "Calculates relief for cells in an input DEM for a specified direction."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        azimuth = arcpy.Parameter(
            displayName="Azimuth",
            name="azimuth",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        azimuth.value = "0.0"

        max_dist = arcpy.Parameter(
            displayName="Maximum Search Distance",
            name="max_dist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, azimuth, max_dist]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        azimuth = parameters[2].valueAsText
        max_dist = parameters[3].valueAsText
        messages.addMessage(wbt.directional_relief(dem, output, azimuth, max_dist))
        return


class DownslopeIndex(object):
    def __init__(self):
        self.label = "Downslope Index"
        self.description = "Calculates the Hjerdt et al. (2004) downslope index."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        drop = arcpy.Parameter(
            displayName="Verical Drop",
            name="drop",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        drop.value = "2.0"

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['tangent', 'degrees', 'radians', 'distance']

        out_type.value = "tangent"

        params = [dem, output, drop, out_type]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        drop = parameters[2].valueAsText
        out_type = parameters[3].valueAsText
        messages.addMessage(wbt.downslope_index(dem, output, drop, out_type))
        return


class DrainagePreservingSmoothing(object):
    def __init__(self):
        self.label = "Drainage Preserving Smoothing"
        self.description = "Reduces short-scale variation in an input DEM while preserving breaks-in-slope and small drainage features using a modified Sun et al. (2007) algorithm."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Size",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "11"

        norm_diff = arcpy.Parameter(
            displayName="Normal Difference Threshold",
            name="norm_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        norm_diff.value = "15.0"

        num_iter = arcpy.Parameter(
            displayName="Iterations",
            name="num_iter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        num_iter.value = "3"

        max_diff = arcpy.Parameter(
            displayName="Maximum Elevation Change",
            name="max_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        reduction = arcpy.Parameter(
            displayName="Max. Smoothing Reduction Factor (%)",
            name="reduction",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        reduction.value = "80.0"

        dfm = arcpy.Parameter(
            displayName="Diff. From Median Threshold",
            name="dfm",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        dfm.value = "0.15"

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, filter, norm_diff, num_iter, max_diff, reduction, dfm, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        norm_diff = parameters[3].valueAsText
        num_iter = parameters[4].valueAsText
        max_diff = parameters[5].valueAsText
        reduction = parameters[6].valueAsText
        dfm = parameters[7].valueAsText
        zfactor = parameters[8].valueAsText
        messages.addMessage(wbt.drainage_preserving_smoothing(dem, output, filter, norm_diff, num_iter, max_diff, reduction, dfm, zfactor))
        return


class EdgeDensity(object):
    def __init__(self):
        self.label = "Edge Density"
        self.description = "Calculates the density of edges, or breaks-in-slope within DEMs."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Size",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "11"

        norm_diff = arcpy.Parameter(
            displayName="Normal Difference Threshold",
            name="norm_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        norm_diff.value = "5.0"

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, filter, norm_diff, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        norm_diff = parameters[3].valueAsText
        zfactor = parameters[4].valueAsText
        messages.addMessage(wbt.edge_density(dem, output, filter, norm_diff, zfactor))
        return


class ElevAbovePit(object):
    def __init__(self):
        self.label = "Elev Above Pit"
        self.description = "Calculate the elevation of each grid cell above the nearest downstream pit cell or grid edge cell."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.elev_above_pit(dem, output))
        return


class ElevPercentile(object):
    def __init__(self):
        self.label = "Elev Percentile"
        self.description = "Calculates the elevation percentile raster from a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        sig_digits = arcpy.Parameter(
            displayName="Number of Significant Digits",
            name="sig_digits",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        sig_digits.value = "2"

        params = [dem, output, filterx, filtery, sig_digits]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        sig_digits = parameters[4].valueAsText
        messages.addMessage(wbt.elev_percentile(dem, output, filterx, filtery, sig_digits))
        return


class ElevRelativeToMinMax(object):
    def __init__(self):
        self.label = "Elev Relative To Min Max"
        self.description = "Calculates the elevation of a location relative to the minimum and maximum elevations in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.elev_relative_to_min_max(dem, output))
        return


class ElevRelativeToWatershedMinMax(object):
    def __init__(self):
        self.label = "Elev Relative To Watershed Min Max"
        self.description = "Calculates the elevation of a location relative to the minimum and maximum elevations in a watershed."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        watersheds = arcpy.Parameter(
            displayName="Input Watersheds File",
            name="watersheds",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, watersheds, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        watersheds = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.elev_relative_to_watershed_min_max(dem, watersheds, output))
        return


class FeaturePreservingDenoise(object):
    def __init__(self):
        self.label = "Feature Preserving Denoise"
        self.description = "Reduces short-scale variation in an input DEM using a modified Sun et al. (2007) algorithm."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Size",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "11"

        norm_diff = arcpy.Parameter(
            displayName="Normal Difference Threshold",
            name="norm_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        norm_diff.value = "15.0"

        num_iter = arcpy.Parameter(
            displayName="Iterations",
            name="num_iter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        num_iter.value = "3"

        max_diff = arcpy.Parameter(
            displayName="Maximum Elevation Change",
            name="max_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, filter, norm_diff, num_iter, max_diff, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        norm_diff = parameters[3].valueAsText
        num_iter = parameters[4].valueAsText
        max_diff = parameters[5].valueAsText
        zfactor = parameters[6].valueAsText
        messages.addMessage(wbt.feature_preserving_denoise(dem, output, filter, norm_diff, num_iter, max_diff, zfactor))
        return


class FetchAnalysis(object):
    def __init__(self):
        self.label = "Fetch Analysis"
        self.description = "Performs an analysis of fetch or upwind distance to an obstacle."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        azimuth = arcpy.Parameter(
            displayName="Azimuth (degrees)",
            name="azimuth",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        azimuth.value = "0.0"

        hgt_inc = arcpy.Parameter(
            displayName="Height Increment Value",
            name="hgt_inc",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        hgt_inc.value = "0.05"

        params = [dem, output, azimuth, hgt_inc]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        azimuth = parameters[2].valueAsText
        hgt_inc = parameters[3].valueAsText
        messages.addMessage(wbt.fetch_analysis(dem, output, azimuth, hgt_inc))
        return


class FillMissingData(object):
    def __init__(self):
        self.label = "Fill Missing Data"
        self.description = "Fills NoData holes in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Dimension",
            name="filter",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        filter.value = "11"

        weight = arcpy.Parameter(
            displayName="IDW Weight (Exponent) Value",
            name="weight",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        weight.value = "2.0"

        params = [input, output, filter, weight]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        weight = parameters[3].valueAsText
        messages.addMessage(wbt.fill_missing_data(input, output, filter, weight))
        return


class FindRidges(object):
    def __init__(self):
        self.label = "Find Ridges"
        self.description = "Identifies potential ridge and peak grid cells."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        line_thin = arcpy.Parameter(
            displayName="Perform line-thinning?",
            name="line_thin",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        line_thin.value = "true"

        params = [dem, output, line_thin]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        line_thin = parameters[2].valueAsText
        messages.addMessage(wbt.find_ridges(dem, output, line_thin))
        return


class Hillshade(object):
    def __init__(self):
        self.label = "Hillshade"
        self.description = "Calculates a hillshade raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        azimuth = arcpy.Parameter(
            displayName="Azimuth (degrees)",
            name="azimuth",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        azimuth.value = "315.0"

        altitude = arcpy.Parameter(
            displayName="Altitude (degrees)",
            name="altitude",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        altitude.value = "30.0"

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, azimuth, altitude, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        azimuth = parameters[2].valueAsText
        altitude = parameters[3].valueAsText
        zfactor = parameters[4].valueAsText
        messages.addMessage(wbt.hillshade(dem, output, azimuth, altitude, zfactor))
        return


class HorizonAngle(object):
    def __init__(self):
        self.label = "Horizon Angle"
        self.description = "Calculates horizon angle (maximum upwind slope) for each grid cell in an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        azimuth = arcpy.Parameter(
            displayName="Azimuth",
            name="azimuth",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        azimuth.value = "0.0"

        max_dist = arcpy.Parameter(
            displayName="Maximum Search Distance",
            name="max_dist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, azimuth, max_dist]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        azimuth = parameters[2].valueAsText
        max_dist = parameters[3].valueAsText
        messages.addMessage(wbt.horizon_angle(dem, output, azimuth, max_dist))
        return


class HypsometricAnalysis(object):
    def __init__(self):
        self.label = "Hypsometric Analysis"
        self.description = "Calculates a hypsometric curve for one or more DEMs."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input DEM Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        watershed = arcpy.Parameter(
            displayName="Input Watershed Files (optional)",
            name="watershed",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [inputs, watershed, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        watershed = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.hypsometric_analysis(inputs, watershed, output))
        return


class MaxAnisotropyDev(object):
    def __init__(self):
        self.label = "Max Anisotropy Dev"
        self.description = "Calculates the maximum anisotropy (directionality) in elevation deviation over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        out_mag = arcpy.Parameter(
            displayName="Output DEVmax Magnitude File",
            name="out_mag",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_mag.filter.list = ["tif"]

        out_scale = arcpy.Parameter(
            displayName="Output DEVmax Scale File",
            name="out_scale",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_scale.filter.list = ["tif"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        min_scale.value = "3"

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step.value = "2"

        params = [dem, out_mag, out_scale, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        out_mag = parameters[1].valueAsText
        out_scale = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.max_anisotropy_dev(dem, out_mag, out_scale, min_scale, max_scale, step))
        return


class MaxAnisotropyDevSignature(object):
    def __init__(self):
        self.label = "Max Anisotropy Dev Signature"
        self.description = "Calculates the anisotropy in deviation from mean for points over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        points = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="points",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        points.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        min_scale.value = "1"

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        step.value = "1"

        params = [dem, points, output, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        points = parameters[1].valueAsText
        output = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.max_anisotropy_dev_signature(dem, points, output, min_scale, max_scale, step))
        return


class MaxBranchLength(object):
    def __init__(self):
        self.label = "Max Branch Length"
        self.description = "Lindsay and Seibert's (2013) branch length index is used to map drainage divides or ridge lines."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        log = arcpy.Parameter(
            displayName="Log-transform the output?",
            name="log",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, log]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        log = parameters[2].valueAsText
        messages.addMessage(wbt.max_branch_length(dem, output, log))
        return


class MaxDifferenceFromMean(object):
    def __init__(self):
        self.label = "Max Difference From Mean"
        self.description = "Calculates the maximum difference from mean elevation over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        out_mag = arcpy.Parameter(
            displayName="Output DIFFmax Magnitude File",
            name="out_mag",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_mag.filter.list = ["tif"]

        out_scale = arcpy.Parameter(
            displayName="Output DIFFmax Scale File",
            name="out_scale",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_scale.filter.list = ["tif"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step.value = "1"

        params = [dem, out_mag, out_scale, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        out_mag = parameters[1].valueAsText
        out_scale = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.max_difference_from_mean(dem, out_mag, out_scale, min_scale, max_scale, step))
        return


class MaxDownslopeElevChange(object):
    def __init__(self):
        self.label = "Max Downslope Elev Change"
        self.description = "Calculates the maximum downslope change in elevation between a grid cell and its eight downslope neighbors."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.max_downslope_elev_change(dem, output))
        return


class MaxElevDevSignature(object):
    def __init__(self):
        self.label = "Max Elev Dev Signature"
        self.description = "Calculates the maximum elevation deviation over a range of spatial scales and for a set of points."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        points = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="points",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        points.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step.value = "10"

        params = [dem, points, output, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        points = parameters[1].valueAsText
        output = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.max_elev_dev_signature(dem, points, output, min_scale, max_scale, step))
        return


class MaxElevationDeviation(object):
    def __init__(self):
        self.label = "Max Elevation Deviation"
        self.description = "Calculates the maximum elevation deviation over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        out_mag = arcpy.Parameter(
            displayName="Output DEVmax Magnitude File",
            name="out_mag",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_mag.filter.list = ["tif"]

        out_scale = arcpy.Parameter(
            displayName="Output DEVmax Scale File",
            name="out_scale",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_scale.filter.list = ["tif"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step.value = "1"

        params = [dem, out_mag, out_scale, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        out_mag = parameters[1].valueAsText
        out_scale = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.max_elevation_deviation(dem, out_mag, out_scale, min_scale, max_scale, step))
        return


class MinDownslopeElevChange(object):
    def __init__(self):
        self.label = "Min Downslope Elev Change"
        self.description = "Calculates the minimum downslope change in elevation between a grid cell and its eight downslope neighbors."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.min_downslope_elev_change(dem, output))
        return


class MultiscaleRoughness(object):
    def __init__(self):
        self.label = "Multiscale Roughness"
        self.description = "Calculates surface roughness over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        out_mag = arcpy.Parameter(
            displayName="Output Roughness Magnitude File",
            name="out_mag",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_mag.filter.list = ["tif"]

        out_scale = arcpy.Parameter(
            displayName="Output Roughness Scale File",
            name="out_scale",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_scale.filter.list = ["tif"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        min_scale.value = "1"

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        step.value = "1"

        params = [dem, out_mag, out_scale, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        out_mag = parameters[1].valueAsText
        out_scale = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.multiscale_roughness(dem, out_mag, out_scale, min_scale, max_scale, step))
        return


class MultiscaleRoughnessSignature(object):
    def __init__(self):
        self.label = "Multiscale Roughness Signature"
        self.description = "Calculates the surface roughness for points over a range of spatial scales."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        points = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="points",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        points.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        min_scale = arcpy.Parameter(
            displayName="Minimum Search Neighbourhood Radius (grid cells)",
            name="min_scale",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        min_scale.value = "1"

        max_scale = arcpy.Parameter(
            displayName="Maximum Search Neighbourhood Radius (grid cells)",
            name="max_scale",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        step = arcpy.Parameter(
            displayName="Step Size",
            name="step",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        step.value = "1"

        params = [dem, points, output, min_scale, max_scale, step]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        points = parameters[1].valueAsText
        output = parameters[2].valueAsText
        min_scale = parameters[3].valueAsText
        max_scale = parameters[4].valueAsText
        step = parameters[5].valueAsText
        messages.addMessage(wbt.multiscale_roughness_signature(dem, points, output, min_scale, max_scale, step))
        return


class MultiscaleTopographicPositionImage(object):
    def __init__(self):
        self.label = "Multiscale Topographic Position Image"
        self.description = "Creates a multiscale topographic position image from three DEVmax rasters of differing spatial scale ranges."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        local = arcpy.Parameter(
            displayName="Input Local-Scale File",
            name="local",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        meso = arcpy.Parameter(
            displayName="Input Meso-Scale File",
            name="meso",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        broad = arcpy.Parameter(
            displayName="Input Broad-Scale File",
            name="broad",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        lightness = arcpy.Parameter(
            displayName="Image Lightness Value",
            name="lightness",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        lightness.value = "1.2"

        params = [local, meso, broad, output, lightness]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        local = parameters[0].valueAsText
        meso = parameters[1].valueAsText
        broad = parameters[2].valueAsText
        output = parameters[3].valueAsText
        lightness = parameters[4].valueAsText
        messages.addMessage(wbt.multiscale_topographic_position_image(local, meso, broad, output, lightness))
        return


class NumDownslopeNeighbours(object):
    def __init__(self):
        self.label = "Num Downslope Neighbours"
        self.description = "Calculates the number of downslope neighbours to each grid cell in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.num_downslope_neighbours(dem, output))
        return


class NumUpslopeNeighbours(object):
    def __init__(self):
        self.label = "Num Upslope Neighbours"
        self.description = "Calculates the number of upslope neighbours to each grid cell in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.num_upslope_neighbours(dem, output))
        return


class PennockLandformClass(object):
    def __init__(self):
        self.label = "Pennock Landform Class"
        self.description = "Classifies hillslope zones based on slope, profile curvature, and plan curvature."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        slope = arcpy.Parameter(
            displayName="Slope Threshold (degrees)",
            name="slope",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        slope.value = "3.0"

        prof = arcpy.Parameter(
            displayName="Profile Curvature Threshold",
            name="prof",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        prof.value = "0.1"

        plan = arcpy.Parameter(
            displayName="Plan Curvature Threshold",
            name="plan",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        plan.value = "0.0"

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, slope, prof, plan, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        slope = parameters[2].valueAsText
        prof = parameters[3].valueAsText
        plan = parameters[4].valueAsText
        zfactor = parameters[5].valueAsText
        messages.addMessage(wbt.pennock_landform_class(dem, output, slope, prof, plan, zfactor))
        return


class PercentElevRange(object):
    def __init__(self):
        self.label = "Percent Elev Range"
        self.description = "Calculates percent of elevation range from a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "3"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "3"

        params = [dem, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.percent_elev_range(dem, output, filterx, filtery))
        return


class PlanCurvature(object):
    def __init__(self):
        self.label = "Plan Curvature"
        self.description = "Calculates a plan (contour) curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.plan_curvature(dem, output, zfactor))
        return


class Profile(object):
    def __init__(self):
        self.label = "Profile"
        self.description = "Plots profiles from digital surface models."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        lines = arcpy.Parameter(
            displayName="Input Vector Line File",
            name="lines",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        lines.filter.list = ["Polyline"]

        surface = arcpy.Parameter(
            displayName="Input Surface File",
            name="surface",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [lines, surface, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        lines = parameters[0].valueAsText
        surface = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.profile(lines, surface, output))
        return


class ProfileCurvature(object):
    def __init__(self):
        self.label = "Profile Curvature"
        self.description = "Calculates a profile curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.profile_curvature(dem, output, zfactor))
        return


class RelativeAspect(object):
    def __init__(self):
        self.label = "Relative Aspect"
        self.description = "Calculates relative aspect (relative to a user-specified direction) from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        azimuth = arcpy.Parameter(
            displayName="Azimuth",
            name="azimuth",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        azimuth.value = "0.0"

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, azimuth, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        azimuth = parameters[2].valueAsText
        zfactor = parameters[3].valueAsText
        messages.addMessage(wbt.relative_aspect(dem, output, azimuth, zfactor))
        return


class RelativeStreamPowerIndex(object):
    def __init__(self):
        self.label = "Relative Stream Power Index"
        self.description = "Calculates the relative stream power index."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        sca = arcpy.Parameter(
            displayName="Input Specific Contributing Area (SCA) File",
            name="sca",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        slope = arcpy.Parameter(
            displayName="Input Slope File",
            name="slope",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        exponent = arcpy.Parameter(
            displayName="Specific Contributing Area (SCA) Exponent",
            name="exponent",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        exponent.value = "1.0"

        params = [sca, slope, output, exponent]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        sca = parameters[0].valueAsText
        slope = parameters[1].valueAsText
        output = parameters[2].valueAsText
        exponent = parameters[3].valueAsText
        messages.addMessage(wbt.relative_stream_power_index(sca, slope, output, exponent))
        return


class RelativeTopographicPosition(object):
    def __init__(self):
        self.label = "Relative Topographic Position"
        self.description = "Calculates the relative topographic position index from a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [dem, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.relative_topographic_position(dem, output, filterx, filtery))
        return


class RemoveOffTerrainObjects(object):
    def __init__(self):
        self.label = "Remove Off Terrain Objects"
        self.description = "Removes off-terrain objects from a raster digital elevation model (DEM)."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Dimension",
            name="filter",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        filter.value = "11"

        slope = arcpy.Parameter(
            displayName="Slope Threshold",
            name="slope",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        slope.value = "15.0"

        params = [dem, output, filter, slope]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        slope = parameters[3].valueAsText
        messages.addMessage(wbt.remove_off_terrain_objects(dem, output, filter, slope))
        return


class RuggednessIndex(object):
    def __init__(self):
        self.label = "Ruggedness Index"
        self.description = "Calculates the Riley et al.'s (1999) terrain ruggedness index from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.ruggedness_index(dem, output, zfactor))
        return


class SedimentTransportIndex(object):
    def __init__(self):
        self.label = "Sediment Transport Index"
        self.description = "Calculates the sediment transport index."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        sca = arcpy.Parameter(
            displayName="Input Specific Contributing Area (SCA) File",
            name="sca",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        slope = arcpy.Parameter(
            displayName="Input Slope File",
            name="slope",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sca_exponent = arcpy.Parameter(
            displayName="Specific Contributing Area (SCA) Exponent",
            name="sca_exponent",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        sca_exponent.value = "0.4"

        slope_exponent = arcpy.Parameter(
            displayName="Slope Exponent",
            name="slope_exponent",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        slope_exponent.value = "1.3"

        params = [sca, slope, output, sca_exponent, slope_exponent]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        sca = parameters[0].valueAsText
        slope = parameters[1].valueAsText
        output = parameters[2].valueAsText
        sca_exponent = parameters[3].valueAsText
        slope_exponent = parameters[4].valueAsText
        messages.addMessage(wbt.sediment_transport_index(sca, slope, output, sca_exponent, slope_exponent))
        return


class Slope(object):
    def __init__(self):
        self.label = "Slope"
        self.description = "Calculates a slope raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.slope(dem, output, zfactor))
        return


class SlopeVsElevationPlot(object):
    def __init__(self):
        self.label = "Slope Vs Elevation Plot"
        self.description = "Creates a slope vs. elevation plot for one or more DEMs."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input DEM Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        watershed = arcpy.Parameter(
            displayName="Input Watershed Files (optional)",
            name="watershed",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [inputs, watershed, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        watershed = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.slope_vs_elevation_plot(inputs, watershed, output))
        return


class StandardDeviationOfSlope(object):
    def __init__(self):
        self.label = "Standard Deviation Of Slope"
        self.description = "Calculates the standard deviation of slope from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, zfactor, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        filterx = parameters[3].valueAsText
        filtery = parameters[4].valueAsText
        messages.addMessage(wbt.standard_deviation_of_slope(input, output, zfactor, filterx, filtery))
        return


class SurfaceAreaRatio(object):
    def __init__(self):
        self.label = "Surface Area Ratio"
        self.description = "Calculates a the surface area ratio of each grid cell in an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.surface_area_ratio(dem, output))
        return


class TangentialCurvature(object):
    def __init__(self):
        self.label = "Tangential Curvature"
        self.description = "Calculates a tangential curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.tangential_curvature(dem, output, zfactor))
        return


class TotalCurvature(object):
    def __init__(self):
        self.label = "Total Curvature"
        self.description = "Calculates a total curvature raster from an input DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zfactor = arcpy.Parameter(
            displayName="Z Conversion Factor",
            name="zfactor",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        zfactor.value = "1.0"

        params = [dem, output, zfactor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zfactor = parameters[2].valueAsText
        messages.addMessage(wbt.total_curvature(dem, output, zfactor))
        return


class Viewshed(object):
    def __init__(self):
        self.label = "Viewshed"
        self.description = "Identifies the viewshed for a point or set of points."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        stations = arcpy.Parameter(
            displayName="Viewing Station Vector File",
            name="stations",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        stations.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        height = arcpy.Parameter(
            displayName="Station Height (in z units)",
            name="height",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        height.value = "2.0"

        params = [dem, stations, output, height]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        stations = parameters[1].valueAsText
        output = parameters[2].valueAsText
        height = parameters[3].valueAsText
        messages.addMessage(wbt.viewshed(dem, stations, output, height))
        return


class VisibilityIndex(object):
    def __init__(self):
        self.label = "Visibility Index"
        self.description = "Estimates the relative visibility of sites in a DEM."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        height = arcpy.Parameter(
            displayName="Station Height (in z units)",
            name="height",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        height.value = "2.0"

        res_factor = arcpy.Parameter(
            displayName="Resolution Factor",
            name="res_factor",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        res_factor.value = "2"

        params = [dem, output, height, res_factor]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        height = parameters[2].valueAsText
        res_factor = parameters[3].valueAsText
        messages.addMessage(wbt.visibility_index(dem, output, height, res_factor))
        return


class WetnessIndex(object):
    def __init__(self):
        self.label = "Wetness Index"
        self.description = "Calculates the topographic wetness index, Ln(A / tan(slope))."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        sca = arcpy.Parameter(
            displayName="Input Specific Contributing Area (SCA) File",
            name="sca",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        slope = arcpy.Parameter(
            displayName="Input Slope File",
            name="slope",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [sca, slope, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        sca = parameters[0].valueAsText
        slope = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.wetness_index(sca, slope, output))
        return


class AverageFlowpathSlope(object):
    def __init__(self):
        self.label = "Average Flowpath Slope"
        self.description = "Measures the average slope gradient from each grid cell to all upslope divide cells."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.average_flowpath_slope(dem, output))
        return


class AverageUpslopeFlowpathLength(object):
    def __init__(self):
        self.label = "Average Upslope Flowpath Length"
        self.description = "Measures the average length of all upslope flowpaths draining each grid cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.average_upslope_flowpath_length(dem, output))
        return


class Basins(object):
    def __init__(self):
        self.label = "Basins"
        self.description = "Identifies drainage basins that drain to the DEM edge."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        output = parameters[1].valueAsText
        esri_pntr = parameters[2].valueAsText
        messages.addMessage(wbt.basins(d8_pntr, output, esri_pntr))
        return


class BreachDepressions(object):
    def __init__(self):
        self.label = "Breach Depressions"
        self.description = "Breaches all of the depressions in a DEM using Lindsay's (2016) algorithm. This should be preferred over depression filling in most cases."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        max_depth = arcpy.Parameter(
            displayName="Maximum Breach Depth (z units)",
            name="max_depth",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        max_length = arcpy.Parameter(
            displayName="Maximum Breach Channel Length (grid cells)",
            name="max_length",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, max_depth, max_length]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        max_depth = parameters[2].valueAsText
        max_length = parameters[3].valueAsText
        messages.addMessage(wbt.breach_depressions(dem, output, max_depth, max_length))
        return


class BreachSingleCellPits(object):
    def __init__(self):
        self.label = "Breach Single Cell Pits"
        self.description = "Removes single-cell pits from an input DEM by breaching."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.breach_single_cell_pits(dem, output))
        return


class D8FlowAccumulation(object):
    def __init__(self):
        self.label = "D8 Flow Accumulation"
        self.description = "Calculates a D8 flow accumulation raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['cells', 'catchment area', 'specific contributing area']

        out_type.value = "cells"

        log = arcpy.Parameter(
            displayName="Log-transform the output?",
            name="log",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        clip = arcpy.Parameter(
            displayName="Clip the upper tail by 1%?",
            name="clip",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, out_type, log, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_type = parameters[2].valueAsText
        log = parameters[3].valueAsText
        clip = parameters[4].valueAsText
        messages.addMessage(wbt.d8_flow_accumulation(dem, output, out_type, log, clip))
        return


class D8MassFlux(object):
    def __init__(self):
        self.label = "D8 Mass Flux"
        self.description = "Performs a D8 mass flux calculation."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        loading = arcpy.Parameter(
            displayName="Input Loading File",
            name="loading",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        efficiency = arcpy.Parameter(
            displayName="Input Efficiency File",
            name="efficiency",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        absorption = arcpy.Parameter(
            displayName="Input Absorption File",
            name="absorption",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, loading, efficiency, absorption, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        loading = parameters[1].valueAsText
        efficiency = parameters[2].valueAsText
        absorption = parameters[3].valueAsText
        output = parameters[4].valueAsText
        messages.addMessage(wbt.d8_mass_flux(dem, loading, efficiency, absorption, output))
        return


class D8Pointer(object):
    def __init__(self):
        self.label = "D8 Pointer"
        self.description = "Calculates a D8 flow pointer raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Should the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [dem, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        esri_pntr = parameters[2].valueAsText
        messages.addMessage(wbt.d8_pointer(dem, output, esri_pntr))
        return


class DInfFlowAccumulation(object):
    def __init__(self):
        self.label = "D Inf Flow Accumulation"
        self.description = "Calculates a D-infinity flow accumulation raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['Cells', 'Specific Contributing Area', 'Catchment Area']

        out_type.value = "Specific Contributing Area"

        threshold = arcpy.Parameter(
            displayName="Convergence Threshold (grid cells; blank for none)",
            name="threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        log = arcpy.Parameter(
            displayName="Log-transform the output?",
            name="log",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        clip = arcpy.Parameter(
            displayName="Clip the upper tail by 1%?",
            name="clip",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, out_type, threshold, log, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_type = parameters[2].valueAsText
        threshold = parameters[3].valueAsText
        log = parameters[4].valueAsText
        clip = parameters[5].valueAsText
        messages.addMessage(wbt.d_inf_flow_accumulation(dem, output, out_type, threshold, log, clip))
        return


class DInfMassFlux(object):
    def __init__(self):
        self.label = "D Inf Mass Flux"
        self.description = "Performs a D-infinity mass flux calculation."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        loading = arcpy.Parameter(
            displayName="Input Loading File",
            name="loading",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        efficiency = arcpy.Parameter(
            displayName="Input Efficiency File",
            name="efficiency",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        absorption = arcpy.Parameter(
            displayName="Input Absorption File",
            name="absorption",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, loading, efficiency, absorption, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        loading = parameters[1].valueAsText
        efficiency = parameters[2].valueAsText
        absorption = parameters[3].valueAsText
        output = parameters[4].valueAsText
        messages.addMessage(wbt.d_inf_mass_flux(dem, loading, efficiency, absorption, output))
        return


class DInfPointer(object):
    def __init__(self):
        self.label = "D Inf Pointer"
        self.description = "Calculates a D-infinity flow pointer (flow direction) raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.d_inf_pointer(dem, output))
        return


class DepthInSink(object):
    def __init__(self):
        self.label = "Depth In Sink"
        self.description = "Measures the depth of sinks (depressions) in a DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zero_background = parameters[2].valueAsText
        messages.addMessage(wbt.depth_in_sink(dem, output, zero_background))
        return


class DownslopeDistanceToStream(object):
    def __init__(self):
        self.label = "Downslope Distance To Stream"
        self.description = "Measures distance to the nearest downslope stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, streams, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.downslope_distance_to_stream(dem, streams, output))
        return


class DownslopeFlowpathLength(object):
    def __init__(self):
        self.label = "Downslope Flowpath Length"
        self.description = "Calculates the downslope flowpath length from each cell to basin outlet."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        watersheds = arcpy.Parameter(
            displayName="Input Watersheds File (optional)",
            name="watersheds",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        weights = arcpy.Parameter(
            displayName="Input Weights File (optional)",
            name="weights",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, watersheds, weights, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        watersheds = parameters[1].valueAsText
        weights = parameters[2].valueAsText
        output = parameters[3].valueAsText
        esri_pntr = parameters[4].valueAsText
        messages.addMessage(wbt.downslope_flowpath_length(d8_pntr, watersheds, weights, output, esri_pntr))
        return


class ElevationAboveStream(object):
    def __init__(self):
        self.label = "Elevation Above Stream"
        self.description = "Calculates the elevation of cells above the nearest downslope stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, streams, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.elevation_above_stream(dem, streams, output))
        return


class ElevationAboveStreamEuclidean(object):
    def __init__(self):
        self.label = "Elevation Above Stream Euclidean"
        self.description = "Calculates the elevation of cells above the nearest (Euclidean distance) stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, streams, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.elevation_above_stream_euclidean(dem, streams, output))
        return


class Fd8FlowAccumulation(object):
    def __init__(self):
        self.label = "Fd8 Flow Accumulation"
        self.description = "Calculates an FD8 flow accumulation raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['cells', 'specific contributing area', 'catchment area']

        out_type.value = "specific contributing area"

        exponent = arcpy.Parameter(
            displayName="Exponent Parameter",
            name="exponent",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        exponent.value = "1.1"

        threshold = arcpy.Parameter(
            displayName="Convergence Threshold (grid cells; blank for none)",
            name="threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        log = arcpy.Parameter(
            displayName="Log-transform the output?",
            name="log",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        clip = arcpy.Parameter(
            displayName="Clip the upper tail by 1%?",
            name="clip",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, out_type, exponent, threshold, log, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_type = parameters[2].valueAsText
        exponent = parameters[3].valueAsText
        threshold = parameters[4].valueAsText
        log = parameters[5].valueAsText
        clip = parameters[6].valueAsText
        messages.addMessage(wbt.fd8_flow_accumulation(dem, output, out_type, exponent, threshold, log, clip))
        return


class Fd8Pointer(object):
    def __init__(self):
        self.label = "Fd8 Pointer"
        self.description = "Calculates an FD8 flow pointer raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.fd8_pointer(dem, output))
        return


class FillBurn(object):
    def __init__(self):
        self.label = "Fill Burn"
        self.description = "Burns streams into a DEM using the FillBurn (Saunders, 1999) method."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Vector Streams File",
            name="streams",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        streams.filter.list = ["Polyline"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, streams, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.fill_burn(dem, streams, output))
        return


class FillDepressions(object):
    def __init__(self):
        self.label = "Fill Depressions"
        self.description = "Fills all of the depressions in a DEM. Depression breaching should be preferred in most cases."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        fix_flats = arcpy.Parameter(
            displayName="Fix flat areas?",
            name="fix_flats",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        fix_flats.value = "true"

        params = [dem, output, fix_flats]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        fix_flats = parameters[2].valueAsText
        messages.addMessage(wbt.fill_depressions(dem, output, fix_flats))
        return


class FillSingleCellPits(object):
    def __init__(self):
        self.label = "Fill Single Cell Pits"
        self.description = "Raises pit cells to the elevation of their lowest neighbour."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.fill_single_cell_pits(dem, output))
        return


class FindNoFlowCells(object):
    def __init__(self):
        self.label = "Find No Flow Cells"
        self.description = "Finds grid cells with no downslope neighbours."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.find_no_flow_cells(dem, output))
        return


class FindParallelFlow(object):
    def __init__(self):
        self.label = "Find Parallel Flow"
        self.description = "Finds areas of parallel flow in D8 flow direction rasters."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [d8_pntr, streams, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.find_parallel_flow(d8_pntr, streams, output))
        return


class FlattenLakes(object):
    def __init__(self):
        self.label = "Flatten Lakes"
        self.description = "Flattens lake polygons in a raster DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        lakes = arcpy.Parameter(
            displayName="Input Lakes Vector Polygon File",
            name="lakes",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        lakes.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, lakes, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        lakes = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.flatten_lakes(dem, lakes, output))
        return


class FloodOrder(object):
    def __init__(self):
        self.label = "Flood Order"
        self.description = "Assigns each DEM grid cell its order in the sequence of inundations that are encountered during a search starting from the edges, moving inward at increasing elevations."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.flood_order(dem, output))
        return


class FlowAccumulationFullWorkflow(object):
    def __init__(self):
        self.label = "Flow Accumulation Full Workflow"
        self.description = "Resolves all of the depressions in a DEM, outputting a breached DEM, an aspect-aligned non-divergent flow pointer, and a flow accumulation raster."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        out_dem = arcpy.Parameter(
            displayName="Output DEM File",
            name="out_dem",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_dem.filter.list = ["tif"]

        out_pntr = arcpy.Parameter(
            displayName="Output Flow Pointer File",
            name="out_pntr",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_pntr.filter.list = ["tif"]

        out_accum = arcpy.Parameter(
            displayName="Output Flow Accumulation File",
            name="out_accum",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_accum.filter.list = ["tif"]

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['Cells', 'Specific Contributing Area', 'Catchment Area']

        out_type.value = "Specific Contributing Area"

        log = arcpy.Parameter(
            displayName="Log-transform the output?",
            name="log",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        clip = arcpy.Parameter(
            displayName="Clip the upper tail by 1%?",
            name="clip",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [dem, out_dem, out_pntr, out_accum, out_type, log, clip, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        out_dem = parameters[1].valueAsText
        out_pntr = parameters[2].valueAsText
        out_accum = parameters[3].valueAsText
        out_type = parameters[4].valueAsText
        log = parameters[5].valueAsText
        clip = parameters[6].valueAsText
        esri_pntr = parameters[7].valueAsText
        messages.addMessage(wbt.flow_accumulation_full_workflow(dem, out_dem, out_pntr, out_accum, out_type, log, clip, esri_pntr))
        return


class FlowLengthDiff(object):
    def __init__(self):
        self.label = "Flow Length Diff"
        self.description = "Calculates the local maximum absolute difference in downslope flowpath length, useful in mapping drainage divides and ridges."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        output = parameters[1].valueAsText
        esri_pntr = parameters[2].valueAsText
        messages.addMessage(wbt.flow_length_diff(d8_pntr, output, esri_pntr))
        return


class Hillslopes(object):
    def __init__(self):
        self.label = "Hillslopes"
        self.description = "Identifies the individual hillslopes draining to each link in a stream network."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, streams, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        messages.addMessage(wbt.hillslopes(d8_pntr, streams, output, esri_pntr))
        return


class ImpoundmentSizeIndex(object):
    def __init__(self):
        self.label = "Impoundment Size Index"
        self.description = "Calculates the impoundment size resulting from damming a DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_type = arcpy.Parameter(
            displayName="Output Type",
            name="out_type",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        out_type.filter.type = "ValueList"
        out_type.filter.list = ['depth', 'volume', 'area']

        out_type.value = "depth"

        damlength = arcpy.Parameter(
            displayName="Max dam length (grid cells)",
            name="damlength",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [dem, output, out_type, damlength]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_type = parameters[2].valueAsText
        damlength = parameters[3].valueAsText
        messages.addMessage(wbt.impoundment_size_index(dem, output, out_type, damlength))
        return


class Isobasins(object):
    def __init__(self):
        self.label = "Isobasins"
        self.description = "Divides a landscape into nearly equal sized drainage basins (i.e. watersheds)."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        size = arcpy.Parameter(
            displayName="Target Basin Size (grid cells)",
            name="size",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        params = [dem, output, size]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        size = parameters[2].valueAsText
        messages.addMessage(wbt.isobasins(dem, output, size))
        return


class JensonSnapPourPoints(object):
    def __init__(self):
        self.label = "Jenson Snap Pour Points"
        self.description = "Moves outlet points used to specify points of interest in a watershedding operation to the nearest stream cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        pour_pts = arcpy.Parameter(
            displayName="Input Pour Points (Outlet) File",
            name="pour_pts",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        pour_pts.filter.list = ["Point"]

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        snap_dist = arcpy.Parameter(
            displayName="Maximum Snap Distance (map units)",
            name="snap_dist",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [pour_pts, streams, output, snap_dist]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        pour_pts = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        snap_dist = parameters[3].valueAsText
        messages.addMessage(wbt.jenson_snap_pour_points(pour_pts, streams, output, snap_dist))
        return


class LongestFlowpath(object):
    def __init__(self):
        self.label = "Longest Flowpath"
        self.description = "Delineates the longest flowpaths for a group of subbasins or watersheds."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        basins = arcpy.Parameter(
            displayName="Basins File",
            name="basins",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polyline"]

        params = [dem, basins, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        basins = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.longest_flowpath(dem, basins, output))
        return


class MaxUpslopeFlowpathLength(object):
    def __init__(self):
        self.label = "Max Upslope Flowpath Length"
        self.description = "Measures the maximum length of all upslope flowpaths draining each grid cell."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.max_upslope_flowpath_length(dem, output))
        return


class NumInflowingNeighbours(object):
    def __init__(self):
        self.label = "Num Inflowing Neighbours"
        self.description = "Computes the number of inflowing neighbours to each cell in an input DEM based on the D8 algorithm."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [dem, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.num_inflowing_neighbours(dem, output))
        return


class RaiseWalls(object):
    def __init__(self):
        self.label = "Raise Walls"
        self.description = "Raises walls in a DEM along a line or around a polygon, e.g. a watershed."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Line or Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        breach = arcpy.Parameter(
            displayName="Input Breach Lines (optional)",
            name="breach",
            datatype="DEShapefile",
            parameterType="Optional",
            direction="Input")
        breach.filter.list = ["Polyline"]

        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        height = arcpy.Parameter(
            displayName="Wall Height",
            name="height",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        height.value = "100.0"

        params = [input, breach, dem, output, height]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        breach = parameters[1].valueAsText
        dem = parameters[2].valueAsText
        output = parameters[3].valueAsText
        height = parameters[4].valueAsText
        messages.addMessage(wbt.raise_walls(input, breach, dem, output, height))
        return


class Rho8Pointer(object):
    def __init__(self):
        self.label = "Rho8 Pointer"
        self.description = "Calculates a stochastic Rho8 flow pointer raster from an input DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Should the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [dem, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        esri_pntr = parameters[2].valueAsText
        messages.addMessage(wbt.rho8_pointer(dem, output, esri_pntr))
        return


class Sink(object):
    def __init__(self):
        self.label = "Sink"
        self.description = "Identifies the depressions in a DEM, giving each feature a unique identifier."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [dem, output, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        zero_background = parameters[2].valueAsText
        messages.addMessage(wbt.sink(dem, output, zero_background))
        return


class SnapPourPoints(object):
    def __init__(self):
        self.label = "Snap Pour Points"
        self.description = "Moves outlet points used to specify points of interest in a watershedding operation to the cell with the highest flow accumulation in its neighbourhood."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        pour_pts = arcpy.Parameter(
            displayName="Input Pour Points (Outlet) File",
            name="pour_pts",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        pour_pts.filter.list = ["Point"]

        flow_accum = arcpy.Parameter(
            displayName="Input D8 Flow Accumulation File",
            name="flow_accum",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Point"]

        snap_dist = arcpy.Parameter(
            displayName="Maximum Snap Distance (map units)",
            name="snap_dist",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [pour_pts, flow_accum, output, snap_dist]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        pour_pts = parameters[0].valueAsText
        flow_accum = parameters[1].valueAsText
        output = parameters[2].valueAsText
        snap_dist = parameters[3].valueAsText
        messages.addMessage(wbt.snap_pour_points(pour_pts, flow_accum, output, snap_dist))
        return


class StochasticDepressionAnalysis(object):
    def __init__(self):
        self.label = "Stochastic Depression Analysis"
        self.description = "Preforms a stochastic analysis of depressions within a DEM."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        rmse = arcpy.Parameter(
            displayName="DEM root-mean-square-error (z units)",
            name="rmse",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        range = arcpy.Parameter(
            displayName="Range of Autocorrelation (map units)",
            name="range",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        iterations = arcpy.Parameter(
            displayName="Iterations",
            name="iterations",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        iterations.value = "1000"

        params = [dem, output, rmse, range, iterations]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        rmse = parameters[2].valueAsText
        range = parameters[3].valueAsText
        iterations = parameters[4].valueAsText
        messages.addMessage(wbt.stochastic_depression_analysis(dem, output, rmse, range, iterations))
        return


class StrahlerOrderBasins(object):
    def __init__(self):
        self.label = "Strahler Order Basins"
        self.description = "Identifies Strahler-order basins from an input stream network."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, streams, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        messages.addMessage(wbt.strahler_order_basins(d8_pntr, streams, output, esri_pntr))
        return


class Subbasins(object):
    def __init__(self):
        self.label = "Subbasins"
        self.description = "Identifies the catchments, or sub-basin, draining to each link in a stream network."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, streams, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        messages.addMessage(wbt.subbasins(d8_pntr, streams, output, esri_pntr))
        return


class TraceDownslopeFlowpaths(object):
    def __init__(self):
        self.label = "Trace Downslope Flowpaths"
        self.description = "Traces downslope flowpaths from one or more target sites (i.e. seed points)."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        seed_pts = arcpy.Parameter(
            displayName="Input Vector Seed Points File",
            name="seed_pts",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        seed_pts.filter.list = ["Point"]

        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [seed_pts, d8_pntr, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        seed_pts = parameters[0].valueAsText
        d8_pntr = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.trace_downslope_flowpaths(seed_pts, d8_pntr, output, esri_pntr, zero_background))
        return


class UnnestBasins(object):
    def __init__(self):
        self.label = "Unnest Basins"
        self.description = "Extract whole watersheds for a set of outlet points."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        pour_pts = arcpy.Parameter(
            displayName="Input Pour Points (Outlet) File",
            name="pour_pts",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        pour_pts.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, pour_pts, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        pour_pts = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        messages.addMessage(wbt.unnest_basins(d8_pntr, pour_pts, output, esri_pntr))
        return


class Watershed(object):
    def __init__(self):
        self.label = "Watershed"
        self.description = "Identifies the watershed, or drainage basin, draining to a set of target cells."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        pour_pts = arcpy.Parameter(
            displayName="Input Pour Points (Outlet) File",
            name="pour_pts",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        pour_pts.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, pour_pts, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        pour_pts = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        messages.addMessage(wbt.watershed(d8_pntr, pour_pts, output, esri_pntr))
        return


class ChangeVectorAnalysis(object):
    def __init__(self):
        self.label = "Change Vector Analysis"
        self.description = "Performs a change vector analysis on a two-date multi-spectral dataset."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        date1 = arcpy.Parameter(
            displayName="Earlier Date Input Files",
            name="date1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        date2 = arcpy.Parameter(
            displayName="Later Date Input Files",
            name="date2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        magnitude = arcpy.Parameter(
            displayName="Output Vector Magnitude File",
            name="magnitude",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        magnitude.filter.list = ["tif"]

        direction = arcpy.Parameter(
            displayName="Output Vector Direction File",
            name="direction",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        direction.filter.list = ["tif"]

        params = [date1, date2, magnitude, direction]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        date1 = parameters[0].valueAsText
        date2 = parameters[1].valueAsText
        magnitude = parameters[2].valueAsText
        direction = parameters[3].valueAsText
        messages.addMessage(wbt.change_vector_analysis(date1, date2, magnitude, direction))
        return


class Closing(object):
    def __init__(self):
        self.label = "Closing"
        self.description = "A closing is a mathematical morphology operation involving an erosion (min filter) of a dilation (max filter) set."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.closing(input, output, filterx, filtery))
        return


class CreateColourComposite(object):
    def __init__(self):
        self.label = "Create Colour Composite"
        self.description = "Creates a colour-composite image from three bands of multispectral imagery."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        red = arcpy.Parameter(
            displayName="Input Red Band Image File",
            name="red",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        green = arcpy.Parameter(
            displayName="Input Green Band Image File",
            name="green",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        blue = arcpy.Parameter(
            displayName="Input Blue Band Image File",
            name="blue",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        opacity = arcpy.Parameter(
            displayName="Input Opacity Band Image File (Optional)",
            name="opacity",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Colour Composite File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        enhance = arcpy.Parameter(
            displayName="Perform balance contrast enhancement?",
            name="enhance",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        enhance.value = "true"

        zeros = arcpy.Parameter(
            displayName="Treat zeros as nodata?",
            name="zeros",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        zeros.value = "false"

        params = [red, green, blue, opacity, output, enhance, zeros]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        red = parameters[0].valueAsText
        green = parameters[1].valueAsText
        blue = parameters[2].valueAsText
        opacity = parameters[3].valueAsText
        output = parameters[4].valueAsText
        enhance = parameters[5].valueAsText
        zeros = parameters[6].valueAsText
        messages.addMessage(wbt.create_colour_composite(red, green, blue, opacity, output, enhance, zeros))
        return


class FlipImage(object):
    def __init__(self):
        self.label = "Flip Image"
        self.description = "Reflects an image in the vertical or horizontal axis."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        direction = arcpy.Parameter(
            displayName="Direction",
            name="direction",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        direction.filter.type = "ValueList"
        direction.filter.list = ['vertical', 'horizontal', 'both']

        direction.value = "vertical"

        params = [input, output, direction]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        direction = parameters[2].valueAsText
        messages.addMessage(wbt.flip_image(input, output, direction))
        return


class IhsToRgb(object):
    def __init__(self):
        self.label = "Ihs To Rgb"
        self.description = "Converts intensity, hue, and saturation (IHS) images into red, green, and blue (RGB) images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        intensity = arcpy.Parameter(
            displayName="Input Intensity File",
            name="intensity",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        hue = arcpy.Parameter(
            displayName="Input Hue File",
            name="hue",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        saturation = arcpy.Parameter(
            displayName="Input Saturation File",
            name="saturation",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        red = arcpy.Parameter(
            displayName="Output Red Band File (optional; only if colour-composite not specified)",
            name="red",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        red.filter.list = ["tif"]

        green = arcpy.Parameter(
            displayName="Output Green Band File (optional; only if colour-composite not specified)",
            name="green",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        green.filter.list = ["tif"]

        blue = arcpy.Parameter(
            displayName="Output Blue Band File (optional; only if colour-composite not specified)",
            name="blue",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        blue.filter.list = ["tif"]

        output = arcpy.Parameter(
            displayName="Output Colour-Composite File (optional; only if individual bands not specified)",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [intensity, hue, saturation, red, green, blue, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        intensity = parameters[0].valueAsText
        hue = parameters[1].valueAsText
        saturation = parameters[2].valueAsText
        red = parameters[3].valueAsText
        green = parameters[4].valueAsText
        blue = parameters[5].valueAsText
        output = parameters[6].valueAsText
        messages.addMessage(wbt.ihs_to_rgb(intensity, hue, saturation, red, green, blue, output))
        return


class ImageStackProfile(object):
    def __init__(self):
        self.label = "Image Stack Profile"
        self.description = "Plots an image stack profile (i.e. signature) for a set of points and multispectral images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        points = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="points",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        points.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [inputs, points, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        points = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.image_stack_profile(inputs, points, output))
        return


class IntegralImage(object):
    def __init__(self):
        self.label = "Integral Image"
        self.description = "Transforms an input image (summed area table) into its integral image equivalent."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.integral_image(input, output))
        return


class KMeansClustering(object):
    def __init__(self):
        self.label = "K Means Clustering"
        self.description = "Performs a k-means clustering operation on a multi-spectral dataset."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Raster File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_html = arcpy.Parameter(
            displayName="Output HTML Report File",
            name="out_html",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_html.filter.list = ["html"]

        classes = arcpy.Parameter(
            displayName="Num. Classes (k)",
            name="classes",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        max_iterations = arcpy.Parameter(
            displayName="Max. Iterations",
            name="max_iterations",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        max_iterations.value = "10"

        class_change = arcpy.Parameter(
            displayName="Percent Class Change Threshold",
            name="class_change",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        class_change.value = "2.0"

        initialize = arcpy.Parameter(
            displayName="How to Initialize Cluster Centres?",
            name="initialize",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        initialize.filter.type = "ValueList"
        initialize.filter.list = ['diagonal', 'random']

        initialize.value = "diagonal"

        min_class_size = arcpy.Parameter(
            displayName="Min. Class Size",
            name="min_class_size",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        min_class_size.value = "10"

        params = [inputs, output, out_html, classes, max_iterations, class_change, initialize, min_class_size]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_html = parameters[2].valueAsText
        classes = parameters[3].valueAsText
        max_iterations = parameters[4].valueAsText
        class_change = parameters[5].valueAsText
        initialize = parameters[6].valueAsText
        min_class_size = parameters[7].valueAsText
        messages.addMessage(wbt.k_means_clustering(inputs, output, out_html, classes, max_iterations, class_change, initialize, min_class_size))
        return


class LineThinning(object):
    def __init__(self):
        self.label = "Line Thinning"
        self.description = "Performs line thinning a on Boolean raster image; intended to be used with the RemoveSpurs tool."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.line_thinning(input, output))
        return


class ModifiedKMeansClustering(object):
    def __init__(self):
        self.label = "Modified K Means Clustering"
        self.description = "Performs a modified k-means clustering operation on a multi-spectral dataset."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Raster File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_html = arcpy.Parameter(
            displayName="Output HTML Report File",
            name="out_html",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_html.filter.list = ["html"]

        start_clusters = arcpy.Parameter(
            displayName="Initial Num. of Clusters",
            name="start_clusters",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        start_clusters.value = "1000"

        merger_dist = arcpy.Parameter(
            displayName="Cluster Merger Distance",
            name="merger_dist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        max_iterations = arcpy.Parameter(
            displayName="Max. Iterations",
            name="max_iterations",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        max_iterations.value = "10"

        class_change = arcpy.Parameter(
            displayName="Percent Class Change Threshold",
            name="class_change",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        class_change.value = "2.0"

        params = [inputs, output, out_html, start_clusters, merger_dist, max_iterations, class_change]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_html = parameters[2].valueAsText
        start_clusters = parameters[3].valueAsText
        merger_dist = parameters[4].valueAsText
        max_iterations = parameters[5].valueAsText
        class_change = parameters[6].valueAsText
        messages.addMessage(wbt.modified_k_means_clustering(inputs, output, out_html, start_clusters, merger_dist, max_iterations, class_change))
        return


class Mosaic(object):
    def __init__(self):
        self.label = "Mosaic"
        self.description = "Mosaics two or more images together."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        method = arcpy.Parameter(
            displayName="Resampling Method",
            name="method",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        method.filter.type = "ValueList"
        method.filter.list = ['nn', 'bilinear', 'cc']

        method.value = "cc"

        params = [inputs, output, method]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        method = parameters[2].valueAsText
        messages.addMessage(wbt.mosaic(inputs, output, method))
        return


class MosaicWithFeathering(object):
    def __init__(self):
        self.label = "Mosaic With Feathering"
        self.description = "Mosaics two images together using a feathering technique in overlapping areas to reduce edge-effects."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File To Modify",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input Reference File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        method = arcpy.Parameter(
            displayName="Resampling Method",
            name="method",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        method.filter.type = "ValueList"
        method.filter.list = ['nn', 'bilinear', 'cc']

        method.value = "cc"

        weight = arcpy.Parameter(
            displayName="Distance Weight",
            name="weight",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        weight.value = "4.0"

        params = [input1, input2, output, method, weight]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        method = parameters[3].valueAsText
        weight = parameters[4].valueAsText
        messages.addMessage(wbt.mosaic_with_feathering(input1, input2, output, method, weight))
        return


class NormalizedDifferenceVegetationIndex(object):
    def __init__(self):
        self.label = "Normalized Difference Vegetation Index"
        self.description = "Calculates the normalized difference vegetation index (NDVI) from near-infrared and red imagery."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        nir = arcpy.Parameter(
            displayName="Input Near-Infrared File",
            name="nir",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        red = arcpy.Parameter(
            displayName="Input Red File",
            name="red",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (%)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        osavi = arcpy.Parameter(
            displayName="Use the optimized soil-adjusted veg index (OSAVI)?",
            name="osavi",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [nir, red, output, clip, osavi]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        nir = parameters[0].valueAsText
        red = parameters[1].valueAsText
        output = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        osavi = parameters[4].valueAsText
        messages.addMessage(wbt.normalized_difference_vegetation_index(nir, red, output, clip, osavi))
        return


class Opening(object):
    def __init__(self):
        self.label = "Opening"
        self.description = "An opening is a mathematical morphology operation involving a dilation (max filter) of an erosion (min filter) set."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.opening(input, output, filterx, filtery))
        return


class RemoveSpurs(object):
    def __init__(self):
        self.label = "Remove Spurs"
        self.description = "Removes the spurs (pruning operation) from a Boolean line image; intended to be used on the output of the LineThinning tool."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        iterations = arcpy.Parameter(
            displayName="Maximum Iterations",
            name="iterations",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        iterations.value = "10"

        params = [input, output, iterations]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        iterations = parameters[2].valueAsText
        messages.addMessage(wbt.remove_spurs(input, output, iterations))
        return


class Resample(object):
    def __init__(self):
        self.label = "Resample"
        self.description = "Resamples one or more input images into a destination image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        destination = arcpy.Parameter(
            displayName="Destination File",
            name="destination",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        method = arcpy.Parameter(
            displayName="Resampling Method",
            name="method",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        method.filter.type = "ValueList"
        method.filter.list = ['nn', 'bilinear', 'cc']

        method.value = "cc"

        params = [inputs, destination, method]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        destination = parameters[1].valueAsText
        method = parameters[2].valueAsText
        messages.addMessage(wbt.resample(inputs, destination, method))
        return


class RgbToIhs(object):
    def __init__(self):
        self.label = "Rgb To Ihs"
        self.description = "Converts red, green, and blue (RGB) images into intensity, hue, and saturation (IHS) images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        red = arcpy.Parameter(
            displayName="Input Red Band File (optional; only if colour-composite not specified)",
            name="red",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        green = arcpy.Parameter(
            displayName="Input Green Band File (optional; only if colour-composite not specified)",
            name="green",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        blue = arcpy.Parameter(
            displayName="Input Blue Band File (optional; only if colour-composite not specified)",
            name="blue",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        composite = arcpy.Parameter(
            displayName="Input Colour-Composite Image File (optional; only if individual bands not specified)",
            name="composite",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        intensity = arcpy.Parameter(
            displayName="Output Intensity File",
            name="intensity",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        intensity.filter.list = ["tif"]

        hue = arcpy.Parameter(
            displayName="Output Hue File",
            name="hue",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        hue.filter.list = ["tif"]

        saturation = arcpy.Parameter(
            displayName="Output Saturation File",
            name="saturation",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        saturation.filter.list = ["tif"]

        params = [red, green, blue, composite, intensity, hue, saturation]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        red = parameters[0].valueAsText
        green = parameters[1].valueAsText
        blue = parameters[2].valueAsText
        composite = parameters[3].valueAsText
        intensity = parameters[4].valueAsText
        hue = parameters[5].valueAsText
        saturation = parameters[6].valueAsText
        messages.addMessage(wbt.rgb_to_ihs(red, green, blue, composite, intensity, hue, saturation))
        return


class SplitColourComposite(object):
    def __init__(self):
        self.label = "Split Colour Composite"
        self.description = "This tool splits an RGB colour composite image into seperate multispectral images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Colour Composite Image File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.split_colour_composite(input, output))
        return


class ThickenRasterLine(object):
    def __init__(self):
        self.label = "Thicken Raster Line"
        self.description = "Thickens single-cell wide lines within a raster image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.thicken_raster_line(input, output))
        return


class TophatTransform(object):
    def __init__(self):
        self.label = "Tophat Transform"
        self.description = "Performs either a white or black top-hat transform on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        variant = arcpy.Parameter(
            displayName="Variant",
            name="variant",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        variant.filter.type = "ValueList"
        variant.filter.list = ['white', 'black']

        variant.value = "white"

        params = [input, output, filterx, filtery, variant]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        variant = parameters[4].valueAsText
        messages.addMessage(wbt.tophat_transform(input, output, filterx, filtery, variant))
        return


class WriteFunctionMemoryInsertion(object):
    def __init__(self):
        self.label = "Write Function Memory Insertion"
        self.description = "Performs a write function memory insertion for single-band multi-date change detection."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="First Date Input File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Second Date Input File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input3 = arcpy.Parameter(
            displayName="Third Date Input File (Optional)",
            name="input3",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, input3, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        input3 = parameters[2].valueAsText
        output = parameters[3].valueAsText
        messages.addMessage(wbt.write_function_memory_insertion(input1, input2, input3, output))
        return


class AdaptiveFilter(object):
    def __init__(self):
        self.label = "Adaptive Filter"
        self.description = "Performs an adaptive filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        threshold = arcpy.Parameter(
            displayName="Difference From Mean Threshold (# Std. Dev.)",
            name="threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        threshold.value = "2.0"

        params = [input, output, filterx, filtery, threshold]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        threshold = parameters[4].valueAsText
        messages.addMessage(wbt.adaptive_filter(input, output, filterx, filtery, threshold))
        return


class BilateralFilter(object):
    def __init__(self):
        self.label = "Bilateral Filter"
        self.description = "A bilateral filter is an edge-preserving smoothing filter introduced by Tomasi and Manduchi (1998)."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sigma_dist = arcpy.Parameter(
            displayName="Distance Standard Deviation (pixels)",
            name="sigma_dist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma_dist.value = "0.75"

        sigma_int = arcpy.Parameter(
            displayName="Intensity Standard Deviation (intensity units)",
            name="sigma_int",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma_int.value = "1.0"

        params = [input, output, sigma_dist, sigma_int]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        sigma_dist = parameters[2].valueAsText
        sigma_int = parameters[3].valueAsText
        messages.addMessage(wbt.bilateral_filter(input, output, sigma_dist, sigma_int))
        return


class ConservativeSmoothingFilter(object):
    def __init__(self):
        self.label = "Conservative Smoothing Filter"
        self.description = "Performs a conservative-smoothing filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "3"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "3"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.conservative_smoothing_filter(input, output, filterx, filtery))
        return


class CornerDetection(object):
    def __init__(self):
        self.label = "Corner Detection"
        self.description = "Identifies corner patterns in boolean images using hit-and-miss pattern mattching."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.corner_detection(input, output))
        return


class DiffOfGaussianFilter(object):
    def __init__(self):
        self.label = "Diff Of Gaussian Filter"
        self.description = "Performs a Difference of Gaussian (DoG) filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sigma1 = arcpy.Parameter(
            displayName="Sigma 1 (pixels)",
            name="sigma1",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma1.value = "2.0"

        sigma2 = arcpy.Parameter(
            displayName="Sigma 2 (pixels)",
            name="sigma2",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma2.value = "4.0"

        params = [input, output, sigma1, sigma2]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        sigma1 = parameters[2].valueAsText
        sigma2 = parameters[3].valueAsText
        messages.addMessage(wbt.diff_of_gaussian_filter(input, output, sigma1, sigma2))
        return


class DiversityFilter(object):
    def __init__(self):
        self.label = "Diversity Filter"
        self.description = "Assigns each cell in the output grid the number of different values in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.diversity_filter(input, output, filterx, filtery))
        return


class EdgePreservingMeanFilter(object):
    def __init__(self):
        self.label = "Edge Preserving Mean Filter"
        self.description = "Performs a simple edge-preserving mean filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filter = arcpy.Parameter(
            displayName="Filter Size",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "11"

        threshold = arcpy.Parameter(
            displayName="Value Difference Threshold",
            name="threshold",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [input, output, filter, threshold]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filter = parameters[2].valueAsText
        threshold = parameters[3].valueAsText
        messages.addMessage(wbt.edge_preserving_mean_filter(input, output, filter, threshold))
        return


class EmbossFilter(object):
    def __init__(self):
        self.label = "Emboss Filter"
        self.description = "Performs an emboss filter on an image, similar to a hillshade operation."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        direction = arcpy.Parameter(
            displayName="Direction",
            name="direction",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        direction.filter.type = "ValueList"
        direction.filter.list = ['n', 's', 'e', 'w', 'ne', 'se', 'nw', 'sw']

        direction.value = "n"

        clip = arcpy.Parameter(
            displayName="Percent to clip the distribution tails",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, direction, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        direction = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        messages.addMessage(wbt.emboss_filter(input, output, direction, clip))
        return


class FastAlmostGaussianFilter(object):
    def __init__(self):
        self.label = "Fast Almost Gaussian Filter"
        self.description = "Performs a fast approximate Gaussian filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sigma = arcpy.Parameter(
            displayName="Standard Deviation (pixels)",
            name="sigma",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        sigma.value = "1.8"

        params = [input, output, sigma]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        sigma = parameters[2].valueAsText
        messages.addMessage(wbt.fast_almost_gaussian_filter(input, output, sigma))
        return


class GaussianFilter(object):
    def __init__(self):
        self.label = "Gaussian Filter"
        self.description = "Performs a Gaussian filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sigma = arcpy.Parameter(
            displayName="Standard Deviation (pixels)",
            name="sigma",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        sigma.value = "0.75"

        params = [input, output, sigma]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        sigma = parameters[2].valueAsText
        messages.addMessage(wbt.gaussian_filter(input, output, sigma))
        return


class HighPassFilter(object):
    def __init__(self):
        self.label = "High Pass Filter"
        self.description = "Performs a high-pass filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.high_pass_filter(input, output, filterx, filtery))
        return


class HighPassMedianFilter(object):
    def __init__(self):
        self.label = "High Pass Median Filter"
        self.description = "Performs a high pass median filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        sig_digits = arcpy.Parameter(
            displayName="Number of Significant Digits",
            name="sig_digits",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        sig_digits.value = "2"

        params = [input, output, filterx, filtery, sig_digits]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        sig_digits = parameters[4].valueAsText
        messages.addMessage(wbt.high_pass_median_filter(input, output, filterx, filtery, sig_digits))
        return


class KNearestMeanFilter(object):
    def __init__(self):
        self.label = "K Nearest Mean Filter"
        self.description = "A k-nearest mean filter is a type of edge-preserving smoothing filter."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        k = arcpy.Parameter(
            displayName="K-value (pixels)",
            name="k",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        k.value = "5"

        params = [input, output, filterx, filtery, k]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        k = parameters[4].valueAsText
        messages.addMessage(wbt.k_nearest_mean_filter(input, output, filterx, filtery, k))
        return


class LaplacianFilter(object):
    def __init__(self):
        self.label = "Laplacian Filter"
        self.description = "Performs a Laplacian filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        variant = arcpy.Parameter(
            displayName="Variant",
            name="variant",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        variant.filter.type = "ValueList"
        variant.filter.list = ['3x3(1)', '3x3(2)', '3x3(3)', '3x3(4)', '5x5(1)', '5x5(2)']

        variant.value = "3x3(1)"

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (%)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, variant, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        variant = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        messages.addMessage(wbt.laplacian_filter(input, output, variant, clip))
        return


class LaplacianOfGaussianFilter(object):
    def __init__(self):
        self.label = "Laplacian Of Gaussian Filter"
        self.description = "Performs a Laplacian-of-Gaussian (LoG) filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sigma = arcpy.Parameter(
            displayName="Standard Deviation (Pixels)",
            name="sigma",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma.value = "0.75"

        params = [input, output, sigma]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        sigma = parameters[2].valueAsText
        messages.addMessage(wbt.laplacian_of_gaussian_filter(input, output, sigma))
        return


class LeeFilter(object):
    def __init__(self):
        self.label = "Lee Filter"
        self.description = "Performs a Lee (Sigma) smoothing filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        sigma = arcpy.Parameter(
            displayName="Sigma",
            name="sigma",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma.value = "10.0"

        m = arcpy.Parameter(
            displayName="M-value",
            name="m",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        m.value = "5.0"

        params = [input, output, filterx, filtery, sigma, m]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        sigma = parameters[4].valueAsText
        m = parameters[5].valueAsText
        messages.addMessage(wbt.lee_filter(input, output, filterx, filtery, sigma, m))
        return


class LineDetectionFilter(object):
    def __init__(self):
        self.label = "Line Detection Filter"
        self.description = "Performs a line-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        variant = arcpy.Parameter(
            displayName="Variant",
            name="variant",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        variant.filter.type = "ValueList"
        variant.filter.list = ['vertical', 'horizontal', '45', '135']

        variant.value = "vertical"

        absvals = arcpy.Parameter(
            displayName="Output absolute values?",
            name="absvals",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (%)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, variant, absvals, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        variant = parameters[2].valueAsText
        absvals = parameters[3].valueAsText
        clip = parameters[4].valueAsText
        messages.addMessage(wbt.line_detection_filter(input, output, variant, absvals, clip))
        return


class MajorityFilter(object):
    def __init__(self):
        self.label = "Majority Filter"
        self.description = "Assigns each cell in the output grid the most frequently occurring value (mode) in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.majority_filter(input, output, filterx, filtery))
        return


class MaximumFilter(object):
    def __init__(self):
        self.label = "Maximum Filter"
        self.description = "Assigns each cell in the output grid the maximum value in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.maximum_filter(input, output, filterx, filtery))
        return


class MeanFilter(object):
    def __init__(self):
        self.label = "Mean Filter"
        self.description = "Performs a mean filter (low-pass filter) on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "3"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "3"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.mean_filter(input, output, filterx, filtery))
        return


class MedianFilter(object):
    def __init__(self):
        self.label = "Median Filter"
        self.description = "Performs a median filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        sig_digits = arcpy.Parameter(
            displayName="Number of Significant Digits",
            name="sig_digits",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        sig_digits.value = "2"

        params = [input, output, filterx, filtery, sig_digits]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        sig_digits = parameters[4].valueAsText
        messages.addMessage(wbt.median_filter(input, output, filterx, filtery, sig_digits))
        return


class MinimumFilter(object):
    def __init__(self):
        self.label = "Minimum Filter"
        self.description = "Assigns each cell in the output grid the minimum value in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.minimum_filter(input, output, filterx, filtery))
        return


class OlympicFilter(object):
    def __init__(self):
        self.label = "Olympic Filter"
        self.description = "Performs an olympic smoothing filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.olympic_filter(input, output, filterx, filtery))
        return


class PercentileFilter(object):
    def __init__(self):
        self.label = "Percentile Filter"
        self.description = "Performs a percentile filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        sig_digits = arcpy.Parameter(
            displayName="Number of Significant Digits",
            name="sig_digits",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        sig_digits.value = "2"

        params = [input, output, filterx, filtery, sig_digits]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        sig_digits = parameters[4].valueAsText
        messages.addMessage(wbt.percentile_filter(input, output, filterx, filtery, sig_digits))
        return


class PrewittFilter(object):
    def __init__(self):
        self.label = "Prewitt Filter"
        self.description = "Performs a Prewitt edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (Percent)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        clip = parameters[2].valueAsText
        messages.addMessage(wbt.prewitt_filter(input, output, clip))
        return


class RangeFilter(object):
    def __init__(self):
        self.label = "Range Filter"
        self.description = "Assigns each cell in the output grid the range of values in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.range_filter(input, output, filterx, filtery))
        return


class RobertsCrossFilter(object):
    def __init__(self):
        self.label = "Roberts Cross Filter"
        self.description = "Performs a Robert's cross edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (Percent)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        clip = parameters[2].valueAsText
        messages.addMessage(wbt.roberts_cross_filter(input, output, clip))
        return


class ScharrFilter(object):
    def __init__(self):
        self.label = "Scharr Filter"
        self.description = "Performs a Scharr edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (Percent)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        clip = parameters[2].valueAsText
        messages.addMessage(wbt.scharr_filter(input, output, clip))
        return


class SobelFilter(object):
    def __init__(self):
        self.label = "Sobel Filter"
        self.description = "Performs a Sobel edge-detection filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        variant = arcpy.Parameter(
            displayName="Variant",
            name="variant",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        variant.filter.type = "ValueList"
        variant.filter.list = ['3x3', '5x5']

        variant.value = "3x3"

        clip = arcpy.Parameter(
            displayName="Clip Tails (%)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        params = [input, output, variant, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        variant = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        messages.addMessage(wbt.sobel_filter(input, output, variant, clip))
        return


class StandardDeviationFilter(object):
    def __init__(self):
        self.label = "Standard Deviation Filter"
        self.description = "Assigns each cell in the output grid the standard deviation of values in a moving window centred on each grid cell in the input raster."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.standard_deviation_filter(input, output, filterx, filtery))
        return


class TotalFilter(object):
    def __init__(self):
        self.label = "Total Filter"
        self.description = "Performs a total filter on an input image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        filterx = arcpy.Parameter(
            displayName="Filter X-Dimension",
            name="filterx",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filterx.value = "11"

        filtery = arcpy.Parameter(
            displayName="Filter Y-Dimension",
            name="filtery",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filtery.value = "11"

        params = [input, output, filterx, filtery]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        filterx = parameters[2].valueAsText
        filtery = parameters[3].valueAsText
        messages.addMessage(wbt.total_filter(input, output, filterx, filtery))
        return


class UnsharpMasking(object):
    def __init__(self):
        self.label = "Unsharp Masking"
        self.description = "An image sharpening technique that enhances edges."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        sigma = arcpy.Parameter(
            displayName="Standard Deviation (pixels)",
            name="sigma",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        sigma.value = "0.75"

        amount = arcpy.Parameter(
            displayName="Amount (%)",
            name="amount",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        amount.value = "100.0"

        threshold = arcpy.Parameter(
            displayName="Threshold",
            name="threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        threshold.value = "0.0"

        params = [input, output, sigma, amount, threshold]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        sigma = parameters[2].valueAsText
        amount = parameters[3].valueAsText
        threshold = parameters[4].valueAsText
        messages.addMessage(wbt.unsharp_masking(input, output, sigma, amount, threshold))
        return


class UserDefinedWeightsFilter(object):
    def __init__(self):
        self.label = "User Defined Weights Filter"
        self.description = "Performs a user-defined weights filter on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        weights = arcpy.Parameter(
            displayName="Input Weights File",
            name="weights",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        weights.filter.list = ["csv"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        center = arcpy.Parameter(
            displayName="Kernel Center",
            name="center",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        center.filter.type = "ValueList"
        center.filter.list = ['center', 'upper-left', 'upper-right', 'lower-left', 'lower-right']

        center.value = "center"

        normalize = arcpy.Parameter(
            displayName="Normalize kernel weights?",
            name="normalize",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        normalize.value = "false"

        params = [input, weights, output, center, normalize]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        weights = parameters[1].valueAsText
        output = parameters[2].valueAsText
        center = parameters[3].valueAsText
        normalize = parameters[4].valueAsText
        messages.addMessage(wbt.user_defined_weights_filter(input, weights, output, center, normalize))
        return


class BalanceContrastEnhancement(object):
    def __init__(self):
        self.label = "Balance Contrast Enhancement"
        self.description = "Performs a balance contrast enhancement on a colour-composite image of multispectral data."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Colour Composite Image File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        band_mean = arcpy.Parameter(
            displayName="Band Mean Value",
            name="band_mean",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        band_mean.value = "100.0"

        params = [input, output, band_mean]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        band_mean = parameters[2].valueAsText
        messages.addMessage(wbt.balance_contrast_enhancement(input, output, band_mean))
        return


class CorrectVignetting(object):
    def __init__(self):
        self.label = "Correct Vignetting"
        self.description = "Corrects the darkening of images towards corners."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        pp = arcpy.Parameter(
            displayName="Input Principal Point File",
            name="pp",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        pp.filter.list = ["Point"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        focal_length = arcpy.Parameter(
            displayName="Camera Focal Length (mm)",
            name="focal_length",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        focal_length.value = "304.8"

        image_width = arcpy.Parameter(
            displayName="Distance Between Left-Right Edges (mm)",
            name="image_width",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        image_width.value = "228.6"

        n = arcpy.Parameter(
            displayName="n Parameter",
            name="n",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        n.value = "4.0"

        params = [input, pp, output, focal_length, image_width, n]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        pp = parameters[1].valueAsText
        output = parameters[2].valueAsText
        focal_length = parameters[3].valueAsText
        image_width = parameters[4].valueAsText
        n = parameters[5].valueAsText
        messages.addMessage(wbt.correct_vignetting(input, pp, output, focal_length, image_width, n))
        return


class DirectDecorrelationStretch(object):
    def __init__(self):
        self.label = "Direct Decorrelation Stretch"
        self.description = "Performs a direct decorrelation stretch enhancement on a colour-composite image of multispectral data."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Colour Composite Image File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        k = arcpy.Parameter(
            displayName="Achromatic Factor (0-1)",
            name="k",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        k.value = "0.5"

        clip = arcpy.Parameter(
            displayName="Percent to clip the upper tail",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "1.0"

        params = [input, output, k, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        k = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        messages.addMessage(wbt.direct_decorrelation_stretch(input, output, k, clip))
        return


class GammaCorrection(object):
    def __init__(self):
        self.label = "Gamma Correction"
        self.description = "Performs a gamma correction on an input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        gamma = arcpy.Parameter(
            displayName="Gamma Value",
            name="gamma",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        gamma.value = "0.5"

        params = [input, output, gamma]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        gamma = parameters[2].valueAsText
        messages.addMessage(wbt.gamma_correction(input, output, gamma))
        return


class GaussianContrastStretch(object):
    def __init__(self):
        self.label = "Gaussian Contrast Stretch"
        self.description = "Performs a Gaussian contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        num_tones = arcpy.Parameter(
            displayName="Number of Tones",
            name="num_tones",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_tones.value = "256"

        params = [input, output, num_tones]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_tones = parameters[2].valueAsText
        messages.addMessage(wbt.gaussian_contrast_stretch(input, output, num_tones))
        return


class HistogramEqualization(object):
    def __init__(self):
        self.label = "Histogram Equalization"
        self.description = "Performs a histogram equalization contrast enhancment on an image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        num_tones = arcpy.Parameter(
            displayName="Number of Tones",
            name="num_tones",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_tones.value = "256"

        params = [input, output, num_tones]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_tones = parameters[2].valueAsText
        messages.addMessage(wbt.histogram_equalization(input, output, num_tones))
        return


class HistogramMatching(object):
    def __init__(self):
        self.label = "Histogram Matching"
        self.description = "Alters the statistical distribution of a raster image matching it to a specified PDF."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        histo_file = arcpy.Parameter(
            displayName="Input Probability Distribution Function (PDF) Text File",
            name="histo_file",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, histo_file, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        histo_file = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.histogram_matching(input, histo_file, output))
        return


class HistogramMatchingTwoImages(object):
    def __init__(self):
        self.label = "Histogram Matching Two Images"
        self.description = "This tool alters the cumulative distribution function of a raster image to that of another image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File To Modify",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input Reference File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.histogram_matching_two_images(input1, input2, output))
        return


class MinMaxContrastStretch(object):
    def __init__(self):
        self.label = "Min Max Contrast Stretch"
        self.description = "Performs a min-max contrast stretch on an input greytone image."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        min_val = arcpy.Parameter(
            displayName="Lower Tail Clip Value",
            name="min_val",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        max_val = arcpy.Parameter(
            displayName="Upper Tail Clip Value",
            name="max_val",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        num_tones = arcpy.Parameter(
            displayName="Number of Tones",
            name="num_tones",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_tones.value = "256"

        params = [input, output, min_val, max_val, num_tones]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        min_val = parameters[2].valueAsText
        max_val = parameters[3].valueAsText
        num_tones = parameters[4].valueAsText
        messages.addMessage(wbt.min_max_contrast_stretch(input, output, min_val, max_val, num_tones))
        return


class PanchromaticSharpening(object):
    def __init__(self):
        self.label = "Panchromatic Sharpening"
        self.description = "Increases the spatial resolution of image data by combining multispectral bands with panchromatic data."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        red = arcpy.Parameter(
            displayName="Input Red Band File (optional; only if colour-composite not specified)",
            name="red",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        green = arcpy.Parameter(
            displayName="Input Green Band File (optional; only if colour-composite not specified)",
            name="green",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        blue = arcpy.Parameter(
            displayName="Input Blue Band File (optional; only if colour-composite not specified)",
            name="blue",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        composite = arcpy.Parameter(
            displayName="Input Colour-Composite Image File (optional; only if individual bands not specified)",
            name="composite",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        pan = arcpy.Parameter(
            displayName="Input Panchromatic Band File",
            name="pan",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Colour Composite File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        method = arcpy.Parameter(
            displayName="Pan-Sharpening Method",
            name="method",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        method.filter.type = "ValueList"
        method.filter.list = ['brovey', 'ihs']

        method.value = "brovey"

        params = [red, green, blue, composite, pan, output, method]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        red = parameters[0].valueAsText
        green = parameters[1].valueAsText
        blue = parameters[2].valueAsText
        composite = parameters[3].valueAsText
        pan = parameters[4].valueAsText
        output = parameters[5].valueAsText
        method = parameters[6].valueAsText
        messages.addMessage(wbt.panchromatic_sharpening(red, green, blue, composite, pan, output, method))
        return


class PercentageContrastStretch(object):
    def __init__(self):
        self.label = "Percentage Contrast Stretch"
        self.description = "Performs a percentage linear contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (%)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "1.0"

        tail = arcpy.Parameter(
            displayName="Tail",
            name="tail",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        tail.filter.type = "ValueList"
        tail.filter.list = ['upper', 'lower', 'both']

        tail.value = "both"

        num_tones = arcpy.Parameter(
            displayName="Number of Tones",
            name="num_tones",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_tones.value = "256"

        params = [input, output, clip, tail, num_tones]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        clip = parameters[2].valueAsText
        tail = parameters[3].valueAsText
        num_tones = parameters[4].valueAsText
        messages.addMessage(wbt.percentage_contrast_stretch(input, output, clip, tail, num_tones))
        return


class SigmoidalContrastStretch(object):
    def __init__(self):
        self.label = "Sigmoidal Contrast Stretch"
        self.description = "Performs a sigmoidal contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        cutoff = arcpy.Parameter(
            displayName="Cutoff Value (0.0 - 0.95)",
            name="cutoff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        cutoff.value = "0.0"

        gain = arcpy.Parameter(
            displayName="Gain Value",
            name="gain",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        gain.value = "1.0"

        num_tones = arcpy.Parameter(
            displayName="Number of Tones",
            name="num_tones",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_tones.value = "256"

        params = [input, output, cutoff, gain, num_tones]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        cutoff = parameters[2].valueAsText
        gain = parameters[3].valueAsText
        num_tones = parameters[4].valueAsText
        messages.addMessage(wbt.sigmoidal_contrast_stretch(input, output, cutoff, gain, num_tones))
        return


class StandardDeviationContrastStretch(object):
    def __init__(self):
        self.label = "Standard Deviation Contrast Stretch"
        self.description = "Performs a standard-deviation contrast stretch on input images."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        stdev = arcpy.Parameter(
            displayName="Standard Deviation Threshold",
            name="stdev",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        stdev.value = "2.0"

        num_tones = arcpy.Parameter(
            displayName="Number of Tones",
            name="num_tones",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_tones.value = "256"

        params = [input, output, stdev, num_tones]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        stdev = parameters[2].valueAsText
        num_tones = parameters[3].valueAsText
        messages.addMessage(wbt.standard_deviation_contrast_stretch(input, output, stdev, num_tones))
        return


class ClassifyOverlapPoints(object):
    def __init__(self):
        self.label = "Classify Overlap Points"
        self.description = "Classifies or filters LAS points in regions of overlapping flight lines."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        resolution = arcpy.Parameter(
            displayName="Sample Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "2.0"

        filter = arcpy.Parameter(
            displayName="Filter out points from overlapping flightlines?",
            name="filter",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        filter.value = "false"

        params = [input, output, resolution, filter]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        resolution = parameters[2].valueAsText
        filter = parameters[3].valueAsText
        messages.addMessage(wbt.classify_overlap_points(input, output, resolution, filter))
        return


class ClipLidarToPolygon(object):
    def __init__(self):
        self.label = "Clip Lidar To Polygon"
        self.description = "Clips a LiDAR point cloud to a vector polygon or polygons."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        polygons = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="polygons",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        polygons.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        params = [input, polygons, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        polygons = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.clip_lidar_to_polygon(input, polygons, output))
        return


class ErasePolygonFromLidar(object):
    def __init__(self):
        self.label = "Erase Polygon From Lidar"
        self.description = "Erases (cuts out) a vector polygon or polygons from a LiDAR point cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        polygons = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="polygons",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        polygons.filter.list = ["Polygon"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        params = [input, polygons, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        polygons = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.erase_polygon_from_lidar(input, polygons, output))
        return


class FilterLidarScanAngles(object):
    def __init__(self):
        self.label = "Filter Lidar Scan Angles"
        self.description = "Removes points in a LAS file with scan angles greater than a threshold."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        threshold = arcpy.Parameter(
            displayName="Threshold (degrees)",
            name="threshold",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        params = [input, output, threshold]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        threshold = parameters[2].valueAsText
        messages.addMessage(wbt.filter_lidar_scan_angles(input, output, threshold))
        return


class FindFlightlineEdgePoints(object):
    def __init__(self):
        self.label = "Find Flightline Edge Points"
        self.description = "Identifies points along a flightline's edge in a LAS file."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.find_flightline_edge_points(input, output))
        return


class FlightlineOverlap(object):
    def __init__(self):
        self.label = "Flightline Overlap"
        self.description = "Reads a LiDAR (LAS) point file and outputs a raster containing the number of overlapping flight lines in each grid cell."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        params = [input, output, resolution]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        resolution = parameters[2].valueAsText
        messages.addMessage(wbt.flightline_overlap(input, output, resolution))
        return


class LasToAscii(object):
    def __init__(self):
        self.label = "Las To Ascii"
        self.description = "Converts one or more LAS files into ASCII text files."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input LiDAR Files",
            name="inputs",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        inputs.filter.list = ["las", "zip"]

        params = [inputs]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        messages.addMessage(wbt.las_to_ascii(inputs))
        return


class LasToMultipointShapefile(object):
    def __init__(self):
        self.label = "Las To Multipoint Shapefile"
        self.description = "Converts one or more LAS files into MultipointZ vector Shapefiles. When the input parameter is not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.las_to_multipoint_shapefile(input))
        return


class LasToShapefile(object):
    def __init__(self):
        self.label = "Las To Shapefile"
        self.description = "Converts one or more LAS files into a vector Shapefile of POINT ShapeType."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.las_to_shapefile(input))
        return


class LidarBlockMaximum(object):
    def __init__(self):
        self.label = "Lidar Block Maximum"
        self.description = "Creates a block-maximum raster from an input LAS file. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        params = [input, output, resolution]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        resolution = parameters[2].valueAsText
        messages.addMessage(wbt.lidar_block_maximum(input, output, resolution))
        return


class LidarBlockMinimum(object):
    def __init__(self):
        self.label = "Lidar Block Minimum"
        self.description = "Creates a block-minimum raster from an input LAS file. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        params = [input, output, resolution]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        resolution = parameters[2].valueAsText
        messages.addMessage(wbt.lidar_block_minimum(input, output, resolution))
        return


class LidarClassifySubset(object):
    def __init__(self):
        self.label = "Lidar Classify Subset"
        self.description = "Classifies the values in one LiDAR point cloud that correpond with points in a subset cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        base = arcpy.Parameter(
            displayName="Input Base LiDAR File",
            name="base",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        base.filter.list = ["las", "zip"]

        subset = arcpy.Parameter(
            displayName="Input Subset LiDAR File",
            name="subset",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        subset.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output LiDAR File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        subset_class = arcpy.Parameter(
            displayName="Subset Point Class Value",
            name="subset_class",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        nonsubset_class = arcpy.Parameter(
            displayName="Non-Subset Point Class Value (Optional)",
            name="nonsubset_class",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [base, subset, output, subset_class, nonsubset_class]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        base = parameters[0].valueAsText
        subset = parameters[1].valueAsText
        output = parameters[2].valueAsText
        subset_class = parameters[3].valueAsText
        nonsubset_class = parameters[4].valueAsText
        messages.addMessage(wbt.lidar_classify_subset(base, subset, output, subset_class, nonsubset_class))
        return


class LidarColourize(object):
    def __init__(self):
        self.label = "Lidar Colourize"
        self.description = "Adds the red-green-blue colour fields of a LiDAR (LAS) file based on an input image."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        in_lidar = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="in_lidar",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        in_lidar.filter.list = ["las", "zip"]

        in_image = arcpy.Parameter(
            displayName="Input Colour Image File",
            name="in_image",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output LiDAR File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        params = [in_lidar, in_image, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        in_lidar = parameters[0].valueAsText
        in_image = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.lidar_colourize(in_lidar, in_image, output))
        return


class LidarConstructVectorTin(object):
    def __init__(self):
        self.label = "Lidar Construct Vector Tin"
        self.description = "Creates a vector triangular irregular network (TIN) fitted to LiDAR points."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        returns = arcpy.Parameter(
            displayName="Point Returns Included",
            name="returns",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        returns.filter.type = "ValueList"
        returns.filter.list = ['all', 'last', 'first']

        returns.value = "all"

        exclude_cls = arcpy.Parameter(
            displayName="Exclusion Classes (0-18, based on LAS spec; e.g. 3,4,5,6,7)",
            name="exclude_cls",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        minz = arcpy.Parameter(
            displayName="Minimum Elevation Value (optional)",
            name="minz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxz = arcpy.Parameter(
            displayName="Maximum Elevation Value (optional)",
            name="maxz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, returns, exclude_cls, minz, maxz]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        returns = parameters[2].valueAsText
        exclude_cls = parameters[3].valueAsText
        minz = parameters[4].valueAsText
        maxz = parameters[5].valueAsText
        messages.addMessage(wbt.lidar_construct_vector_tin(input, output, returns, exclude_cls, minz, maxz))
        return


class LidarElevationSlice(object):
    def __init__(self):
        self.label = "Lidar Elevation Slice"
        self.description = "Outputs all of the points within a LiDAR (LAS) point file that lie between a specified elevation range."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        minz = arcpy.Parameter(
            displayName="Minimum Elevation Value",
            name="minz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxz = arcpy.Parameter(
            displayName="Maximum Elevation Value",
            name="maxz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        class1 = arcpy.Parameter(
            displayName="Retain but reclass points outside the specified elevation range?",
            name="class1",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        inclassval = arcpy.Parameter(
            displayName="Class Value Assigned to Points Within Range (Optional)",
            name="inclassval",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        inclassval.value = "2"

        outclassval = arcpy.Parameter(
            displayName="Class Value Assigned to Points Outside Range (Optional)",
            name="outclassval",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        outclassval.value = "1"

        params = [input, output, minz, maxz, class1, inclassval, outclassval]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        minz = parameters[2].valueAsText
        maxz = parameters[3].valueAsText
        class1 = parameters[4].valueAsText
        inclassval = parameters[5].valueAsText
        outclassval = parameters[6].valueAsText
        messages.addMessage(wbt.lidar_elevation_slice(input, output, minz, maxz, class1, inclassval, outclassval))
        return


class LidarGroundPointFilter(object):
    def __init__(self):
        self.label = "Lidar Ground Point Filter"
        self.description = "Identifies ground points within LiDAR dataset using a slope-based method."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        radius.value = "2.0"

        min_neighbours = arcpy.Parameter(
            displayName="Minimum Number of Neighbours",
            name="min_neighbours",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        min_neighbours.value = "0"

        slope_threshold = arcpy.Parameter(
            displayName="Inter-point Slope Threshold",
            name="slope_threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        slope_threshold.value = "45.0"

        height_threshold = arcpy.Parameter(
            displayName="Off-terrain Point Height Threshold",
            name="height_threshold",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        height_threshold.value = "1.0"

        classify = arcpy.Parameter(
            displayName="Classify Points",
            name="classify",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        classify.value = "true"

        slope_norm = arcpy.Parameter(
            displayName="Perform initial ground slope normalization?",
            name="slope_norm",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        slope_norm.value = "true"

        params = [input, output, radius, min_neighbours, slope_threshold, height_threshold, classify, slope_norm]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        radius = parameters[2].valueAsText
        min_neighbours = parameters[3].valueAsText
        slope_threshold = parameters[4].valueAsText
        height_threshold = parameters[5].valueAsText
        classify = parameters[6].valueAsText
        slope_norm = parameters[7].valueAsText
        messages.addMessage(wbt.lidar_ground_point_filter(input, output, radius, min_neighbours, slope_threshold, height_threshold, classify, slope_norm))
        return


class LidarHexBinning(object):
    def __init__(self):
        self.label = "Lidar Hex Binning"
        self.description = "Hex-bins a set of LiDAR points."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Base File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        width = arcpy.Parameter(
            displayName="Hexagon Width",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        orientation = arcpy.Parameter(
            displayName="Grid Orientation",
            name="orientation",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        orientation.filter.type = "ValueList"
        orientation.filter.list = ['horizontal', 'vertical']

        orientation.value = "horizontal"

        params = [input, output, width, orientation]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        width = parameters[2].valueAsText
        orientation = parameters[3].valueAsText
        messages.addMessage(wbt.lidar_hex_binning(input, output, width, orientation))
        return


class LidarHillshade(object):
    def __init__(self):
        self.label = "Lidar Hillshade"
        self.description = "Calculates a hillshade value for points within a LAS file and stores these data in the RGB field."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        azimuth = arcpy.Parameter(
            displayName="Azimuth (degrees)",
            name="azimuth",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        azimuth.value = "315.0"

        altitude = arcpy.Parameter(
            displayName="Altitude (degrees)",
            name="altitude",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        altitude.value = "30.0"

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        radius.value = "1.0"

        params = [input, output, azimuth, altitude, radius]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        azimuth = parameters[2].valueAsText
        altitude = parameters[3].valueAsText
        radius = parameters[4].valueAsText
        messages.addMessage(wbt.lidar_hillshade(input, output, azimuth, altitude, radius))
        return


class LidarHistogram(object):
    def __init__(self):
        self.label = "Lidar Histogram"
        self.description = "Creates a histogram of LiDAR data."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        parameter = arcpy.Parameter(
            displayName="Parameter",
            name="parameter",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        parameter.filter.type = "ValueList"
        parameter.filter.list = ['elevation', 'intensity', 'scan angle', 'class']

        parameter.value = "elevation"

        clip = arcpy.Parameter(
            displayName="Tail Clip Percent",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "1.0"

        params = [input, output, parameter, clip]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        parameter = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        messages.addMessage(wbt.lidar_histogram(input, output, parameter, clip))
        return


class LidarIdwInterpolation(object):
    def __init__(self):
        self.label = "Lidar Idw Interpolation"
        self.description = "Interpolates LAS files using an inverse-distance weighted (IDW) scheme. When the input/output parameters are not specified, the tool interpolates all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        parameter = arcpy.Parameter(
            displayName="Interpolation Parameter",
            name="parameter",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        parameter.filter.type = "ValueList"
        parameter.filter.list = ['elevation', 'intensity', 'class', 'scan angle', 'user data']

        parameter.value = "elevation"

        returns = arcpy.Parameter(
            displayName="Point Returns Included",
            name="returns",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        returns.filter.type = "ValueList"
        returns.filter.list = ['all', 'last', 'first']

        returns.value = "all"

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        weight = arcpy.Parameter(
            displayName="IDW Weight (Exponent) Value",
            name="weight",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        weight.value = "1.0"

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        radius.value = "2.5"

        exclude_cls = arcpy.Parameter(
            displayName="Exclusion Classes (0-18, based on LAS spec; e.g. 3,4,5,6,7)",
            name="exclude_cls",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        minz = arcpy.Parameter(
            displayName="Minimum Elevation Value (optional)",
            name="minz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxz = arcpy.Parameter(
            displayName="Maximum Elevation Value (optional)",
            name="maxz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, parameter, returns, resolution, weight, radius, exclude_cls, minz, maxz]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        parameter = parameters[2].valueAsText
        returns = parameters[3].valueAsText
        resolution = parameters[4].valueAsText
        weight = parameters[5].valueAsText
        radius = parameters[6].valueAsText
        exclude_cls = parameters[7].valueAsText
        minz = parameters[8].valueAsText
        maxz = parameters[9].valueAsText
        messages.addMessage(wbt.lidar_idw_interpolation(input, output, parameter, returns, resolution, weight, radius, exclude_cls, minz, maxz))
        return


class LidarInfo(object):
    def __init__(self):
        self.label = "Lidar Info"
        self.description = "Prints information about a LiDAR (LAS) dataset, including header, point return frequency, and classification data and information about the variable length records (VLRs) and geokeys."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output Summary Report File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        vlr = arcpy.Parameter(
            displayName="Print the variable length records (VLRs)?",
            name="vlr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        geokeys = arcpy.Parameter(
            displayName="Print the geokeys?",
            name="geokeys",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input, output, vlr, geokeys]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        vlr = parameters[2].valueAsText
        geokeys = parameters[3].valueAsText
        messages.addMessage(wbt.lidar_info(input, output, vlr, geokeys))
        return


class LidarJoin(object):
    def __init__(self):
        self.label = "Lidar Join"
        self.description = "Joins multiple LiDAR (LAS) files into a single LAS file."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input LiDAR Files",
            name="inputs",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        inputs.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.lidar_join(inputs, output))
        return


class LidarKappaIndex(object):
    def __init__(self):
        self.label = "Lidar Kappa Index"
        self.description = "Performs a kappa index of agreement (KIA) analysis on the classifications of two LAS files."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input LiDAR File (Classification)",
            name="input1",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input1.filter.list = ["las", "zip"]

        input2 = arcpy.Parameter(
            displayName="Input LiDAR File (Reference)",
            name="input2",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input2.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        class_accuracy = arcpy.Parameter(
            displayName="Output Class Accuracy Raster File",
            name="class_accuracy",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        class_accuracy.filter.list = ["tif"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        params = [input1, input2, output, class_accuracy, resolution]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        class_accuracy = parameters[3].valueAsText
        resolution = parameters[4].valueAsText
        messages.addMessage(wbt.lidar_kappa_index(input1, input2, output, class_accuracy, resolution))
        return


class LidarNearestNeighbourGridding(object):
    def __init__(self):
        self.label = "Lidar Nearest Neighbour Gridding"
        self.description = "Grids LAS files using nearest-neighbour scheme. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        parameter = arcpy.Parameter(
            displayName="Interpolation Parameter",
            name="parameter",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        parameter.filter.type = "ValueList"
        parameter.filter.list = ['elevation', 'intensity', 'class', 'scan angle', 'user data']

        parameter.value = "elevation"

        returns = arcpy.Parameter(
            displayName="Point Returns Included",
            name="returns",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        returns.filter.type = "ValueList"
        returns.filter.list = ['all', 'last', 'first']

        returns.value = "all"

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        radius.value = "2.5"

        exclude_cls = arcpy.Parameter(
            displayName="Exclusion Classes (0-18, based on LAS spec; e.g. 3,4,5,6,7)",
            name="exclude_cls",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        minz = arcpy.Parameter(
            displayName="Minimum Elevation Value (optional)",
            name="minz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxz = arcpy.Parameter(
            displayName="Maximum Elevation Value (optional)",
            name="maxz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, parameter, returns, resolution, radius, exclude_cls, minz, maxz]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        parameter = parameters[2].valueAsText
        returns = parameters[3].valueAsText
        resolution = parameters[4].valueAsText
        radius = parameters[5].valueAsText
        exclude_cls = parameters[6].valueAsText
        minz = parameters[7].valueAsText
        maxz = parameters[8].valueAsText
        messages.addMessage(wbt.lidar_nearest_neighbour_gridding(input, output, parameter, returns, resolution, radius, exclude_cls, minz, maxz))
        return


class LidarPointDensity(object):
    def __init__(self):
        self.label = "Lidar Point Density"
        self.description = "Calculates the spatial pattern of point density for a LiDAR data set. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        returns = arcpy.Parameter(
            displayName="Point Returns Included",
            name="returns",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        returns.filter.type = "ValueList"
        returns.filter.list = ['all', 'last', 'first']

        returns.value = "all"

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        radius.value = "2.5"

        exclude_cls = arcpy.Parameter(
            displayName="Exclusion Classes (0-18, based on LAS spec; e.g. 3,4,5,6,7)",
            name="exclude_cls",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        minz = arcpy.Parameter(
            displayName="Minimum Elevation Value (optional)",
            name="minz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxz = arcpy.Parameter(
            displayName="Maximum Elevation Value (optional)",
            name="maxz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, returns, resolution, radius, exclude_cls, minz, maxz]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        returns = parameters[2].valueAsText
        resolution = parameters[3].valueAsText
        radius = parameters[4].valueAsText
        exclude_cls = parameters[5].valueAsText
        minz = parameters[6].valueAsText
        maxz = parameters[7].valueAsText
        messages.addMessage(wbt.lidar_point_density(input, output, returns, resolution, radius, exclude_cls, minz, maxz))
        return


class LidarPointStats(object):
    def __init__(self):
        self.label = "Lidar Point Stats"
        self.description = "Creates several rasters summarizing the distribution of LAS point data. When the input/output parameters are not specified, the tool works on all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        num_points = arcpy.Parameter(
            displayName="Output number of points?",
            name="num_points",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        num_points.value = "True"

        num_pulses = arcpy.Parameter(
            displayName="Output number of pulses?",
            name="num_pulses",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        z_range = arcpy.Parameter(
            displayName="Output elevation range?",
            name="z_range",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        intensity_range = arcpy.Parameter(
            displayName="Output intensity range?",
            name="intensity_range",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        predom_class = arcpy.Parameter(
            displayName="Output predominant class?",
            name="predom_class",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input, resolution, num_points, num_pulses, z_range, intensity_range, predom_class]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        resolution = parameters[1].valueAsText
        num_points = parameters[2].valueAsText
        num_pulses = parameters[3].valueAsText
        z_range = parameters[4].valueAsText
        intensity_range = parameters[5].valueAsText
        predom_class = parameters[6].valueAsText
        messages.addMessage(wbt.lidar_point_stats(input, resolution, num_points, num_pulses, z_range, intensity_range, predom_class))
        return


class LidarRemoveDuplicates(object):
    def __init__(self):
        self.label = "Lidar Remove Duplicates"
        self.description = "Removes duplicate points from a LiDAR data set."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        include_z = arcpy.Parameter(
            displayName="Include z-values in point comparison?",
            name="include_z",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        include_z.value = "false"

        params = [input, output, include_z]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        include_z = parameters[2].valueAsText
        messages.addMessage(wbt.lidar_remove_duplicates(input, output, include_z))
        return


class LidarRemoveOutliers(object):
    def __init__(self):
        self.label = "Lidar Remove Outliers"
        self.description = "Removes outliers (high and low points) in a LiDAR point cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        radius.value = "2.0"

        elev_diff = arcpy.Parameter(
            displayName="Max. Elevation Difference",
            name="elev_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        elev_diff.value = "50.0"

        params = [input, output, radius, elev_diff]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        radius = parameters[2].valueAsText
        elev_diff = parameters[3].valueAsText
        messages.addMessage(wbt.lidar_remove_outliers(input, output, radius, elev_diff))
        return


class LidarSegmentation(object):
    def __init__(self):
        self.label = "Lidar Segmentation"
        self.description = "Segments a LiDAR point cloud based on normal vectors."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        radius.value = "5.0"

        norm_diff = arcpy.Parameter(
            displayName="Normal Difference Threshold",
            name="norm_diff",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        norm_diff.value = "10.0"

        maxzdiff = arcpy.Parameter(
            displayName="Maximum Elevation Difference Between Points",
            name="maxzdiff",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        maxzdiff.value = "1.0"

        params = [input, output, radius, norm_diff, maxzdiff]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        radius = parameters[2].valueAsText
        norm_diff = parameters[3].valueAsText
        maxzdiff = parameters[4].valueAsText
        messages.addMessage(wbt.lidar_segmentation(input, output, radius, norm_diff, maxzdiff))
        return


class LidarSegmentationBasedFilter(object):
    def __init__(self):
        self.label = "Lidar Segmentation Based Filter"
        self.description = "Identifies ground points within LiDAR point clouds using a segmentation based approach."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        radius.value = "5.0"

        norm_diff = arcpy.Parameter(
            displayName="Normal Difference Threshold",
            name="norm_diff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        norm_diff.value = "2.0"

        maxzdiff = arcpy.Parameter(
            displayName="Maximum Elevation Difference Between Points",
            name="maxzdiff",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxzdiff.value = "1.0"

        classify = arcpy.Parameter(
            displayName="Classify Points",
            name="classify",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input, output, radius, norm_diff, maxzdiff, classify]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        radius = parameters[2].valueAsText
        norm_diff = parameters[3].valueAsText
        maxzdiff = parameters[4].valueAsText
        classify = parameters[5].valueAsText
        messages.addMessage(wbt.lidar_segmentation_based_filter(input, output, radius, norm_diff, maxzdiff, classify))
        return


class LidarThin(object):
    def __init__(self):
        self.label = "Lidar Thin"
        self.description = "Thins a LiDAR point cloud, reducing point density."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        resolution = arcpy.Parameter(
            displayName="Sample Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "2.0"

        method = arcpy.Parameter(
            displayName="Point Selection Method",
            name="method",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        method.filter.type = "ValueList"
        method.filter.list = ['first', 'last', 'lowest', 'highest', 'nearest']

        method.value = "lowest"

        save_filtered = arcpy.Parameter(
            displayName="Save filtered points to seperate file?",
            name="save_filtered",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        save_filtered.value = "false"

        params = [input, output, resolution, method, save_filtered]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        resolution = parameters[2].valueAsText
        method = parameters[3].valueAsText
        save_filtered = parameters[4].valueAsText
        messages.addMessage(wbt.lidar_thin(input, output, resolution, method, save_filtered))
        return


class LidarThinHighDensity(object):
    def __init__(self):
        self.label = "Lidar Thin High Density"
        self.description = "Thins points from high density areas within a LiDAR point cloud."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        density = arcpy.Parameter(
            displayName="Max. Point Density (pts/m^2)",
            name="density",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        save_filtered = arcpy.Parameter(
            displayName="Save filtered points to seperate file?",
            name="save_filtered",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        save_filtered.value = "false"

        params = [input, output, resolution, density, save_filtered]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        resolution = parameters[2].valueAsText
        density = parameters[3].valueAsText
        save_filtered = parameters[4].valueAsText
        messages.addMessage(wbt.lidar_thin_high_density(input, output, resolution, density, save_filtered))
        return


class LidarTile(object):
    def __init__(self):
        self.label = "Lidar Tile"
        self.description = "Tiles a LiDAR LAS file into multiple LAS files."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        width = arcpy.Parameter(
            displayName="Tile Width",
            name="width",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        width.value = "1000.0"

        height = arcpy.Parameter(
            displayName="Tile Height",
            name="height",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        height.value = "1000.0"

        origin_x = arcpy.Parameter(
            displayName="Origin Point X-Coordinate",
            name="origin_x",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        origin_x.value = "0.0"

        origin_y = arcpy.Parameter(
            displayName="Origin Point Y-Coordinate",
            name="origin_y",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        origin_y.value = "0.0"

        min_points = arcpy.Parameter(
            displayName="Minimum Number of Tile Points",
            name="min_points",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        min_points.value = "2"

        params = [input, width, height, origin_x, origin_y, min_points]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        width = parameters[1].valueAsText
        height = parameters[2].valueAsText
        origin_x = parameters[3].valueAsText
        origin_y = parameters[4].valueAsText
        min_points = parameters[5].valueAsText
        messages.addMessage(wbt.lidar_tile(input, width, height, origin_x, origin_y, min_points))
        return


class LidarTileFootprint(object):
    def __init__(self):
        self.label = "Lidar Tile Footprint"
        self.description = "Creates a vector polygon of the convex hull of a LiDAR point cloud. When the input/output parameters are not specified, the tool works with all LAS files contained within the working directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input LiDAR File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output Polygon File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polygon"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.lidar_tile_footprint(input, output))
        return


class LidarTinGridding(object):
    def __init__(self):
        self.label = "Lidar Tin Gridding"
        self.description = "Creates a raster grid based on a Delaunay triangular irregular network (TIN) fitted to LiDAR points."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Optional",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        parameter = arcpy.Parameter(
            displayName="Interpolation Parameter",
            name="parameter",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        parameter.filter.type = "ValueList"
        parameter.filter.list = ['elevation', 'intensity', 'class', 'scan angle', 'user data']

        parameter.value = "elevation"

        returns = arcpy.Parameter(
            displayName="Point Returns Included",
            name="returns",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        returns.filter.type = "ValueList"
        returns.filter.list = ['all', 'last', 'first']

        returns.value = "all"

        resolution = arcpy.Parameter(
            displayName="Grid Resolution",
            name="resolution",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        resolution.value = "1.0"

        exclude_cls = arcpy.Parameter(
            displayName="Exclusion Classes (0-18, based on LAS spec; e.g. 3,4,5,6,7)",
            name="exclude_cls",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        minz = arcpy.Parameter(
            displayName="Minimum Elevation Value (optional)",
            name="minz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxz = arcpy.Parameter(
            displayName="Maximum Elevation Value (optional)",
            name="maxz",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, parameter, returns, resolution, exclude_cls, minz, maxz]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        parameter = parameters[2].valueAsText
        returns = parameters[3].valueAsText
        resolution = parameters[4].valueAsText
        exclude_cls = parameters[5].valueAsText
        minz = parameters[6].valueAsText
        maxz = parameters[7].valueAsText
        messages.addMessage(wbt.lidar_tin_gridding(input, output, parameter, returns, resolution, exclude_cls, minz, maxz))
        return


class LidarTophatTransform(object):
    def __init__(self):
        self.label = "Lidar Tophat Transform"
        self.description = "Performs a white top-hat transform on a Lidar dataset; as an estimate of height above ground, this is useful for modelling the vegetation canopy."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        radius.value = "1.0"

        params = [input, output, radius]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        radius = parameters[2].valueAsText
        messages.addMessage(wbt.lidar_tophat_transform(input, output, radius))
        return


class NormalVectors(object):
    def __init__(self):
        self.label = "Normal Vectors"
        self.description = "Calculates normal vectors for points within a LAS file and stores these data (XYZ vector components) in the RGB field."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["las", "zip"]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["las", "zip"]

        radius = arcpy.Parameter(
            displayName="Search Radius",
            name="radius",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        radius.value = "1.0"

        params = [input, output, radius]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        radius = parameters[2].valueAsText
        messages.addMessage(wbt.normal_vectors(input, output, radius))
        return


class SelectTilesByPolygon(object):
    def __init__(self):
        self.label = "Select Tiles By Polygon"
        self.description = "Copies LiDAR tiles overlapping with a polygon into an output directory."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        indir = arcpy.Parameter(
            displayName="Input Directory",
            name="indir",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        outdir = arcpy.Parameter(
            displayName="Output Directory",
            name="outdir",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        polygons = arcpy.Parameter(
            displayName="Input Vector Polygon File",
            name="polygons",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        polygons.filter.list = ["Polygon"]

        params = [indir, outdir, polygons]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        indir = parameters[0].valueAsText
        outdir = parameters[1].valueAsText
        polygons = parameters[2].valueAsText
        messages.addMessage(wbt.select_tiles_by_polygon(indir, outdir, polygons))
        return


class And(object):
    def __init__(self):
        self.label = "And"
        self.description = "Performs a logical AND operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.And(input1, input2, output))
        return


class Not(object):
    def __init__(self):
        self.label = "Not"
        self.description = "Performs a logical NOT operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.Not(input1, input2, output))
        return


class Or(object):
    def __init__(self):
        self.label = "Or"
        self.description = "Performs a logical OR operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.Or(input1, input2, output))
        return


class AbsoluteValue(object):
    def __init__(self):
        self.label = "Absolute Value"
        self.description = "Calculates the absolute value of every cell in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.absolute_value(input, output))
        return


class Add(object):
    def __init__(self):
        self.label = "Add"
        self.description = "Performs an addition operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.add(input1, input2, output))
        return


class Anova(object):
    def __init__(self):
        self.label = "Anova"
        self.description = "Performs an analysis of variance (ANOVA) test on a raster dataset."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        features = arcpy.Parameter(
            displayName="Feature Definition (Class) File",
            name="features",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input, features, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        features = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.anova(input, features, output))
        return


class ArcCos(object):
    def __init__(self):
        self.label = "Arc Cos"
        self.description = "Returns the inverse cosine (arccos) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.arc_cos(input, output))
        return


class ArcSin(object):
    def __init__(self):
        self.label = "Arc Sin"
        self.description = "Returns the inverse sine (arcsin) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.arc_sin(input, output))
        return


class ArcTan(object):
    def __init__(self):
        self.label = "Arc Tan"
        self.description = "Returns the inverse tangent (arctan) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.arc_tan(input, output))
        return


class Atan2(object):
    def __init__(self):
        self.label = "Atan2"
        self.description = "Returns the 2-argument inverse tangent (atan2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input_y = arcpy.Parameter(
            displayName="Input Y File Or Constant Value",
            name="input_y",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input_x = arcpy.Parameter(
            displayName="Input X File Or Constant Value",
            name="input_x",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input_y, input_x, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input_y = parameters[0].valueAsText
        input_x = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.atan2(input_y, input_x, output))
        return


class AttributeCorrelation(object):
    def __init__(self):
        self.label = "Attribute Correlation"
        self.description = "Performs a correlation analysis on attribute fields from a vector database."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.attribute_correlation(input, output))
        return


class AttributeHistogram(object):
    def __init__(self):
        self.label = "Attribute Histogram"
        self.description = "Creates a histogram for the field values of a vector's attribute table."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input, field, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.attribute_histogram(input, field, output))
        return


class AttributeScattergram(object):
    def __init__(self):
        self.label = "Attribute Scattergram"
        self.description = "Creates a scattergram for two field values of a vector's attribute table."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        fieldx = arcpy.Parameter(
            displayName="Field Name X",
            name="fieldx",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        fieldx.parameterDependencies = [input.name]

        fieldy = arcpy.Parameter(
            displayName="Field Name Y",
            name="fieldy",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        fieldy.parameterDependencies = [input.name]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        trendline = arcpy.Parameter(
            displayName="Draw the trendline?",
            name="trendline",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        trendline.value = "false"

        params = [input, fieldx, fieldy, output, trendline]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        fieldx = parameters[1].valueAsText
        fieldy = parameters[2].valueAsText
        output = parameters[3].valueAsText
        trendline = parameters[4].valueAsText
        messages.addMessage(wbt.attribute_scattergram(input, fieldx, fieldy, output, trendline))
        return


class Ceil(object):
    def __init__(self):
        self.label = "Ceil"
        self.description = "Returns the smallest (closest to negative infinity) value that is greater than or equal to the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.ceil(input, output))
        return


class Cos(object):
    def __init__(self):
        self.label = "Cos"
        self.description = "Returns the cosine (cos) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.cos(input, output))
        return


class Cosh(object):
    def __init__(self):
        self.label = "Cosh"
        self.description = "Returns the hyperbolic cosine (cosh) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.cosh(input, output))
        return


class CrispnessIndex(object):
    def __init__(self):
        self.label = "Crispness Index"
        self.description = "Calculates the Crispness Index, which is used to quantify how crisp (or conversely how fuzzy) a probability image is."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.crispness_index(input, output))
        return


class CrossTabulation(object):
    def __init__(self):
        self.label = "Cross Tabulation"
        self.description = "Performs a cross-tabulation on two categorical images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File 1",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File 2",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.cross_tabulation(input1, input2, output))
        return


class CumulativeDistribution(object):
    def __init__(self):
        self.label = "Cumulative Distribution"
        self.description = "Converts a raster image to its cumulative distribution function."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.cumulative_distribution(input, output))
        return


class Decrement(object):
    def __init__(self):
        self.label = "Decrement"
        self.description = "Decreases the values of each grid cell in an input raster by 1.0 (see also InPlaceSubtract)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.decrement(input, output))
        return


class Divide(object):
    def __init__(self):
        self.label = "Divide"
        self.description = "Performs a division operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.divide(input1, input2, output))
        return


class EqualTo(object):
    def __init__(self):
        self.label = "Equal To"
        self.description = "Performs a equal-to comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.equal_to(input1, input2, output))
        return


class Exp(object):
    def __init__(self):
        self.label = "Exp"
        self.description = "Returns the exponential (base e) of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.exp(input, output))
        return


class Exp2(object):
    def __init__(self):
        self.label = "Exp2"
        self.description = "Returns the exponential (base 2) of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.exp2(input, output))
        return


class ExtractRasterStatistics(object):
    def __init__(self):
        self.label = "Extract Raster Statistics"
        self.description = "Extracts descriptive statistics for a group of patches in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Data File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        features = arcpy.Parameter(
            displayName="Input Feature Definition File",
            name="features",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Raster File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        stat = arcpy.Parameter(
            displayName="Statistic Type",
            name="stat",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        stat.filter.type = "ValueList"
        stat.filter.list = ['average', 'minimum', 'maximum', 'range', 'standard deviation', 'total']

        stat.value = "average"

        out_table = arcpy.Parameter(
            displayName="Output HTML Table File",
            name="out_table",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_table.filter.list = ["html"]

        params = [input, features, output, stat, out_table]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        features = parameters[1].valueAsText
        output = parameters[2].valueAsText
        stat = parameters[3].valueAsText
        out_table = parameters[4].valueAsText
        messages.addMessage(wbt.extract_raster_statistics(input, features, output, stat, out_table))
        return


class Floor(object):
    def __init__(self):
        self.label = "Floor"
        self.description = "Returns the largest (closest to positive infinity) value that is less than or equal to the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.floor(input, output))
        return


class GreaterThan(object):
    def __init__(self):
        self.label = "Greater Than"
        self.description = "Performs a greater-than comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        incl_equals = arcpy.Parameter(
            displayName="Perform a greater-than-OR-EQUAL-TO operation?",
            name="incl_equals",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input1, input2, output, incl_equals]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        incl_equals = parameters[3].valueAsText
        messages.addMessage(wbt.greater_than(input1, input2, output, incl_equals))
        return


class ImageAutocorrelation(object):
    def __init__(self):
        self.label = "Image Autocorrelation"
        self.description = "Performs Moran's I analysis on two or more input images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        contiguity = arcpy.Parameter(
            displayName="Contiguity Type",
            name="contiguity",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        contiguity.filter.type = "ValueList"
        contiguity.filter.list = ['Rook', 'King', 'Bishop']

        contiguity.value = "Rook"

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [inputs, contiguity, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        contiguity = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.image_autocorrelation(inputs, contiguity, output))
        return


class ImageCorrelation(object):
    def __init__(self):
        self.label = "Image Correlation"
        self.description = "Performs image correlation on two or more input images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [inputs, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.image_correlation(inputs, output))
        return


class ImageRegression(object):
    def __init__(self):
        self.label = "Image Regression"
        self.description = "Performs image regression analysis on two input images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Independent Variable (X).",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Dependent Variable (Y).",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output Summary Report File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        out_residuals = arcpy.Parameter(
            displayName="Optional Residuals Output File",
            name="out_residuals",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        out_residuals.filter.list = ["tif"]

        standardize = arcpy.Parameter(
            displayName="Standardize the residuals map?",
            name="standardize",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input1, input2, output, out_residuals, standardize]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        out_residuals = parameters[3].valueAsText
        standardize = parameters[4].valueAsText
        messages.addMessage(wbt.image_regression(input1, input2, output, out_residuals, standardize))
        return


class InPlaceAdd(object):
    def __init__(self):
        self.label = "In Place Add"
        self.description = "Performs an in-place addition operation (input1 += input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Raster File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input1, input2]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        messages.addMessage(wbt.in_place_add(input1, input2))
        return


class InPlaceDivide(object):
    def __init__(self):
        self.label = "In Place Divide"
        self.description = "Performs an in-place division operation (input1 /= input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Raster File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input1, input2]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        messages.addMessage(wbt.in_place_divide(input1, input2))
        return


class InPlaceMultiply(object):
    def __init__(self):
        self.label = "In Place Multiply"
        self.description = "Performs an in-place multiplication operation (input1 *= input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Raster File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input1, input2]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        messages.addMessage(wbt.in_place_multiply(input1, input2))
        return


class InPlaceSubtract(object):
    def __init__(self):
        self.label = "In Place Subtract"
        self.description = "Performs an in-place subtraction operation (input1 -= input2)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Raster File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input1, input2]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        messages.addMessage(wbt.in_place_subtract(input1, input2))
        return


class Increment(object):
    def __init__(self):
        self.label = "Increment"
        self.description = "Increases the values of each grid cell in an input raster by 1.0. (see also InPlaceAdd)."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.increment(input, output))
        return


class IntegerDivision(object):
    def __init__(self):
        self.label = "Integer Division"
        self.description = "Performs an integer division operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.integer_division(input1, input2, output))
        return


class IsNoData(object):
    def __init__(self):
        self.label = "Is No Data"
        self.description = "Identifies NoData valued pixels in an image."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.is_no_data(input, output))
        return


class KappaIndex(object):
    def __init__(self):
        self.label = "Kappa Index"
        self.description = "Performs a kappa index of agreement (KIA) analysis on two categorical raster files."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input Classification File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input Reference File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.kappa_index(input1, input2, output))
        return


class KsTestForNormality(object):
    def __init__(self):
        self.label = "Ks Test For Normality"
        self.description = "Evaluates whether the values in a raster are normally distributed."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        num_samples = arcpy.Parameter(
            displayName="Num. Samples (blank for while image)",
            name="num_samples",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        params = [input, output, num_samples]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_samples = parameters[2].valueAsText
        messages.addMessage(wbt.ks_test_for_normality(input, output, num_samples))
        return


class LessThan(object):
    def __init__(self):
        self.label = "Less Than"
        self.description = "Performs a less-than comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        incl_equals = arcpy.Parameter(
            displayName="Perform a less-than-OR-EQUAL-TO operation?",
            name="incl_equals",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [input1, input2, output, incl_equals]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        incl_equals = parameters[3].valueAsText
        messages.addMessage(wbt.less_than(input1, input2, output, incl_equals))
        return


class ListUniqueValues(object):
    def __init__(self):
        self.label = "List Unique Values"
        self.description = "Lists the unique values contained in a field witin a vector's attribute table."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input, field, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.list_unique_values(input, field, output))
        return


class Ln(object):
    def __init__(self):
        self.label = "Ln"
        self.description = "Returns the natural logarithm of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.ln(input, output))
        return


class Log10(object):
    def __init__(self):
        self.label = "Log10"
        self.description = "Returns the base-10 logarithm of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.log10(input, output))
        return


class Log2(object):
    def __init__(self):
        self.label = "Log2"
        self.description = "Returns the base-2 logarithm of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.log2(input, output))
        return


class Max(object):
    def __init__(self):
        self.label = "Max"
        self.description = "Performs a MAX operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.max(input1, input2, output))
        return


class Min(object):
    def __init__(self):
        self.label = "Min"
        self.description = "Performs a MIN operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.min(input1, input2, output))
        return


class Modulo(object):
    def __init__(self):
        self.label = "Modulo"
        self.description = "Performs a modulo operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.modulo(input1, input2, output))
        return


class Multiply(object):
    def __init__(self):
        self.label = "Multiply"
        self.description = "Performs a multiplication operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.multiply(input1, input2, output))
        return


class Negate(object):
    def __init__(self):
        self.label = "Negate"
        self.description = "Changes the sign of values in a raster or the 0-1 values of a Boolean raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.negate(input, output))
        return


class NotEqualTo(object):
    def __init__(self):
        self.label = "Not Equal To"
        self.description = "Performs a not-equal-to comparison operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.not_equal_to(input1, input2, output))
        return


class Power(object):
    def __init__(self):
        self.label = "Power"
        self.description = "Raises the values in grid cells of one rasters, or a constant value, by values in another raster or constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.power(input1, input2, output))
        return


class PrincipalComponentAnalysis(object):
    def __init__(self):
        self.label = "Principal Component Analysis"
        self.description = "Performs a principal component analysis (PCA) on a multi-spectral dataset."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML Report File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        num_comp = arcpy.Parameter(
            displayName="Num. of Component Images (blank for all)",
            name="num_comp",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        standardized = arcpy.Parameter(
            displayName="Perform Standaradized PCA?",
            name="standardized",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [inputs, output, num_comp, standardized]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_comp = parameters[2].valueAsText
        standardized = parameters[3].valueAsText
        messages.addMessage(wbt.principal_component_analysis(inputs, output, num_comp, standardized))
        return


class Quantiles(object):
    def __init__(self):
        self.label = "Quantiles"
        self.description = "Transforms raster values into quantiles."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        num_quantiles = arcpy.Parameter(
            displayName="Number of Quantiles",
            name="num_quantiles",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_quantiles.value = "5"

        params = [input, output, num_quantiles]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_quantiles = parameters[2].valueAsText
        messages.addMessage(wbt.quantiles(input, output, num_quantiles))
        return


class RandomField(object):
    def __init__(self):
        self.label = "Random Field"
        self.description = "Creates an image containing random values."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        base = arcpy.Parameter(
            displayName="Input Base File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [base, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        base = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.random_field(base, output))
        return


class RandomSample(object):
    def __init__(self):
        self.label = "Random Sample"
        self.description = "Creates an image containing randomly located sample grid cells with unique IDs."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        base = arcpy.Parameter(
            displayName="Input Base File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        num_samples = arcpy.Parameter(
            displayName="Num. Samples",
            name="num_samples",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        num_samples.value = "1000"

        params = [base, output, num_samples]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        base = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_samples = parameters[2].valueAsText
        messages.addMessage(wbt.random_sample(base, output, num_samples))
        return


class RasterHistogram(object):
    def __init__(self):
        self.label = "Raster Histogram"
        self.description = "Creates a histogram from raster values."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.raster_histogram(input, output))
        return


class RasterSummaryStats(object):
    def __init__(self):
        self.label = "Raster Summary Stats"
        self.description = "Measures a rasters min, max, average, standard deviation, num. non-nodata cells, and total."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        messages.addMessage(wbt.raster_summary_stats(input))
        return


class Reciprocal(object):
    def __init__(self):
        self.label = "Reciprocal"
        self.description = "Returns the reciprocal (i.e. 1 / z) of values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.reciprocal(input, output))
        return


class RescaleValueRange(object):
    def __init__(self):
        self.label = "Rescale Value Range"
        self.description = "Performs a min-max contrast stretch on an input greytone image."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        out_min_val = arcpy.Parameter(
            displayName="Output Raster Minimum Value",
            name="out_min_val",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        out_max_val = arcpy.Parameter(
            displayName="Output Raster Maximum Value",
            name="out_max_val",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        clip_min = arcpy.Parameter(
            displayName="Lower-Tail Clip Value (optional)",
            name="clip_min",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip_max = arcpy.Parameter(
            displayName="Upper-Tail Clip Value (optional)",
            name="clip_max",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, out_min_val, out_max_val, clip_min, clip_max]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        out_min_val = parameters[2].valueAsText
        out_max_val = parameters[3].valueAsText
        clip_min = parameters[4].valueAsText
        clip_max = parameters[5].valueAsText
        messages.addMessage(wbt.rescale_value_range(input, output, out_min_val, out_max_val, clip_min, clip_max))
        return


class RootMeanSquareError(object):
    def __init__(self):
        self.label = "Root Mean Square Error"
        self.description = "Calculates the RMSE and other accuracy statistics."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Input Base File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        params = [input, base]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        base = parameters[1].valueAsText
        messages.addMessage(wbt.root_mean_square_error(input, base))
        return


class Round(object):
    def __init__(self):
        self.label = "Round"
        self.description = "Rounds the values in an input raster to the nearest integer value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.round(input, output))
        return


class Sin(object):
    def __init__(self):
        self.label = "Sin"
        self.description = "Returns the sine (sin) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.sin(input, output))
        return


class Sinh(object):
    def __init__(self):
        self.label = "Sinh"
        self.description = "Returns the hyperbolic sine (sinh) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.sinh(input, output))
        return


class Square(object):
    def __init__(self):
        self.label = "Square"
        self.description = "Squares the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.square(input, output))
        return


class SquareRoot(object):
    def __init__(self):
        self.label = "Square Root"
        self.description = "Returns the square root of the values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.square_root(input, output))
        return


class Subtract(object):
    def __init__(self):
        self.label = "Subtract"
        self.description = "Performs a differencing operation on two rasters or a raster and a constant value."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File Or Constant Value",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.subtract(input1, input2, output))
        return


class Tan(object):
    def __init__(self):
        self.label = "Tan"
        self.description = "Returns the tangent (tan) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.tan(input, output))
        return


class Tanh(object):
    def __init__(self):
        self.label = "Tanh"
        self.description = "Returns the hyperbolic tangent (tanh) of each values in a raster."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.tanh(input, output))
        return


class ToDegrees(object):
    def __init__(self):
        self.label = "To Degrees"
        self.description = "Converts a raster from radians to degrees."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.to_degrees(input, output))
        return


class ToRadians(object):
    def __init__(self):
        self.label = "To Radians"
        self.description = "Converts a raster from degrees to radians."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.to_radians(input, output))
        return


class TrendSurface(object):
    def __init__(self):
        self.label = "Trend Surface"
        self.description = "Estimates the trend surface of an input raster file."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        order = arcpy.Parameter(
            displayName="Polynomial Order",
            name="order",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        order.value = "1"

        params = [input, output, order]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        order = parameters[2].valueAsText
        messages.addMessage(wbt.trend_surface(input, output, order))
        return


class TrendSurfaceVectorPoints(object):
    def __init__(self):
        self.label = "Trend Surface Vector Points"
        self.description = "Estimates a trend surface from vector points."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        input.filter.list = ["Point"]

        field = arcpy.Parameter(
            displayName="Field Name",
            name="field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        field.parameterDependencies = [input.name]

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        order = arcpy.Parameter(
            displayName="Polynomial Order",
            name="order",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        order.value = "1"

        cell_size = arcpy.Parameter(
            displayName="Cell Size (optional)",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, field, output, order, cell_size]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        field = parameters[1].valueAsText
        output = parameters[2].valueAsText
        order = parameters[3].valueAsText
        cell_size = parameters[4].valueAsText
        messages.addMessage(wbt.trend_surface_vector_points(input, field, output, order, cell_size))
        return


class Truncate(object):
    def __init__(self):
        self.label = "Truncate"
        self.description = "Truncates the values in a raster to the desired number of decimal places."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        num_decimals = arcpy.Parameter(
            displayName="Number of Decimals After Truncation",
            name="num_decimals",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        params = [input, output, num_decimals]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        num_decimals = parameters[2].valueAsText
        messages.addMessage(wbt.truncate(input, output, num_decimals))
        return


class TurningBandsSimulation(object):
    def __init__(self):
        self.label = "Turning Bands Simulation"
        self.description = "Creates an image containing random values based on a turning-bands simulation."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        base = arcpy.Parameter(
            displayName="Input Base File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        range = arcpy.Parameter(
            displayName="Range of Autocorrelation (map units)",
            name="range",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        iterations = arcpy.Parameter(
            displayName="Iterations",
            name="iterations",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        iterations.value = "1000"

        params = [base, output, range, iterations]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        base = parameters[0].valueAsText
        output = parameters[1].valueAsText
        range = parameters[2].valueAsText
        iterations = parameters[3].valueAsText
        messages.addMessage(wbt.turning_bands_simulation(base, output, range, iterations))
        return


class Xor(object):
    def __init__(self):
        self.label = "Xor"
        self.description = "Performs a logical XOR operator on two Boolean raster images."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input File",
            name="input2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input1, input2, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        messages.addMessage(wbt.xor(input1, input2, output))
        return


class ZScores(object):
    def __init__(self):
        self.label = "Z Scores"
        self.description = "Standardizes the values in an input raster by converting to z-scores."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        params = [input, output]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        messages.addMessage(wbt.z_scores(input, output))
        return


class DistanceToOutlet(object):
    def __init__(self):
        self.label = "Distance To Outlet"
        self.description = "Calculates the distance of stream grid cells to the channel network outlet cell."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.distance_to_outlet(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class ExtractStreams(object):
    def __init__(self):
        self.label = "Extract Streams"
        self.description = "Extracts stream grid cells from a flow accumulation raster."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        flow_accum = arcpy.Parameter(
            displayName="Input D8 Flow Accumulation File",
            name="flow_accum",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        threshold = arcpy.Parameter(
            displayName="Channelization Threshold",
            name="threshold",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [flow_accum, output, threshold, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        flow_accum = parameters[0].valueAsText
        output = parameters[1].valueAsText
        threshold = parameters[2].valueAsText
        zero_background = parameters[3].valueAsText
        messages.addMessage(wbt.extract_streams(flow_accum, output, threshold, zero_background))
        return


class ExtractValleys(object):
    def __init__(self):
        self.label = "Extract Valleys"
        self.description = "Identifies potential valley bottom grid cells based on local topolography alone."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        variant = arcpy.Parameter(
            displayName="Variant",
            name="variant",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        variant.filter.type = "ValueList"
        variant.filter.list = ['Lower Quartile', 'Johnston and Rosenfeld', 'Peucker and Douglas']

        variant.value = "Lower Quartile"

        line_thin = arcpy.Parameter(
            displayName="Perform line-thinning?",
            name="line_thin",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        line_thin.value = "true"

        filter = arcpy.Parameter(
            displayName="Filter Size (Only For Lower Quartile)",
            name="filter",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")

        filter.value = "5"

        params = [dem, output, variant, line_thin, filter]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        dem = parameters[0].valueAsText
        output = parameters[1].valueAsText
        variant = parameters[2].valueAsText
        line_thin = parameters[3].valueAsText
        filter = parameters[4].valueAsText
        messages.addMessage(wbt.extract_valleys(dem, output, variant, line_thin, filter))
        return


class FarthestChannelHead(object):
    def __init__(self):
        self.label = "Farthest Channel Head"
        self.description = "Calculates the distance to the furthest upstream channel head for each stream cell."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.farthest_channel_head(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class FindMainStem(object):
    def __init__(self):
        self.label = "Find Main Stem"
        self.description = "Finds the main stem, based on stream lengths, of each stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.find_main_stem(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class HackStreamOrder(object):
    def __init__(self):
        self.label = "Hack Stream Order"
        self.description = "Assigns the Hack stream order to each tributary in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.hack_stream_order(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class HortonStreamOrder(object):
    def __init__(self):
        self.label = "Horton Stream Order"
        self.description = "Assigns the Horton stream order to each tributary in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.horton_stream_order(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class LengthOfUpstreamChannels(object):
    def __init__(self):
        self.label = "Length Of Upstream Channels"
        self.description = "Calculates the total length of channels upstream."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.length_of_upstream_channels(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class LongProfile(object):
    def __init__(self):
        self.label = "Long Profile"
        self.description = "Plots the stream longitudinal profiles for one or more rivers."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, streams, dem, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        dem = parameters[2].valueAsText
        output = parameters[3].valueAsText
        esri_pntr = parameters[4].valueAsText
        messages.addMessage(wbt.long_profile(d8_pntr, streams, dem, output, esri_pntr))
        return


class LongProfileFromPoints(object):
    def __init__(self):
        self.label = "Long Profile From Points"
        self.description = "Plots the longitudinal profiles from flow-paths initiating from a set of vector points."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        points = arcpy.Parameter(
            displayName="Input Vector Points File",
            name="points",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        points.filter.list = ["Point"]

        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output HTML File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["html"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, points, dem, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        points = parameters[1].valueAsText
        dem = parameters[2].valueAsText
        output = parameters[3].valueAsText
        esri_pntr = parameters[4].valueAsText
        messages.addMessage(wbt.long_profile_from_points(d8_pntr, points, dem, output, esri_pntr))
        return


class RasterStreamsToVector(object):
    def __init__(self):
        self.label = "Raster Streams To Vector"
        self.description = "Converts a raster stream file into a vector file."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["Polyline"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [streams, d8_pntr, output, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        streams = parameters[0].valueAsText
        d8_pntr = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        messages.addMessage(wbt.raster_streams_to_vector(streams, d8_pntr, output, esri_pntr))
        return


class RasterizeStreams(object):
    def __init__(self):
        self.label = "Rasterize Streams"
        self.description = "Rasterizes vector streams based on Lindsay (2016) method."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        streams = arcpy.Parameter(
            displayName="Input Vector Streams File",
            name="streams",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        streams.filter.list = ["Polyline"]

        base = arcpy.Parameter(
            displayName="Input Base Raster File",
            name="base",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        nodata = arcpy.Parameter(
            displayName="Use NoData value for background?",
            name="nodata",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        nodata.value = "true"

        feature_id = arcpy.Parameter(
            displayName="Use feature number as output value?",
            name="feature_id",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        feature_id.value = "false"

        params = [streams, base, output, nodata, feature_id]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        streams = parameters[0].valueAsText
        base = parameters[1].valueAsText
        output = parameters[2].valueAsText
        nodata = parameters[3].valueAsText
        feature_id = parameters[4].valueAsText
        messages.addMessage(wbt.rasterize_streams(streams, base, output, nodata, feature_id))
        return


class RemoveShortStreams(object):
    def __init__(self):
        self.label = "Remove Short Streams"
        self.description = "Removes short first-order streams from a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        min_length = arcpy.Parameter(
            displayName="Minimum Tributary Length (map units)",
            name="min_length",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        params = [d8_pntr, streams, output, min_length, esri_pntr]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        min_length = parameters[3].valueAsText
        esri_pntr = parameters[4].valueAsText
        messages.addMessage(wbt.remove_short_streams(d8_pntr, streams, output, min_length, esri_pntr))
        return


class ShreveStreamMagnitude(object):
    def __init__(self):
        self.label = "Shreve Stream Magnitude"
        self.description = "Assigns the Shreve stream magnitude to each link in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.shreve_stream_magnitude(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class StrahlerStreamOrder(object):
    def __init__(self):
        self.label = "Strahler Stream Order"
        self.description = "Assigns the Strahler stream order to each link in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.strahler_stream_order(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class StreamLinkClass(object):
    def __init__(self):
        self.label = "Stream Link Class"
        self.description = "Identifies the exterior/interior links and nodes in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.stream_link_class(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class StreamLinkIdentifier(object):
    def __init__(self):
        self.label = "Stream Link Identifier"
        self.description = "Assigns a unique identifier to each link in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.stream_link_identifier(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class StreamLinkLength(object):
    def __init__(self):
        self.label = "Stream Link Length"
        self.description = "Estimates the length of each link (or tributary) in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        linkid = arcpy.Parameter(
            displayName="Input Stream Link (Tributary) ID File",
            name="linkid",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, linkid, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        linkid = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.stream_link_length(d8_pntr, linkid, output, esri_pntr, zero_background))
        return


class StreamLinkSlope(object):
    def __init__(self):
        self.label = "Stream Link Slope"
        self.description = "Estimates the average slope of each link (or tributary) in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        linkid = arcpy.Parameter(
            displayName="Input Stream Link (Tributary) ID File",
            name="linkid",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, linkid, dem, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        linkid = parameters[1].valueAsText
        dem = parameters[2].valueAsText
        output = parameters[3].valueAsText
        esri_pntr = parameters[4].valueAsText
        zero_background = parameters[5].valueAsText
        messages.addMessage(wbt.stream_link_slope(d8_pntr, linkid, dem, output, esri_pntr, zero_background))
        return


class StreamSlopeContinuous(object):
    def __init__(self):
        self.label = "Stream Slope Continuous"
        self.description = "Estimates the slope of each grid cell in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, dem, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        dem = parameters[2].valueAsText
        output = parameters[3].valueAsText
        esri_pntr = parameters[4].valueAsText
        zero_background = parameters[5].valueAsText
        messages.addMessage(wbt.stream_slope_continuous(d8_pntr, streams, dem, output, esri_pntr, zero_background))
        return


class TopologicalStreamOrder(object):
    def __init__(self):
        self.label = "Topological Stream Order"
        self.description = "Assigns each link in a stream network its topological order."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.topological_stream_order(d8_pntr, streams, output, esri_pntr, zero_background))
        return


class TributaryIdentifier(object):
    def __init__(self):
        self.label = "Tributary Identifier"
        self.description = "Assigns a unique identifier to each tributary in a stream network."
        self.category = "Stream Network Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        streams = arcpy.Parameter(
            displayName="Input Streams File",
            name="streams",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        output = arcpy.Parameter(
            displayName="Output File",
            name="output",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        output.filter.list = ["tif"]

        esri_pntr = arcpy.Parameter(
            displayName="Does the pointer file use the ESRI pointer scheme?",
            name="esri_pntr",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        esri_pntr.value = "false"

        zero_background = arcpy.Parameter(
            displayName="Should a background value of zero be used?",
            name="zero_background",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        params = [d8_pntr, streams, output, esri_pntr, zero_background]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        d8_pntr = parameters[0].valueAsText
        streams = parameters[1].valueAsText
        output = parameters[2].valueAsText
        esri_pntr = parameters[3].valueAsText
        zero_background = parameters[4].valueAsText
        messages.addMessage(wbt.tributary_identifier(d8_pntr, streams, output, esri_pntr, zero_background))
        return



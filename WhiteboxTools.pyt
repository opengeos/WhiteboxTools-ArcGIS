import arcpy
import os
from WBT.whitebox_tools import WhiteboxTools
if sys.version_info < (3, 0):
    from StringIO import StringIO
else:
    from io import StringIO

wbt = WhiteboxTools()
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
tool_labels.append("Average Normal Vector Angular Deviation")
tool_labels.append("Average Overlay")
tool_labels.append("Average Upslope Flowpath Length")
tool_labels.append("Balance Contrast Enhancement")
tool_labels.append("Basins")
tool_labels.append("Bilateral Filter")
tool_labels.append("Block Maximum Gridding")
tool_labels.append("Block Minimum Gridding")
tool_labels.append("Boundary Shape Complexity")
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
tool_labels.append("Merge Line Segments")
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
tool_labels.append("Narrowness Index")
tool_labels.append("Nearest Neighbour Gridding")
tool_labels.append("Negate")
tool_labels.append("New Raster From Base")
tool_labels.append("Normal Vectors")
tool_labels.append("Normalized Difference Index")
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
tool_labels.append("Shape Complexity Index Raster")
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
tool_labels.append("Spherical Std Dev Of Normals")
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
        self.label = "WhiteboxTools Toolbox"
        self.alias = "WBT"

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
        tools.append(MergeLineSegments)
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
        tools.append(BoundaryShapeComplexity)
        tools.append(CompactnessRatio)
        tools.append(EdgeProportion)
        tools.append(ElongationRatio)
        tools.append(FindPatchOrClassEdgeCells)
        tools.append(HoleProportion)
        tools.append(LinearityIndex)
        tools.append(NarrownessIndex)
        tools.append(PatchOrientation)
        tools.append(PerimeterAreaRatio)
        tools.append(RadiusOfGyration)
        tools.append(RelatedCircumscribingCircle)
        tools.append(ShapeComplexityIndex)
        tools.append(ShapeComplexityIndexRaster)
        tools.append(Aspect)
        tools.append(AverageNormalVectorAngularDeviation)
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
        tools.append(SphericalStdDevOfNormals)
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
        tools.append(NormalizedDifferenceIndex)
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

        for index, tool in enumerate(sorted(tools)):
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

        tool_name.value = "Breach Depressions"
        tool_name.filter.type = "ValueList"
        tool_name.filter.list = tool_labels

        args = arcpy.Parameter(
            displayName="Arguments",
            name="agrs",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        args.value = '--dem="/path/to/DEM.tif"  --output="/path/to/output.tif"'

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
        dir_path = os.path.dirname(os.path.realpath(__file__))
        exe_path = os.path.join(dir_path, "WBT/whitebox_tools.exe")
        cmd = '{} --run={} {}'.format(exe_path, tool_name, args)
        if "-v" not in cmd:
            cmd = cmd + ' -v'
        messages.addMessage(cmd)
        messages.addMessage(os.popen(cmd).read().rstrip())
        return


class AddPointCoordinatesToTable(object):
    def __init__(self):
        self.label = "Add Point Coordinates To Table"
        self.description = "Modifies the attribute table of a point vector by adding fields containing each point's X and Y coordinates. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#AddPointCoordinatesToTable' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/add_point_coordinates_to_table.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.add_point_coordinates_to_table(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ConvertNodataToZero(object):
    def __init__(self):
        self.label = "Convert Nodata To Zero"
        self.description = "Converts nodata values in a raster to zero. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#ConvertNodataToZero' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/convert_nodata_to_zero.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.convert_nodata_to_zero(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ConvertRasterFormat(object):
    def __init__(self):
        self.label = "Convert Raster Format"
        self.description = "Converts raster data from one format to another. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#ConvertRasterFormat' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/convert_raster_format.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.convert_raster_format(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExportTableToCsv(object):
    def __init__(self):
        self.label = "Export Table To Csv"
        self.description = "Exports an attribute table to a CSV text file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#ExportTableToCsv' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/export_table_to_csv.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.export_table_to_csv(input, output=output, headers=headers)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class JoinTables(object):
    def __init__(self):
        self.label = "Join Tables"
        self.description = "Merge a vector's attribute table with another table based on a common field. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#JoinTables' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/join_tables.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.join_tables(input1, pkey=pkey, input2=input2, fkey=fkey, import_field=import_field)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LinesToPolygons(object):
    def __init__(self):
        self.label = "Lines To Polygons"
        self.description = "Converts vector polylines to polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#LinesToPolygons' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/lines_to_polygons.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lines_to_polygons(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MergeTableWithCsv(object):
    def __init__(self):
        self.label = "Merge Table With Csv"
        self.description = "Merge a vector's attribute table with a table contained within a CSV text file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#MergeTableWithCsv' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/merge_table_with_csv.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.merge_table_with_csv(input, pkey=pkey, csv=csv, fkey=fkey, import_field=import_field)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MergeVectors(object):
    def __init__(self):
        self.label = "Merge Vectors"
        self.description = "Combines two or more input vectors of the same ShapeType creating a single, new output vector. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#MergeVectors' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/merge_vectors.rs' target='_blank'>GitHub</a>."
        self.category = "Data Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Vector Files",
            name="inputs",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.merge_vectors(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MultiPartToSinglePart(object):
    def __init__(self):
        self.label = "Multi Part To Single Part"
        self.description = "Converts a vector file containing multi-part features into a vector containing only single-part features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#MultiPartToSinglePart' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/multipart_to_singlepart.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.multi_part_to_single_part(input, output=output, exclude_holes=exclude_holes)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NewRasterFromBase(object):
    def __init__(self):
        self.label = "New Raster From Base"
        self.description = "Creates a new raster using a base image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#NewRasterFromBase' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/new_raster.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.new_raster_from_base(base, output=output, value=value, data_type=data_type)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PolygonsToLines(object):
    def __init__(self):
        self.label = "Polygons To Lines"
        self.description = "Converts vector polygons to polylines. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#PolygonsToLines' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/polygons_to_lines.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.polygons_to_lines(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PrintGeoTiffTags(object):
    def __init__(self):
        self.label = "Print Geo Tiff Tags"
        self.description = "Prints the tags within a GeoTIFF. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#PrintGeoTiffTags' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/print_geotiff_tags.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.print_geo_tiff_tags(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterToVectorLines(object):
    def __init__(self):
        self.label = "Raster To Vector Lines"
        self.description = "Converts a raster lines features into a vector of the POLYLINE shapetype. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#RasterToVectorLines' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/raster_to_vector_lines.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_to_vector_lines(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterToVectorPoints(object):
    def __init__(self):
        self.label = "Raster To Vector Points"
        self.description = "Converts a raster dataset to a vector of the POINT shapetype. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#RasterToVectorPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/raster_to_vector_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_to_vector_points(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ReinitializeAttributeTable(object):
    def __init__(self):
        self.label = "Reinitialize Attribute Table"
        self.description = "Reinitializes a vector's attribute table deleting all fields but the feature ID (FID). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#ReinitializeAttributeTable' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/reinitialize_attribute_table.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.reinitialize_attribute_table(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RemovePolygonHoles(object):
    def __init__(self):
        self.label = "Remove Polygon Holes"
        self.description = "Removes holes within the features of a vector polygon file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#RemovePolygonHoles' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/remove_polygon_holes.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.remove_polygon_holes(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SetNodataValue(object):
    def __init__(self):
        self.label = "Set Nodata Value"
        self.description = "Assign a specified value in an input image to the NoData value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#SetNodataValue' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/set_nodata_value.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.set_nodata_value(input, output=output, back_value=back_value)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SinglePartToMultiPart(object):
    def __init__(self):
        self.label = "Single Part To Multi Part"
        self.description = "Converts a vector file containing multi-part features into a vector containing only single-part features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#SinglePartToMultiPart' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/singlepart_to_multipart.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.single_part_to_multi_part(input, field=field, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class VectorLinesToRaster(object):
    def __init__(self):
        self.label = "Vector Lines To Raster"
        self.description = "Converts a vector containing polylines into a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#VectorLinesToRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/vector_lines_to_raster.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.vector_lines_to_raster(input, field=field, output=output, nodata=nodata, cell_size=cell_size, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class VectorPointsToRaster(object):
    def __init__(self):
        self.label = "Vector Points To Raster"
        self.description = "Converts a vector containing points into a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#VectorPointsToRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/vector_points_to_raster.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.vector_points_to_raster(input, field=field, output=output, assign=assign, nodata=nodata, cell_size=cell_size, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class VectorPolygonsToRaster(object):
    def __init__(self):
        self.label = "Vector Polygons To Raster"
        self.description = "Converts a vector containing polygons into a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/data_tools.html#VectorPolygonsToRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/data_tools/vector_polygons_to_raster.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.vector_polygons_to_raster(input, field=field, output=output, nodata=nodata, cell_size=cell_size, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AggregateRaster(object):
    def __init__(self):
        self.label = "Aggregate Raster"
        self.description = "Aggregates a raster to a lower resolution. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#AggregateRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/aggregate_raster.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.aggregate_raster(input, output=output, agg_factor=agg_factor, type=type)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BlockMaximumGridding(object):
    def __init__(self):
        self.label = "Block Maximum Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using a block maximum scheme. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#BlockMaximumGridding' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/block_maximum.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.block_maximum_gridding(input, field=field, use_z=use_z, output=output, cell_size=cell_size, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BlockMinimumGridding(object):
    def __init__(self):
        self.label = "Block Minimum Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using a block minimum scheme. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#BlockMinimumGridding' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/block_minimum.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.block_minimum_gridding(input, field=field, use_z=use_z, output=output, cell_size=cell_size, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Centroid(object):
    def __init__(self):
        self.label = "Centroid"
        self.description = "Calculates the centroid, or average location, of raster polygon objects. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Centroid' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/centroid.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.centroid(input, output=output, text_output=text_output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CentroidVector(object):
    def __init__(self):
        self.label = "Centroid Vector"
        self.description = "Identifes the centroid point of a vector polyline or polygon feature or a group of vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CentroidVector' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/centroid_vector.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.centroid_vector(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Clump(object):
    def __init__(self):
        self.label = "Clump"
        self.description = "Groups cells that form discrete areas, assigning them unique identifiers. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Clump' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/clump.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.clump(input, output=output, diag=diag, zero_back=zero_back)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ConstructVectorTin(object):
    def __init__(self):
        self.label = "Construct Vector Tin"
        self.description = "Creates a vector triangular irregular network (TIN) for a set of vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ConstructVectorTin' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/construct_vector_tin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.construct_vector_tin(input, field=field, use_z=use_z, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CreateHexagonalVectorGrid(object):
    def __init__(self):
        self.label = "Create Hexagonal Vector Grid"
        self.description = "Creates a hexagonal vector grid. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CreateHexagonalVectorGrid' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/create_hexagonal_vector_grid.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.create_hexagonal_vector_grid(input, output=output, width=width, orientation=orientation)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CreatePlane(object):
    def __init__(self):
        self.label = "Create Plane"
        self.description = "Creates a raster image based on the equation for a simple plane. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CreatePlane' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/create_plane.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.create_plane(base, output=output, gradient=gradient, aspect=aspect, constant=constant)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CreateRectangularVectorGrid(object):
    def __init__(self):
        self.label = "Create Rectangular Vector Grid"
        self.description = "Creates a rectangular vector grid. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CreateRectangularVectorGrid' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/create_rectangular_vector_grid.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.create_rectangular_vector_grid(input, output=output, width=width, height=height, xorig=xorig, yorig=yorig)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Dissolve(object):
    def __init__(self):
        self.label = "Dissolve"
        self.description = "Removes the interior, or shared, boundaries within a vector polygon coverage. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Dissolve' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/dissolve.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.dissolve(input, field=field, output=output, snap=snap)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EliminateCoincidentPoints(object):
    def __init__(self):
        self.label = "Eliminate Coincident Points"
        self.description = "Removes any coincident, or nearly coincident, points from a vector points file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#EliminateCoincidentPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/eliminate_coincident_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.eliminate_coincident_points(input, output=output, tolerance=tolerance)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExtendVectorLines(object):
    def __init__(self):
        self.label = "Extend Vector Lines"
        self.description = "Extends vector lines by a specified distance. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ExtendVectorLines' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/extend_vector_lines.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.extend_vector_lines(input, output=output, dist=dist, extend=extend)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExtractNodes(object):
    def __init__(self):
        self.label = "Extract Nodes"
        self.description = "Converts vector lines or polygons into vertex points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ExtractNodes' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/extract_nodes.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.extract_nodes(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExtractRasterValuesAtPoints(object):
    def __init__(self):
        self.label = "Extract Raster Values At Points"
        self.description = "Extracts the values of raster(s) at vector point locations. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ExtractRasterValuesAtPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/extract_raster_values_at_points.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.extract_raster_values_at_points(inputs, points=points, out_text=out_text)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindLowestOrHighestPoints(object):
    def __init__(self):
        self.label = "Find Lowest Or Highest Points"
        self.description = "Locates the lowest and/or highest valued cells in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#FindLowestOrHighestPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/find_lowest_or_highest_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_lowest_or_highest_points(input, output=output, out_type=out_type)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class IdwInterpolation(object):
    def __init__(self):
        self.label = "Idw Interpolation"
        self.description = "Interpolates vector points into a raster surface using an inverse-distance weighted scheme. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#IdwInterpolation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/idw_interpolation.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.idw_interpolation(input, field=field, use_z=use_z, output=output, weight=weight, radius=radius, min_points=min_points, cell_size=cell_size, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LayerFootprint(object):
    def __init__(self):
        self.label = "Layer Footprint"
        self.description = "Creates a vector polygon footprint of the area covered by a raster grid or vector layer. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#LayerFootprint' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/layer_footprint.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.layer_footprint(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Medoid(object):
    def __init__(self):
        self.label = "Medoid"
        self.description = "Calculates the medoid for a series of vector features contained in a shapefile. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Medoid' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/medoid.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.medoid(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinimumBoundingBox(object):
    def __init__(self):
        self.label = "Minimum Bounding Box"
        self.description = "Creates a vector minimum bounding rectangle around vector features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MinimumBoundingBox' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/minimum_bounding_box.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.minimum_bounding_box(input, output=output, criterion=criterion, features=features)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinimumBoundingCircle(object):
    def __init__(self):
        self.label = "Minimum Bounding Circle"
        self.description = "Delineates the minimum bounding circle (i.e. smallest enclosing circle) for a group of vectors. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MinimumBoundingCircle' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/minimum_bounding_circle.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.minimum_bounding_circle(input, output=output, features=features)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinimumBoundingEnvelope(object):
    def __init__(self):
        self.label = "Minimum Bounding Envelope"
        self.description = "Creates a vector axis-aligned minimum bounding rectangle (envelope) around vector features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MinimumBoundingEnvelope' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/minimum_bounding_envelope.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.minimum_bounding_envelope(input, output=output, features=features)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinimumConvexHull(object):
    def __init__(self):
        self.label = "Minimum Convex Hull"
        self.description = "Creates a vector convex polygon around vector features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MinimumConvexHull' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/minimum_convex_hull.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.minimum_convex_hull(input, output=output, features=features)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NearestNeighbourGridding(object):
    def __init__(self):
        self.label = "Nearest Neighbour Gridding"
        self.description = "Creates a raster grid based on a set of vector points and assigns grid values using the nearest neighbour. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#NearestNeighbourGridding' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/nearest_neighbour_gridding.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        base = arcpy.Parameter(
            displayName="Base Raster File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.nearest_neighbour_gridding(input, field=field, use_z=use_z, output=output, cell_size=cell_size, base=base, max_dist=max_dist)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PolygonArea(object):
    def __init__(self):
        self.label = "Polygon Area"
        self.description = "Calculates the area of vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PolygonArea' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/polygon_area.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.polygon_area(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PolygonLongAxis(object):
    def __init__(self):
        self.label = "Polygon Long Axis"
        self.description = "This tool can be used to map the long axis of polygon features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PolygonLongAxis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/polygon_long_axis.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.polygon_long_axis(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PolygonPerimeter(object):
    def __init__(self):
        self.label = "Polygon Perimeter"
        self.description = "Calculates the perimeter of vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PolygonPerimeter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/polygon_perimeter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.polygon_perimeter(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PolygonShortAxis(object):
    def __init__(self):
        self.label = "Polygon Short Axis"
        self.description = "This tool can be used to map the short axis of polygon features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PolygonShortAxis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/polygon_short_axis.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.polygon_short_axis(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterArea(object):
    def __init__(self):
        self.label = "Raster Area"
        self.description = "Calculates the area of polygons or classes within a raster image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#RasterArea' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/raster_area.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_area(input, output=output, out_text=out_text, units=units, zero_back=zero_back)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterCellAssignment(object):
    def __init__(self):
        self.label = "Raster Cell Assignment"
        self.description = "Assign row or column number to cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#RasterCellAssignment' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/raster_cell_assignment.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_cell_assignment(input, output=output, assign=assign)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Reclass(object):
    def __init__(self):
        self.label = "Reclass"
        self.description = "Reclassifies the values in a raster image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Reclass' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/reclass.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.reclass(input, output=output, reclass_vals=reclass_vals, assign_mode=assign_mode)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ReclassEqualInterval(object):
    def __init__(self):
        self.label = "Reclass Equal Interval"
        self.description = "Reclassifies the values in a raster image based on equal-ranges. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ReclassEqualInterval' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/reclass_equal_interval.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.reclass_equal_interval(input, output=output, interval=interval, start_val=start_val, end_val=end_val)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ReclassFromFile(object):
    def __init__(self):
        self.label = "Reclass From File"
        self.description = "Reclassifies the values in a raster image using reclass ranges in a text file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ReclassFromFile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/reclass_from_file.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.reclass_from_file(input, reclass_file=reclass_file, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SmoothVectors(object):
    def __init__(self):
        self.label = "Smooth Vectors"
        self.description = "Smooths a vector coverage of either a POLYLINE or POLYGON base ShapeType. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#SmoothVectors' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/smooth_vectors.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.smooth_vectors(input, output=output, filter=filter)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TinGridding(object):
    def __init__(self):
        self.label = "Tin Gridding"
        self.description = "Creates a raster grid based on a triangular irregular network (TIN) fitted to vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#TinGridding' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/tin_gridding.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.tin_gridding(input, field=field, use_z=use_z, output=output, resolution=resolution)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class VectorHexBinning(object):
    def __init__(self):
        self.label = "Vector Hex Binning"
        self.description = "Hex-bins a set of vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#VectorHexBinning' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/vector_hex_bin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.vector_hex_binning(input, output=output, width=width, orientation=orientation)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class VoronoiDiagram(object):
    def __init__(self):
        self.label = "Voronoi Diagram"
        self.description = "Creates a vector Voronoi diagram for a set of vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#VoronoiDiagram' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/voronoi_diagram.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.voronoi_diagram(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BufferRaster(object):
    def __init__(self):
        self.label = "Buffer Raster"
        self.description = "Maps a distance-based buffer around each non-background (non-zero/non-nodata) grid cell in an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#BufferRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/buffer_raster.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.buffer_raster(input, output=output, size=size, gridcells=gridcells)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CostAllocation(object):
    def __init__(self):
        self.label = "Cost Allocation"
        self.description = "Identifies the source cell to which each grid cell is connected by a least-cost pathway in a cost-distance analysis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CostAllocation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/cost_allocation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cost_allocation(source, backlink=backlink, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CostDistance(object):
    def __init__(self):
        self.label = "Cost Distance"
        self.description = "Performs cost-distance accumulation on a cost surface and a group of source cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CostDistance' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/cost_distance.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cost_distance(source, cost=cost, out_accum=out_accum, out_backlink=out_backlink)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CostPathway(object):
    def __init__(self):
        self.label = "Cost Pathway"
        self.description = "Performs cost-distance pathway analysis using a series of destination grid cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CostPathway' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/cost_pathway.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cost_pathway(destination, backlink=backlink, output=output, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EuclideanAllocation(object):
    def __init__(self):
        self.label = "Euclidean Allocation"
        self.description = "Assigns grid cells in the output raster the value of the nearest target cell in the input image, measured by the Shih and Wu (2004) Euclidean distance transform. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#EuclideanAllocation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/euclidean_allocation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.euclidean_allocation(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EuclideanDistance(object):
    def __init__(self):
        self.label = "Euclidean Distance"
        self.description = "Calculates the Shih and Wu (2004) Euclidean distance transform. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#EuclideanDistance' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/euclidean_distance.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.euclidean_distance(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AverageOverlay(object):
    def __init__(self):
        self.label = "Average Overlay"
        self.description = "Calculates the average for each grid cell from a group of raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#AverageOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/average_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.average_overlay(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Clip(object):
    def __init__(self):
        self.label = "Clip"
        self.description = "Extract all the features, or parts of features, that overlap with the features of the clip vector. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Clip' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/clip.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.clip(input, clip=clip, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ClipRasterToPolygon(object):
    def __init__(self):
        self.label = "Clip Raster To Polygon"
        self.description = "Clips a raster to a vector polygon. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ClipRasterToPolygon' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/clip_raster_to_polygon.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.clip_raster_to_polygon(input, polygons=polygons, output=output, maintain_dimensions=maintain_dimensions)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CountIf(object):
    def __init__(self):
        self.label = "Count If"
        self.description = "Counts the number of occurrences of a specified value in a cell-stack of rasters. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CountIf' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/count_if.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.count_if(inputs, output=output, value=value)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Difference(object):
    def __init__(self):
        self.label = "Difference"
        self.description = "Outputs the features that occur in one of the two vector inputs but not both, i.e. no overlapping features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Difference' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/difference.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.difference(input, overlay=overlay, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Erase(object):
    def __init__(self):
        self.label = "Erase"
        self.description = "Removes all the features, or parts of features, that overlap with the features of the erase vector polygon. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Erase' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/erase.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.erase(input, erase=erase, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ErasePolygonFromRaster(object):
    def __init__(self):
        self.label = "Erase Polygon From Raster"
        self.description = "Erases (cuts out) a vector polygon from a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ErasePolygonFromRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/erase_polygon_from_raster.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.erase_polygon_from_raster(input, polygons=polygons, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HighestPosition(object):
    def __init__(self):
        self.label = "Highest Position"
        self.description = "Identifies the stack position of the maximum value within a raster stack on a cell-by-cell basis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#HighestPosition' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/highest_pos.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.highest_position(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Intersect(object):
    def __init__(self):
        self.label = "Intersect"
        self.description = "Identifies the parts of features in common between two input vector layers. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Intersect' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/intersect.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.intersect(input, overlay=overlay, output=output, snap=snap)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LineIntersections(object):
    def __init__(self):
        self.label = "Line Intersections"
        self.description = "Identifies points where the features of two vector line layers intersect. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#LineIntersections' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/line_intersections.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.line_intersections(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LowestPosition(object):
    def __init__(self):
        self.label = "Lowest Position"
        self.description = "Identifies the stack position of the minimum value within a raster stack on a cell-by-cell basis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#LowestPosition' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/lowest_pos.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lowest_position(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxAbsoluteOverlay(object):
    def __init__(self):
        self.label = "Max Absolute Overlay"
        self.description = "Evaluates the maximum absolute value for each grid cell from a stack of input rasters. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MaxAbsoluteOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/max_abs_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_absolute_overlay(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxOverlay(object):
    def __init__(self):
        self.label = "Max Overlay"
        self.description = "Evaluates the maximum value for each grid cell from a stack of input rasters. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MaxOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/max_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_overlay(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MergeLineSegments(object):
    def __init__(self):
        self.label = "Merge Line Segments"
        self.description = "Merges vector line segments into larger features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MergeLineSegments' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/merge_line_segments.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector File",
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

        snap = arcpy.Parameter(
            displayName="Snap Tolerance",
            name="snap",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        snap.value = "0.0"

        params = [input, output, snap]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        snap = parameters[2].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.merge_line_segments(input, output=output, snap=snap)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinAbsoluteOverlay(object):
    def __init__(self):
        self.label = "Min Absolute Overlay"
        self.description = "Evaluates the minimum absolute value for each grid cell from a stack of input rasters. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MinAbsoluteOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/min_abs_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.min_absolute_overlay(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinOverlay(object):
    def __init__(self):
        self.label = "Min Overlay"
        self.description = "Evaluates the minimum value for each grid cell from a stack of input rasters. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#MinOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/min_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.min_overlay(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PercentEqualTo(object):
    def __init__(self):
        self.label = "Percent Equal To"
        self.description = "Calculates the percentage of a raster stack that have cell values equal to an input on a cell-by-cell basis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PercentEqualTo' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/percent_equal_to.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.percent_equal_to(inputs, comparison=comparison, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PercentGreaterThan(object):
    def __init__(self):
        self.label = "Percent Greater Than"
        self.description = "Calculates the percentage of a raster stack that have cell values greather than an input on a cell-by-cell basis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PercentGreaterThan' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/percent_greater_than.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.percent_greater_than(inputs, comparison=comparison, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PercentLessThan(object):
    def __init__(self):
        self.label = "Percent Less Than"
        self.description = "Calculates the percentage of a raster stack that have cell values less than an input on a cell-by-cell basis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PercentLessThan' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/percent_less_than.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.percent_less_than(inputs, comparison=comparison, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PickFromList(object):
    def __init__(self):
        self.label = "Pick From List"
        self.description = "Outputs the value from a raster stack specified by a position raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PickFromList' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/pick_from_list.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.pick_from_list(inputs, pos_input=pos_input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Polygonize(object):
    def __init__(self):
        self.label = "Polygonize"
        self.description = "Creates a polygon layer from two or more intersecting line features contained in one or more input vector line files. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Polygonize' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/polygonize.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Vector Lines File",
            name="inputs",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.polygonize(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SplitWithLines(object):
    def __init__(self):
        self.label = "Split With Lines"
        self.description = "Splits the lines or polygons in one layer using the lines in another layer. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#SplitWithLines' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/split_with_lines.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.split_with_lines(input, split=split, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SumOverlay(object):
    def __init__(self):
        self.label = "Sum Overlay"
        self.description = "Calculates the sum for each grid cell from a group of raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#SumOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/sum_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sum_overlay(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SymmetricalDifference(object):
    def __init__(self):
        self.label = "Symmetrical Difference"
        self.description = "Outputs the features that occur in one of the two vector inputs but not both, i.e. no overlapping features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#SymmetricalDifference' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/symmetrical_difference.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.symmetrical_difference(input, overlay=overlay, output=output, snap=snap)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Union(object):
    def __init__(self):
        self.label = "Union"
        self.description = "Splits vector layers at their overlaps, creating a layer containing all the portions from both input and overlay layers. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#Union' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/union.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.union(input, overlay=overlay, output=output, snap=snap)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class WeightedOverlay(object):
    def __init__(self):
        self.label = "Weighted Overlay"
        self.description = "Performs a weighted sum on multiple input rasters after converting each image to a common scale. The tool performs a multi-criteria evaluation (MCE). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#WeightedOverlay' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/weighted_overlay.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        factors = arcpy.Parameter(
            displayName="Input Factor Files",
            name="factors",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        factors.multiValue = True

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
        constraints.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.weighted_overlay(factors, weights=weights, cost=cost, constraints=constraints, output=output, scale_max=scale_max)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class WeightedSum(object):
    def __init__(self):
        self.label = "Weighted Sum"
        self.description = "Performs a weighted-sum overlay on multiple input raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#WeightedSum' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/weighted_sum.rs' target='_blank'>GitHub</a>."
        self.category = "GIS Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.weighted_sum(inputs, weights=weights, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BoundaryShapeComplexity(object):
    def __init__(self):
        self.label = "Boundary Shape Complexity"
        self.description = "Calculates the complexity of the boundaries of raster polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#BoundaryShapeComplexity' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/boundary_shape_complexity.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.boundary_shape_complexity(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CompactnessRatio(object):
    def __init__(self):
        self.label = "Compactness Ratio"
        self.description = "Calculates the compactness ratio (A/P), a measure of shape complexity, for vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#CompactnessRatio' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/compactness_ratio.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.compactness_ratio(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EdgeProportion(object):
    def __init__(self):
        self.label = "Edge Proportion"
        self.description = "Calculate the proportion of cells in a raster polygon that are edge cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#EdgeProportion' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/edge_proportion.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.edge_proportion(input, output=output, output_text=output_text)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElongationRatio(object):
    def __init__(self):
        self.label = "Elongation Ratio"
        self.description = "Calculates the elongation ratio for vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ElongationRatio' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/elongation_ratio.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elongation_ratio(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindPatchOrClassEdgeCells(object):
    def __init__(self):
        self.label = "Find Patch Or Class Edge Cells"
        self.description = "Finds all cells located on the edge of patch or class features. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#FindPatchOrClassEdgeCells' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/find_patch_edge_cells.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_patch_or_class_edge_cells(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HoleProportion(object):
    def __init__(self):
        self.label = "Hole Proportion"
        self.description = "Calculates the proportion of the total area of a polygon's holes relative to the area of the polygon's hull. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#HoleProportion' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/hole_proportion.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.hole_proportion(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LinearityIndex(object):
    def __init__(self):
        self.label = "Linearity Index"
        self.description = "Calculates the linearity index for vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#LinearityIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/linearity_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.linearity_index(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NarrownessIndex(object):
    def __init__(self):
        self.label = "Narrowness Index"
        self.description = "Calculates the narrowness of raster polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#NarrownessIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/narrowness_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.narrowness_index(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PatchOrientation(object):
    def __init__(self):
        self.label = "Patch Orientation"
        self.description = "Calculates the orientation of vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PatchOrientation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/patch_orientation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.patch_orientation(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PerimeterAreaRatio(object):
    def __init__(self):
        self.label = "Perimeter Area Ratio"
        self.description = "Calculates the perimeter-area ratio of vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#PerimeterAreaRatio' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/perimeter_area_ratio.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.perimeter_area_ratio(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RadiusOfGyration(object):
    def __init__(self):
        self.label = "Radius Of Gyration"
        self.description = "Calculates the distance of cells from their polygon's centroid. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#RadiusOfGyration' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/radius_of_gyration.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.radius_of_gyration(input, output=output, text_output=text_output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RelatedCircumscribingCircle(object):
    def __init__(self):
        self.label = "Related Circumscribing Circle"
        self.description = "Calculates the related circumscribing circle of vector polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#RelatedCircumscribingCircle' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/related_circumscribing_circle.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.related_circumscribing_circle(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ShapeComplexityIndex(object):
    def __init__(self):
        self.label = "Shape Complexity Index"
        self.description = "Calculates overall polygon shape complexity or irregularity. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ShapeComplexityIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/shape_complexity_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.shape_complexity_index(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ShapeComplexityIndexRaster(object):
    def __init__(self):
        self.label = "Shape Complexity Index Raster"
        self.description = "Calculates the complexity of raster polygons or classes. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html#ShapeComplexityIndexRaster' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/gis_analysis/shape_complexity_raster.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.shape_complexity_index_raster(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Aspect(object):
    def __init__(self):
        self.label = "Aspect"
        self.description = "Calculates an aspect raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#Aspect' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/aspect.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.aspect(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AverageNormalVectorAngularDeviation(object):
    def __init__(self):
        self.label = "Average Normal Vector Angular Deviation"
        self.description = "Calculates the circular variance of aspect at a scale for a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#AverageNormalVectorAngularDeviation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/average_normal_vector_angular_deviation.rs' target='_blank'>GitHub</a>."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.average_normal_vector_angular_deviation(dem, output=output, filter=filter)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CircularVarianceOfAspect(object):
    def __init__(self):
        self.label = "Circular Variance Of Aspect"
        self.description = "Calculates the circular variance of aspect at a scale for a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#CircularVarianceOfAspect' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/circular_variance_of_aspect.rs' target='_blank'>GitHub</a>."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.circular_variance_of_aspect(dem, output=output, filter=filter)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DevFromMeanElev(object):
    def __init__(self):
        self.label = "Dev From Mean Elev"
        self.description = "Calculates deviation from mean elevation. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#DevFromMeanElev' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/dev_from_mean_elev.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.dev_from_mean_elev(dem, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DiffFromMeanElev(object):
    def __init__(self):
        self.label = "Diff From Mean Elev"
        self.description = "Calculates difference from mean elevation (equivalent to a high-pass filter). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#DiffFromMeanElev' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/diff_from_mean_elev.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.diff_from_mean_elev(dem, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DirectionalRelief(object):
    def __init__(self):
        self.label = "Directional Relief"
        self.description = "Calculates relief for cells in an input DEM for a specified direction. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#DirectionalRelief' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/directional_relief.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.directional_relief(dem, output=output, azimuth=azimuth, max_dist=max_dist)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DownslopeIndex(object):
    def __init__(self):
        self.label = "Downslope Index"
        self.description = "Calculates the Hjerdt et al. (2004) downslope index. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#DownslopeIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/downslope_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.downslope_index(dem, output=output, drop=drop, out_type=out_type)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DrainagePreservingSmoothing(object):
    def __init__(self):
        self.label = "Drainage Preserving Smoothing"
        self.description = "Reduces short-scale variation in an input DEM while preserving breaks-in-slope and small drainage features using a modified Sun et al. (2007) algorithm. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#DrainagePreservingSmoothing' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/drainage_preserving_smoothing.rs' target='_blank'>GitHub</a>."
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

        max_diff.value = "0.5"

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.drainage_preserving_smoothing(dem, output=output, filter=filter, norm_diff=norm_diff, num_iter=num_iter, max_diff=max_diff, reduction=reduction, dfm=dfm, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EdgeDensity(object):
    def __init__(self):
        self.label = "Edge Density"
        self.description = "Calculates the density of edges, or breaks-in-slope within DEMs. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#EdgeDensity' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/edge_density.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.edge_density(dem, output=output, filter=filter, norm_diff=norm_diff, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElevAbovePit(object):
    def __init__(self):
        self.label = "Elev Above Pit"
        self.description = "Calculate the elevation of each grid cell above the nearest downstream pit cell or grid edge cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#ElevAbovePit' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/elev_above_pit.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elev_above_pit(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElevPercentile(object):
    def __init__(self):
        self.label = "Elev Percentile"
        self.description = "Calculates the elevation percentile raster from a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#ElevPercentile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/elev_percentile.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elev_percentile(dem, output=output, filterx=filterx, filtery=filtery, sig_digits=sig_digits)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElevRelativeToMinMax(object):
    def __init__(self):
        self.label = "Elev Relative To Min Max"
        self.description = "Calculates the elevation of a location relative to the minimum and maximum elevations in a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#ElevRelativeToMinMax' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/elev_relative_to_min_max.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elev_relative_to_min_max(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElevRelativeToWatershedMinMax(object):
    def __init__(self):
        self.label = "Elev Relative To Watershed Min Max"
        self.description = "Calculates the elevation of a location relative to the minimum and maximum elevations in a watershed. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#ElevRelativeToWatershedMinMax' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/elev_relative_to_watershed_min_max.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elev_relative_to_watershed_min_max(dem, watersheds=watersheds, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FeaturePreservingDenoise(object):
    def __init__(self):
        self.label = "Feature Preserving Denoise"
        self.description = "Reduces short-scale variation in an input DEM using a modified Sun et al. (2007) algorithm. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#FeaturePreservingDenoise' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/feature_preserving_denoise.rs' target='_blank'>GitHub</a>."
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

        max_diff.value = "0.5"

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.feature_preserving_denoise(dem, output=output, filter=filter, norm_diff=norm_diff, num_iter=num_iter, max_diff=max_diff, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FetchAnalysis(object):
    def __init__(self):
        self.label = "Fetch Analysis"
        self.description = "Performs an analysis of fetch or upwind distance to an obstacle. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#FetchAnalysis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/fetch_analysis.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fetch_analysis(dem, output=output, azimuth=azimuth, hgt_inc=hgt_inc)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FillMissingData(object):
    def __init__(self):
        self.label = "Fill Missing Data"
        self.description = "Fills NoData holes in a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#FillMissingData' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/fill_missing_data.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fill_missing_data(input, output=output, filter=filter, weight=weight)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindRidges(object):
    def __init__(self):
        self.label = "Find Ridges"
        self.description = "Identifies potential ridge and peak grid cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#FindRidges' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/find_ridges.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_ridges(dem, output=output, line_thin=line_thin)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Hillshade(object):
    def __init__(self):
        self.label = "Hillshade"
        self.description = "Calculates a hillshade raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#Hillshade' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/hillshade.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.hillshade(dem, output=output, azimuth=azimuth, altitude=altitude, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HorizonAngle(object):
    def __init__(self):
        self.label = "Horizon Angle"
        self.description = "Calculates horizon angle (maximum upwind slope) for each grid cell in an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#HorizonAngle' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/horizon_angle.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.horizon_angle(dem, output=output, azimuth=azimuth, max_dist=max_dist)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HypsometricAnalysis(object):
    def __init__(self):
        self.label = "Hypsometric Analysis"
        self.description = "Calculates a hypsometric curve for one or more DEMs. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#HypsometricAnalysis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/hypsometric_analysis.rs' target='_blank'>GitHub</a>."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input DEM Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

        watershed = arcpy.Parameter(
            displayName="Input Watershed Files",
            name="watershed",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")
        watershed.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.hypsometric_analysis(inputs, watershed=watershed, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxAnisotropyDev(object):
    def __init__(self):
        self.label = "Max Anisotropy Dev"
        self.description = "Calculates the maximum anisotropy (directionality) in elevation deviation over a range of spatial scales. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxAnisotropyDev' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_anisotropy_dev.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_anisotropy_dev(dem, out_mag=out_mag, out_scale=out_scale, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxAnisotropyDevSignature(object):
    def __init__(self):
        self.label = "Max Anisotropy Dev Signature"
        self.description = "Calculates the anisotropy in deviation from mean for points over a range of spatial scales. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxAnisotropyDevSignature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_anisotropy_dev_signature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_anisotropy_dev_signature(dem, points=points, output=output, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxBranchLength(object):
    def __init__(self):
        self.label = "Max Branch Length"
        self.description = "Lindsay and Seibert's (2013) branch length index is used to map drainage divides or ridge lines. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxBranchLength' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_branch_length.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_branch_length(dem, output=output, log=log)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxDifferenceFromMean(object):
    def __init__(self):
        self.label = "Max Difference From Mean"
        self.description = "Calculates the maximum difference from mean elevation over a range of spatial scales. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxDifferenceFromMean' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_diff_from_mean.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_difference_from_mean(dem, out_mag=out_mag, out_scale=out_scale, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxDownslopeElevChange(object):
    def __init__(self):
        self.label = "Max Downslope Elev Change"
        self.description = "Calculates the maximum downslope change in elevation between a grid cell and its eight downslope neighbors. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxDownslopeElevChange' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_downslope_elev_change.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_downslope_elev_change(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxElevDevSignature(object):
    def __init__(self):
        self.label = "Max Elev Dev Signature"
        self.description = "Calculates the maximum elevation deviation over a range of spatial scales and for a set of points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxElevDevSignature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_elev_dev_signature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_elev_dev_signature(dem, points=points, output=output, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxElevationDeviation(object):
    def __init__(self):
        self.label = "Max Elevation Deviation"
        self.description = "Calculates the maximum elevation deviation over a range of spatial scales. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MaxElevationDeviation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/max_elev_deviation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_elevation_deviation(dem, out_mag=out_mag, out_scale=out_scale, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinDownslopeElevChange(object):
    def __init__(self):
        self.label = "Min Downslope Elev Change"
        self.description = "Calculates the minimum downslope change in elevation between a grid cell and its eight downslope neighbors. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MinDownslopeElevChange' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/min_downslope_elev_change.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.min_downslope_elev_change(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MultiscaleRoughness(object):
    def __init__(self):
        self.label = "Multiscale Roughness"
        self.description = "Calculates surface roughness over a range of spatial scales. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MultiscaleRoughness' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/multiscale_roughness.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.multiscale_roughness(dem, out_mag=out_mag, out_scale=out_scale, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MultiscaleRoughnessSignature(object):
    def __init__(self):
        self.label = "Multiscale Roughness Signature"
        self.description = "Calculates the surface roughness for points over a range of spatial scales. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MultiscaleRoughnessSignature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/multiscale_roughness_signature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.multiscale_roughness_signature(dem, points=points, output=output, min_scale=min_scale, max_scale=max_scale, step=step)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MultiscaleTopographicPositionImage(object):
    def __init__(self):
        self.label = "Multiscale Topographic Position Image"
        self.description = "Creates a multiscale topographic position image from three DEVmax rasters of differing spatial scale ranges. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#MultiscaleTopographicPositionImage' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/multiscale_topographic_position_image.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.multiscale_topographic_position_image(local, meso=meso, broad=broad, output=output, lightness=lightness)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NumDownslopeNeighbours(object):
    def __init__(self):
        self.label = "Num Downslope Neighbours"
        self.description = "Calculates the number of downslope neighbours to each grid cell in a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#NumDownslopeNeighbours' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/num_downslope_neighbours.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.num_downslope_neighbours(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NumUpslopeNeighbours(object):
    def __init__(self):
        self.label = "Num Upslope Neighbours"
        self.description = "Calculates the number of upslope neighbours to each grid cell in a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#NumUpslopeNeighbours' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/num_upslope_neighbours.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.num_upslope_neighbours(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PennockLandformClass(object):
    def __init__(self):
        self.label = "Pennock Landform Class"
        self.description = "Classifies hillslope zones based on slope, profile curvature, and plan curvature. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#PennockLandformClass' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/pennock_landform_class.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.pennock_landform_class(dem, output=output, slope=slope, prof=prof, plan=plan, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PercentElevRange(object):
    def __init__(self):
        self.label = "Percent Elev Range"
        self.description = "Calculates percent of elevation range from a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#PercentElevRange' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/percent_elev_range.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.percent_elev_range(dem, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PlanCurvature(object):
    def __init__(self):
        self.label = "Plan Curvature"
        self.description = "Calculates a plan (contour) curvature raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#PlanCurvature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/plan_curvature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.plan_curvature(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Profile(object):
    def __init__(self):
        self.label = "Profile"
        self.description = "Plots profiles from digital surface models. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#Profile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/profile.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.profile(lines, surface=surface, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ProfileCurvature(object):
    def __init__(self):
        self.label = "Profile Curvature"
        self.description = "Calculates a profile curvature raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#ProfileCurvature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/prof_curvature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.profile_curvature(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RelativeAspect(object):
    def __init__(self):
        self.label = "Relative Aspect"
        self.description = "Calculates relative aspect (relative to a user-specified direction) from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#RelativeAspect' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/relative_aspect.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.relative_aspect(dem, output=output, azimuth=azimuth, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RelativeStreamPowerIndex(object):
    def __init__(self):
        self.label = "Relative Stream Power Index"
        self.description = "Calculates the relative stream power index. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#RelativeStreamPowerIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/relative_stream_power_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.relative_stream_power_index(sca, slope=slope, output=output, exponent=exponent)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RelativeTopographicPosition(object):
    def __init__(self):
        self.label = "Relative Topographic Position"
        self.description = "Calculates the relative topographic position index from a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#RelativeTopographicPosition' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/relative_topographic_position.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.relative_topographic_position(dem, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RemoveOffTerrainObjects(object):
    def __init__(self):
        self.label = "Remove Off Terrain Objects"
        self.description = "Removes off-terrain objects from a raster digital elevation model (DEM). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#RemoveOffTerrainObjects' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/remove_off_terrain_objects.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.remove_off_terrain_objects(dem, output=output, filter=filter, slope=slope)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RuggednessIndex(object):
    def __init__(self):
        self.label = "Ruggedness Index"
        self.description = "Calculates the Riley et al.'s (1999) terrain ruggedness index from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#RuggednessIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/ruggedness_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.ruggedness_index(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SedimentTransportIndex(object):
    def __init__(self):
        self.label = "Sediment Transport Index"
        self.description = "Calculates the sediment transport index. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#SedimentTransportIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/sediment_transport_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sediment_transport_index(sca, slope=slope, output=output, sca_exponent=sca_exponent, slope_exponent=slope_exponent)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Slope(object):
    def __init__(self):
        self.label = "Slope"
        self.description = "Calculates a slope raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#Slope' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/slope.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.slope(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SlopeVsElevationPlot(object):
    def __init__(self):
        self.label = "Slope Vs Elevation Plot"
        self.description = "Creates a slope vs. elevation plot for one or more DEMs. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#SlopeVsElevationPlot' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/slope_vs_elev_plot.rs' target='_blank'>GitHub</a>."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input DEM Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

        watershed = arcpy.Parameter(
            displayName="Input Watershed Files",
            name="watershed",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")
        watershed.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.slope_vs_elevation_plot(inputs, watershed=watershed, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SphericalStdDevOfNormals(object):
    def __init__(self):
        self.label = "Spherical Std Dev Of Normals"
        self.description = "Calculates the spherical standard deviation of surface normals for a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#SphericalStdDevOfNormals' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/spherical_std_dev_of_normals.rs' target='_blank'>GitHub</a>."
        self.category = "Geomorphometric Analysis"

    def getParameterInfo(self):
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.spherical_std_dev_of_normals(dem, output=output, filter=filter)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StandardDeviationOfSlope(object):
    def __init__(self):
        self.label = "Standard Deviation Of Slope"
        self.description = "Calculates the standard deviation of slope from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#StandardDeviationOfSlope' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/standard_deviation_of_slope.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.standard_deviation_of_slope(input, output=output, zfactor=zfactor, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SurfaceAreaRatio(object):
    def __init__(self):
        self.label = "Surface Area Ratio"
        self.description = "Calculates a the surface area ratio of each grid cell in an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#SurfaceAreaRatio' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/surface_area_ratio.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.surface_area_ratio(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TangentialCurvature(object):
    def __init__(self):
        self.label = "Tangential Curvature"
        self.description = "Calculates a tangential curvature raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#TangentialCurvature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/tan_curvature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.tangential_curvature(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TotalCurvature(object):
    def __init__(self):
        self.label = "Total Curvature"
        self.description = "Calculates a total curvature raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#TotalCurvature' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/total_curvature.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.total_curvature(dem, output=output, zfactor=zfactor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Viewshed(object):
    def __init__(self):
        self.label = "Viewshed"
        self.description = "Identifies the viewshed for a point or set of points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#Viewshed' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/viewshed.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.viewshed(dem, stations=stations, output=output, height=height)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class VisibilityIndex(object):
    def __init__(self):
        self.label = "Visibility Index"
        self.description = "Estimates the relative visibility of sites in a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#VisibilityIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/visibility_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.visibility_index(dem, output=output, height=height, res_factor=res_factor)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class WetnessIndex(object):
    def __init__(self):
        self.label = "Wetness Index"
        self.description = "Calculates the topographic wetness index, Ln(A / tan(slope)). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html#WetnessIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/terrain_analysis/wetness_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.wetness_index(sca, slope=slope, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AverageFlowpathSlope(object):
    def __init__(self):
        self.label = "Average Flowpath Slope"
        self.description = "Measures the average slope gradient from each grid cell to all upslope divide cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#AverageFlowpathSlope' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/average_flowpath_slope.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.average_flowpath_slope(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AverageUpslopeFlowpathLength(object):
    def __init__(self):
        self.label = "Average Upslope Flowpath Length"
        self.description = "Measures the average length of all upslope flowpaths draining each grid cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#AverageUpslopeFlowpathLength' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/average_upslope_flowpath_length.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.average_upslope_flowpath_length(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Basins(object):
    def __init__(self):
        self.label = "Basins"
        self.description = "Identifies drainage basins that drain to the DEM edge. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Basins' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/basins.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.basins(d8_pntr, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BreachDepressions(object):
    def __init__(self):
        self.label = "Breach Depressions"
        self.description = "Breaches all of the depressions in a DEM using Lindsay's (2016) algorithm. This should be preferred over depression filling in most cases. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#BreachDepressions' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/breach_depressions.rs' target='_blank'>GitHub</a>."
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

        flat_increment = arcpy.Parameter(
            displayName="Flat increment value (z units)",
            name="flat_increment",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        fill_pits = arcpy.Parameter(
            displayName="Fill single-cell pits?",
            name="fill_pits",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        fill_pits.value = "false"

        params = [dem, output, max_depth, max_length, flat_increment, fill_pits]

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
        flat_increment = parameters[4].valueAsText
        fill_pits = parameters[5].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.breach_depressions(dem, output=output, max_depth=max_depth, max_length=max_length, flat_increment=flat_increment, fill_pits=fill_pits)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BreachSingleCellPits(object):
    def __init__(self):
        self.label = "Breach Single Cell Pits"
        self.description = "Removes single-cell pits from an input DEM by breaching. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#BreachSingleCellPits' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/breach_pits.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.breach_single_cell_pits(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class D8FlowAccumulation(object):
    def __init__(self):
        self.label = "D8 Flow Accumulation"
        self.description = "Calculates a D8 flow accumulation raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#D8FlowAccumulation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/d8_flow_accum.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.d8_flow_accumulation(dem, output=output, out_type=out_type, log=log, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class D8MassFlux(object):
    def __init__(self):
        self.label = "D8 Mass Flux"
        self.description = "Performs a D8 mass flux calculation. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#D8MassFlux' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/d8_mass_flux.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.d8_mass_flux(dem, loading=loading, efficiency=efficiency, absorption=absorption, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class D8Pointer(object):
    def __init__(self):
        self.label = "D8 Pointer"
        self.description = "Calculates a D8 flow pointer raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#D8Pointer' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/d8_pointer.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.d8_pointer(dem, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DInfFlowAccumulation(object):
    def __init__(self):
        self.label = "D Inf Flow Accumulation"
        self.description = "Calculates a D-infinity flow accumulation raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#DInfFlowAccumulation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/dinf_flow_accum.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.d_inf_flow_accumulation(dem, output=output, out_type=out_type, threshold=threshold, log=log, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DInfMassFlux(object):
    def __init__(self):
        self.label = "D Inf Mass Flux"
        self.description = "Performs a D-infinity mass flux calculation. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#DInfMassFlux' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/dinf_mass_flux.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.d_inf_mass_flux(dem, loading=loading, efficiency=efficiency, absorption=absorption, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DInfPointer(object):
    def __init__(self):
        self.label = "D Inf Pointer"
        self.description = "Calculates a D-infinity flow pointer (flow direction) raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#DInfPointer' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/dinf_pointer.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.d_inf_pointer(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DepthInSink(object):
    def __init__(self):
        self.label = "Depth In Sink"
        self.description = "Measures the depth of sinks (depressions) in a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#DepthInSink' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/depth_in_sink.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.depth_in_sink(dem, output=output, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DownslopeDistanceToStream(object):
    def __init__(self):
        self.label = "Downslope Distance To Stream"
        self.description = "Measures distance to the nearest downslope stream cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#DownslopeDistanceToStream' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/downslope_distance_to_stream.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.downslope_distance_to_stream(dem, streams=streams, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DownslopeFlowpathLength(object):
    def __init__(self):
        self.label = "Downslope Flowpath Length"
        self.description = "Calculates the downslope flowpath length from each cell to basin outlet. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#DownslopeFlowpathLength' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/downslope_flowpath_length.rs' target='_blank'>GitHub</a>."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        d8_pntr = arcpy.Parameter(
            displayName="Input D8 Pointer File",
            name="d8_pntr",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        watersheds = arcpy.Parameter(
            displayName="Input Watersheds File",
            name="watersheds",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        weights = arcpy.Parameter(
            displayName="Input Weights File",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.downslope_flowpath_length(d8_pntr, watersheds=watersheds, weights=weights, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElevationAboveStream(object):
    def __init__(self):
        self.label = "Elevation Above Stream"
        self.description = "Calculates the elevation of cells above the nearest downslope stream cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#ElevationAboveStream' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/elevation_above_stream.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elevation_above_stream(dem, streams=streams, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ElevationAboveStreamEuclidean(object):
    def __init__(self):
        self.label = "Elevation Above Stream Euclidean"
        self.description = "Calculates the elevation of cells above the nearest (Euclidean distance) stream cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#ElevationAboveStreamEuclidean' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/elevation_above_stream_euclidean.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.elevation_above_stream_euclidean(dem, streams=streams, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Fd8FlowAccumulation(object):
    def __init__(self):
        self.label = "Fd8 Flow Accumulation"
        self.description = "Calculates an FD8 flow accumulation raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Fd8FlowAccumulation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/fd8_flow_accum.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fd8_flow_accumulation(dem, output=output, out_type=out_type, exponent=exponent, threshold=threshold, log=log, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Fd8Pointer(object):
    def __init__(self):
        self.label = "Fd8 Pointer"
        self.description = "Calculates an FD8 flow pointer raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Fd8Pointer' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/fd8_pointer.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fd8_pointer(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FillBurn(object):
    def __init__(self):
        self.label = "Fill Burn"
        self.description = "Burns streams into a DEM using the FillBurn (Saunders, 1999) method. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FillBurn' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/fill_burn.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fill_burn(dem, streams=streams, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FillDepressions(object):
    def __init__(self):
        self.label = "Fill Depressions"
        self.description = "Fills all of the depressions in a DEM. Depression breaching should be preferred in most cases. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FillDepressions' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/fill_depressions.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fill_depressions(dem, output=output, fix_flats=fix_flats)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FillSingleCellPits(object):
    def __init__(self):
        self.label = "Fill Single Cell Pits"
        self.description = "Raises pit cells to the elevation of their lowest neighbour. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FillSingleCellPits' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/fill_pits.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fill_single_cell_pits(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindNoFlowCells(object):
    def __init__(self):
        self.label = "Find No Flow Cells"
        self.description = "Finds grid cells with no downslope neighbours. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FindNoFlowCells' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/find_noflow_cells.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_no_flow_cells(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindParallelFlow(object):
    def __init__(self):
        self.label = "Find Parallel Flow"
        self.description = "Finds areas of parallel flow in D8 flow direction rasters. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FindParallelFlow' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/find_parallel_flow.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_parallel_flow(d8_pntr, streams=streams, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FlattenLakes(object):
    def __init__(self):
        self.label = "Flatten Lakes"
        self.description = "Flattens lake polygons in a raster DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FlattenLakes' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/flatten_lakes.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.flatten_lakes(dem, lakes=lakes, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FloodOrder(object):
    def __init__(self):
        self.label = "Flood Order"
        self.description = "Assigns each DEM grid cell its order in the sequence of inundations that are encountered during a search starting from the edges, moving inward at increasing elevations. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FloodOrder' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/flood_order.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.flood_order(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FlowAccumulationFullWorkflow(object):
    def __init__(self):
        self.label = "Flow Accumulation Full Workflow"
        self.description = "Resolves all of the depressions in a DEM, outputting a breached DEM, an aspect-aligned non-divergent flow pointer, and a flow accumulation raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FlowAccumulationFullWorkflow' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/flow_accum_full_workflow.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.flow_accumulation_full_workflow(dem, out_dem=out_dem, out_pntr=out_pntr, out_accum=out_accum, out_type=out_type, log=log, clip=clip, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FlowLengthDiff(object):
    def __init__(self):
        self.label = "Flow Length Diff"
        self.description = "Calculates the local maximum absolute difference in downslope flowpath length, useful in mapping drainage divides and ridges. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#FlowLengthDiff' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/flow_length_diff.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.flow_length_diff(d8_pntr, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Hillslopes(object):
    def __init__(self):
        self.label = "Hillslopes"
        self.description = "Identifies the individual hillslopes draining to each link in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Hillslopes' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/hillslopes.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.hillslopes(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ImpoundmentSizeIndex(object):
    def __init__(self):
        self.label = "Impoundment Size Index"
        self.description = "Calculates the impoundment size resulting from damming a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#ImpoundmentSizeIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/impoundment_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.impoundment_size_index(dem, output=output, out_type=out_type, damlength=damlength)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Isobasins(object):
    def __init__(self):
        self.label = "Isobasins"
        self.description = "Divides a landscape into nearly equal sized drainage basins (i.e. watersheds). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Isobasins' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/isobasins.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.isobasins(dem, output=output, size=size)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class JensonSnapPourPoints(object):
    def __init__(self):
        self.label = "Jenson Snap Pour Points"
        self.description = "Moves outlet points used to specify points of interest in a watershedding operation to the nearest stream cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#JensonSnapPourPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/jenson_snap_pour_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.jenson_snap_pour_points(pour_pts, streams=streams, output=output, snap_dist=snap_dist)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LongestFlowpath(object):
    def __init__(self):
        self.label = "Longest Flowpath"
        self.description = "Delineates the longest flowpaths for a group of subbasins or watersheds. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#LongestFlowpath' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/longest_flowpath.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.longest_flowpath(dem, basins=basins, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaxUpslopeFlowpathLength(object):
    def __init__(self):
        self.label = "Max Upslope Flowpath Length"
        self.description = "Measures the maximum length of all upslope flowpaths draining each grid cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#MaxUpslopeFlowpathLength' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/max_upslope_flowpath.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max_upslope_flowpath_length(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NumInflowingNeighbours(object):
    def __init__(self):
        self.label = "Num Inflowing Neighbours"
        self.description = "Computes the number of inflowing neighbours to each cell in an input DEM based on the D8 algorithm. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#NumInflowingNeighbours' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/num_inflowing_neighbours.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.num_inflowing_neighbours(dem, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RaiseWalls(object):
    def __init__(self):
        self.label = "Raise Walls"
        self.description = "Raises walls in a DEM along a line or around a polygon, e.g. a watershed. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#RaiseWalls' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/raise_walls.rs' target='_blank'>GitHub</a>."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Vector Line or Polygon File",
            name="input",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        breach = arcpy.Parameter(
            displayName="Input Breach Lines",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raise_walls(input, breach=breach, dem=dem, output=output, height=height)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Rho8Pointer(object):
    def __init__(self):
        self.label = "Rho8 Pointer"
        self.description = "Calculates a stochastic Rho8 flow pointer raster from an input DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Rho8Pointer' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/rho8_pointer.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.rho8_pointer(dem, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Sink(object):
    def __init__(self):
        self.label = "Sink"
        self.description = "Identifies the depressions in a DEM, giving each feature a unique identifier. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Sink' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/sink.rs' target='_blank'>GitHub</a>."
        self.category = "Hydrological Analysis"

    def getParameterInfo(self):
        # 9/8/19 - Parameter type set to GPRasterLayer, originally it was DERasterDataset
        dem = arcpy.Parameter(
            displayName="Input DEM File",
            name="dem",
            datatype="GPRasterLayer",
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
        dem = ConvertToFullPath(dem)
        output = parameters[1].valueAsText
        zero_background = parameters[2].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sink(dem, output=output, zero_background=zero_background)

        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SnapPourPoints(object):
    def __init__(self):
        self.label = "Snap Pour Points"
        self.description = "Moves outlet points used to specify points of interest in a watershedding operation to the cell with the highest flow accumulation in its neighbourhood. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#SnapPourPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/snap_pour_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.snap_pour_points(pour_pts, flow_accum=flow_accum, output=output, snap_dist=snap_dist)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StochasticDepressionAnalysis(object):
    def __init__(self):
        self.label = "Stochastic Depression Analysis"
        self.description = "Preforms a stochastic analysis of depressions within a DEM. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#StochasticDepressionAnalysis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/stochastic_depression_analysis.rs' target='_blank'>GitHub</a>."
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

        iterations.value = "100"

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.stochastic_depression_analysis(dem, output=output, rmse=rmse, range=range, iterations=iterations)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StrahlerOrderBasins(object):
    def __init__(self):
        self.label = "Strahler Order Basins"
        self.description = "Identifies Strahler-order basins from an input stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#StrahlerOrderBasins' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/strahler_basins.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.strahler_order_basins(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Subbasins(object):
    def __init__(self):
        self.label = "Subbasins"
        self.description = "Identifies the catchments, or sub-basin, draining to each link in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Subbasins' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/subbasins.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.subbasins(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TraceDownslopeFlowpaths(object):
    def __init__(self):
        self.label = "Trace Downslope Flowpaths"
        self.description = "Traces downslope flowpaths from one or more target sites (i.e. seed points). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#TraceDownslopeFlowpaths' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/trace_downslope_flowpaths.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.trace_downslope_flowpaths(seed_pts, d8_pntr=d8_pntr, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class UnnestBasins(object):
    def __init__(self):
        self.label = "Unnest Basins"
        self.description = "Extract whole watersheds for a set of outlet points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#UnnestBasins' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/unnest_basins.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.unnest_basins(d8_pntr, pour_pts=pour_pts, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Watershed(object):
    def __init__(self):
        self.label = "Watershed"
        self.description = "Identifies the watershed, or drainage basin, draining to a set of target cells. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html#Watershed' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/hydro_analysis/watershed.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.watershed(d8_pntr, pour_pts=pour_pts, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ChangeVectorAnalysis(object):
    def __init__(self):
        self.label = "Change Vector Analysis"
        self.description = "Performs a change vector analysis on a two-date multi-spectral dataset. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#ChangeVectorAnalysis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/change_vector_analysis.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        date1 = arcpy.Parameter(
            displayName="Earlier Date Input Files",
            name="date1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        date1.multiValue = True

        date2 = arcpy.Parameter(
            displayName="Later Date Input Files",
            name="date2",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        date2.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.change_vector_analysis(date1, date2=date2, magnitude=magnitude, direction=direction)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Closing(object):
    def __init__(self):
        self.label = "Closing"
        self.description = "A closing is a mathematical morphology operation involving an erosion (min filter) of a dilation (max filter) set. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#Closing' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/closing.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.closing(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CreateColourComposite(object):
    def __init__(self):
        self.label = "Create Colour Composite"
        self.description = "Creates a colour-composite image from three bands of multispectral imagery. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#CreateColourComposite' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/create_colour_composite.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.create_colour_composite(red, green=green, blue=blue, opacity=opacity, output=output, enhance=enhance, zeros=zeros)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FlipImage(object):
    def __init__(self):
        self.label = "Flip Image"
        self.description = "Reflects an image in the vertical or horizontal axis. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#FlipImage' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/flip_image.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.flip_image(input, output=output, direction=direction)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class IhsToRgb(object):
    def __init__(self):
        self.label = "Ihs To Rgb"
        self.description = "Converts intensity, hue, and saturation (IHS) images into red, green, and blue (RGB) images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#IhsToRgb' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/ihs_to_rgb.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.ihs_to_rgb(intensity, hue=hue, saturation=saturation, red=red, green=green, blue=blue, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ImageStackProfile(object):
    def __init__(self):
        self.label = "Image Stack Profile"
        self.description = "Plots an image stack profile (i.e. signature) for a set of points and multispectral images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#ImageStackProfile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/image_stack_profile.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.image_stack_profile(inputs, points=points, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class IntegralImage(object):
    def __init__(self):
        self.label = "Integral Image"
        self.description = "Transforms an input image (summed area table) into its integral image equivalent. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#IntegralImage' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/integral_image.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.integral_image(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class KMeansClustering(object):
    def __init__(self):
        self.label = "K Means Clustering"
        self.description = "Performs a k-means clustering operation on a multi-spectral dataset. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#KMeansClustering' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/k_means_clustering.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.k_means_clustering(inputs, output=output, out_html=out_html, classes=classes, max_iterations=max_iterations, class_change=class_change, initialize=initialize, min_class_size=min_class_size)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LineThinning(object):
    def __init__(self):
        self.label = "Line Thinning"
        self.description = "Performs line thinning a on Boolean raster image; intended to be used with the RemoveSpurs tool. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#LineThinning' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/line_thin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.line_thinning(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ModifiedKMeansClustering(object):
    def __init__(self):
        self.label = "Modified K Means Clustering"
        self.description = "Performs a modified k-means clustering operation on a multi-spectral dataset. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#ModifiedKMeansClustering' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/modified_k_means_clustering.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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

        merge_dist = arcpy.Parameter(
            displayName="Cluster Merger Distance",
            name="merge_dist",
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

        params = [inputs, output, out_html, start_clusters, merge_dist, max_iterations, class_change]

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
        merge_dist = parameters[4].valueAsText
        max_iterations = parameters[5].valueAsText
        class_change = parameters[6].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.modified_k_means_clustering(inputs, output=output, out_html=out_html, start_clusters=start_clusters, merge_dist=merge_dist, max_iterations=max_iterations, class_change=class_change)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Mosaic(object):
    def __init__(self):
        self.label = "Mosaic"
        self.description = "Mosaics two or more images together. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#Mosaic' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/mosaic.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.mosaic(inputs, output=output, method=method)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MosaicWithFeathering(object):
    def __init__(self):
        self.label = "Mosaic With Feathering"
        self.description = "Mosaics two images together using a feathering technique in overlapping areas to reduce edge-effects. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MosaicWithFeathering' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/mosaic_with_feathering.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.mosaic_with_feathering(input1, input2=input2, output=output, method=method, weight=weight)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NormalizedDifferenceIndex(object):
    def __init__(self):
        self.label = "Normalized Difference Index"
        self.description = "Calculate a normalized-difference index (NDI) from two bands of multispectral image data. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#NormalizedDifferenceIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/normalized_difference_index.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input1 = arcpy.Parameter(
            displayName="Input 1 File",
            name="input1",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        input2 = arcpy.Parameter(
            displayName="Input 2 File",
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

        clip = arcpy.Parameter(
            displayName="Distribution Tail Clip Amount (%)",
            name="clip",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip.value = "0.0"

        correction = arcpy.Parameter(
            displayName="Correction value",
            name="correction",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        correction.value = "0.0"

        params = [input1, input2, output, clip, correction]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input1 = parameters[0].valueAsText
        input2 = parameters[1].valueAsText
        output = parameters[2].valueAsText
        clip = parameters[3].valueAsText
        correction = parameters[4].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.normalized_difference_index(input1, input2=input2, output=output, clip=clip, correction=correction)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Opening(object):
    def __init__(self):
        self.label = "Opening"
        self.description = "An opening is a mathematical morphology operation involving a dilation (max filter) of an erosion (min filter) set. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#Opening' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/opening.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.opening(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RemoveSpurs(object):
    def __init__(self):
        self.label = "Remove Spurs"
        self.description = "Removes the spurs (pruning operation) from a Boolean line image; intended to be used on the output of the LineThinning tool. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#RemoveSpurs' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/remove_spurs.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.remove_spurs(input, output=output, iterations=iterations)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Resample(object):
    def __init__(self):
        self.label = "Resample"
        self.description = "Resamples one or more input images into a destination image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#Resample' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/resample.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.resample(inputs, destination=destination, method=method)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RgbToIhs(object):
    def __init__(self):
        self.label = "Rgb To Ihs"
        self.description = "Converts red, green, and blue (RGB) images into intensity, hue, and saturation (IHS) images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#RgbToIhs' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/rgb_to_ihs.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.rgb_to_ihs(red, green=green, blue=blue, composite=composite, intensity=intensity, hue=hue, saturation=saturation)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SplitColourComposite(object):
    def __init__(self):
        self.label = "Split Colour Composite"
        self.description = "This tool splits an RGB colour composite image into seperate multispectral images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#SplitColourComposite' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/split_colour_composite.rs' target='_blank'>GitHub</a>."
        self.category = "Image Processing Tools"

    def getParameterInfo(self):
        input = arcpy.Parameter(
            displayName="Input Colour Composite Image File",
            name="input",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        red = arcpy.Parameter(
            displayName="Output Red Band File",
            name="red",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        red.filter.list = ["tif"]

        green = arcpy.Parameter(
            displayName="Output Green Band File",
            name="green",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        green.filter.list = ["tif"]

        blue = arcpy.Parameter(
            displayName="Output Blue Band File",
            name="blue",
            datatype="DEFile",
            parameterType="Required",
            direction="Output")
        blue.filter.list = ["tif"]

        params = [input, red, green, blue]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        red = parameters[1].valueAsText
        green = parameters[2].valueAsText
        blue = parameters[3].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.split_colour_composite(input, red=red, green=green, blue=blue)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ThickenRasterLine(object):
    def __init__(self):
        self.label = "Thicken Raster Line"
        self.description = "Thickens single-cell wide lines within a raster image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#ThickenRasterLine' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/thicken_line.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.thicken_raster_line(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TophatTransform(object):
    def __init__(self):
        self.label = "Tophat Transform"
        self.description = "Performs either a white or black top-hat transform on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#TophatTransform' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/tophat.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.tophat_transform(input, output=output, filterx=filterx, filtery=filtery, variant=variant)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class WriteFunctionMemoryInsertion(object):
    def __init__(self):
        self.label = "Write Function Memory Insertion"
        self.description = "Performs a write function memory insertion for single-band multi-date change detection. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#WriteFunctionMemoryInsertion' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/write_func_memory_insertion.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.write_function_memory_insertion(input1, input2=input2, input3=input3, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AdaptiveFilter(object):
    def __init__(self):
        self.label = "Adaptive Filter"
        self.description = "Performs an adaptive filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#AdaptiveFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/adaptive_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.adaptive_filter(input, output=output, filterx=filterx, filtery=filtery, threshold=threshold)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BilateralFilter(object):
    def __init__(self):
        self.label = "Bilateral Filter"
        self.description = "A bilateral filter is an edge-preserving smoothing filter introduced by Tomasi and Manduchi (1998). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#BilateralFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/bilateral_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.bilateral_filter(input, output=output, sigma_dist=sigma_dist, sigma_int=sigma_int)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ConservativeSmoothingFilter(object):
    def __init__(self):
        self.label = "Conservative Smoothing Filter"
        self.description = "Performs a conservative-smoothing filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#ConservativeSmoothingFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/conservative_smoothing_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.conservative_smoothing_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CornerDetection(object):
    def __init__(self):
        self.label = "Corner Detection"
        self.description = "Identifies corner patterns in boolean images using hit-and-miss pattern matching. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#CornerDetection' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/corner_detection.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.corner_detection(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DiffOfGaussianFilter(object):
    def __init__(self):
        self.label = "Diff Of Gaussian Filter"
        self.description = "Performs a Difference of Gaussian (DoG) filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#DiffOfGaussianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/dog_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.diff_of_gaussian_filter(input, output=output, sigma1=sigma1, sigma2=sigma2)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DiversityFilter(object):
    def __init__(self):
        self.label = "Diversity Filter"
        self.description = "Assigns each cell in the output grid the number of different values in a moving window centred on each grid cell in the input raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#DiversityFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/diversity_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.diversity_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EdgePreservingMeanFilter(object):
    def __init__(self):
        self.label = "Edge Preserving Mean Filter"
        self.description = "Performs a simple edge-preserving mean filter on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#EdgePreservingMeanFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/edge_preserving_mean_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.edge_preserving_mean_filter(input, output=output, filter=filter, threshold=threshold)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EmbossFilter(object):
    def __init__(self):
        self.label = "Emboss Filter"
        self.description = "Performs an emboss filter on an image, similar to a hillshade operation. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#EmbossFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/emboss_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.emboss_filter(input, output=output, direction=direction, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FastAlmostGaussianFilter(object):
    def __init__(self):
        self.label = "Fast Almost Gaussian Filter"
        self.description = "Performs a fast approximate Gaussian filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#FastAlmostGaussianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/fast_almost_gaussian_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.fast_almost_gaussian_filter(input, output=output, sigma=sigma)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class GaussianFilter(object):
    def __init__(self):
        self.label = "Gaussian Filter"
        self.description = "Performs a Gaussian filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#GaussianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/gaussian_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.gaussian_filter(input, output=output, sigma=sigma)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HighPassFilter(object):
    def __init__(self):
        self.label = "High Pass Filter"
        self.description = "Performs a high-pass filter on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#HighPassFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/highpass_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.high_pass_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HighPassMedianFilter(object):
    def __init__(self):
        self.label = "High Pass Median Filter"
        self.description = "Performs a high pass median filter on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#HighPassMedianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/highpass_median_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.high_pass_median_filter(input, output=output, filterx=filterx, filtery=filtery, sig_digits=sig_digits)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class KNearestMeanFilter(object):
    def __init__(self):
        self.label = "K Nearest Mean Filter"
        self.description = "A k-nearest mean filter is a type of edge-preserving smoothing filter. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#KNearestMeanFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/k_nearest_mean_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.k_nearest_mean_filter(input, output=output, filterx=filterx, filtery=filtery, k=k)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LaplacianFilter(object):
    def __init__(self):
        self.label = "Laplacian Filter"
        self.description = "Performs a Laplacian filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#LaplacianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/laplacian_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.laplacian_filter(input, output=output, variant=variant, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LaplacianOfGaussianFilter(object):
    def __init__(self):
        self.label = "Laplacian Of Gaussian Filter"
        self.description = "Performs a Laplacian-of-Gaussian (LoG) filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#LaplacianOfGaussianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/log_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.laplacian_of_gaussian_filter(input, output=output, sigma=sigma)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LeeFilter(object):
    def __init__(self):
        self.label = "Lee Filter"
        self.description = "Performs a Lee (Sigma) smoothing filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#LeeFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/lee_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lee_filter(input, output=output, filterx=filterx, filtery=filtery, sigma=sigma, m=m)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LineDetectionFilter(object):
    def __init__(self):
        self.label = "Line Detection Filter"
        self.description = "Performs a line-detection filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#LineDetectionFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/line_detection_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.line_detection_filter(input, output=output, variant=variant, absvals=absvals, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MajorityFilter(object):
    def __init__(self):
        self.label = "Majority Filter"
        self.description = "Assigns each cell in the output grid the most frequently occurring value (mode) in a moving window centred on each grid cell in the input raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MajorityFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/majority_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.majority_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MaximumFilter(object):
    def __init__(self):
        self.label = "Maximum Filter"
        self.description = "Assigns each cell in the output grid the maximum value in a moving window centred on each grid cell in the input raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MaximumFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/max_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.maximum_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MeanFilter(object):
    def __init__(self):
        self.label = "Mean Filter"
        self.description = "Performs a mean filter (low-pass filter) on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MeanFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/mean_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.mean_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MedianFilter(object):
    def __init__(self):
        self.label = "Median Filter"
        self.description = "Performs a median filter on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MedianFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/median_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.median_filter(input, output=output, filterx=filterx, filtery=filtery, sig_digits=sig_digits)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinimumFilter(object):
    def __init__(self):
        self.label = "Minimum Filter"
        self.description = "Assigns each cell in the output grid the minimum value in a moving window centred on each grid cell in the input raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MinimumFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/min_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.minimum_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class OlympicFilter(object):
    def __init__(self):
        self.label = "Olympic Filter"
        self.description = "Performs an olympic smoothing filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#OlympicFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/olympic_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.olympic_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PercentileFilter(object):
    def __init__(self):
        self.label = "Percentile Filter"
        self.description = "Performs a percentile filter on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#PercentileFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/percentile_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.percentile_filter(input, output=output, filterx=filterx, filtery=filtery, sig_digits=sig_digits)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PrewittFilter(object):
    def __init__(self):
        self.label = "Prewitt Filter"
        self.description = "Performs a Prewitt edge-detection filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#PrewittFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/prewitt_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.prewitt_filter(input, output=output, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RangeFilter(object):
    def __init__(self):
        self.label = "Range Filter"
        self.description = "Assigns each cell in the output grid the range of values in a moving window centred on each grid cell in the input raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#RangeFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/range_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.range_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RobertsCrossFilter(object):
    def __init__(self):
        self.label = "Roberts Cross Filter"
        self.description = "Performs a Robert's cross edge-detection filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#RobertsCrossFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/roberts_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.roberts_cross_filter(input, output=output, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ScharrFilter(object):
    def __init__(self):
        self.label = "Scharr Filter"
        self.description = "Performs a Scharr edge-detection filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#ScharrFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/scharr_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.scharr_filter(input, output=output, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SobelFilter(object):
    def __init__(self):
        self.label = "Sobel Filter"
        self.description = "Performs a Sobel edge-detection filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#SobelFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/sobel_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sobel_filter(input, output=output, variant=variant, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StandardDeviationFilter(object):
    def __init__(self):
        self.label = "Standard Deviation Filter"
        self.description = "Assigns each cell in the output grid the standard deviation of values in a moving window centred on each grid cell in the input raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#StandardDeviationFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/stdev_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.standard_deviation_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TotalFilter(object):
    def __init__(self):
        self.label = "Total Filter"
        self.description = "Performs a total filter on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#TotalFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/total_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.total_filter(input, output=output, filterx=filterx, filtery=filtery)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class UnsharpMasking(object):
    def __init__(self):
        self.label = "Unsharp Masking"
        self.description = "An image sharpening technique that enhances edges. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#UnsharpMasking' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/unsharp_masking.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.unsharp_masking(input, output=output, sigma=sigma, amount=amount, threshold=threshold)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class UserDefinedWeightsFilter(object):
    def __init__(self):
        self.label = "User Defined Weights Filter"
        self.description = "Performs a user-defined weights filter on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#UserDefinedWeightsFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/user_defined_weights_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.user_defined_weights_filter(input, weights=weights, output=output, center=center, normalize=normalize)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class BalanceContrastEnhancement(object):
    def __init__(self):
        self.label = "Balance Contrast Enhancement"
        self.description = "Performs a balance contrast enhancement on a colour-composite image of multispectral data. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#BalanceContrastEnhancement' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/balance_contrast_enhancement.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.balance_contrast_enhancement(input, output=output, band_mean=band_mean)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CorrectVignetting(object):
    def __init__(self):
        self.label = "Correct Vignetting"
        self.description = "Corrects the darkening of images towards corners. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#CorrectVignetting' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/correct_vignetting.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.correct_vignetting(input, pp=pp, output=output, focal_length=focal_length, image_width=image_width, n=n)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DirectDecorrelationStretch(object):
    def __init__(self):
        self.label = "Direct Decorrelation Stretch"
        self.description = "Performs a direct decorrelation stretch enhancement on a colour-composite image of multispectral data. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#DirectDecorrelationStretch' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/direct_decorrelation_stretch.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.direct_decorrelation_stretch(input, output=output, k=k, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class GammaCorrection(object):
    def __init__(self):
        self.label = "Gamma Correction"
        self.description = "Performs a gamma correction on an input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#GammaCorrection' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/gamma_correction.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.gamma_correction(input, output=output, gamma=gamma)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class GaussianContrastStretch(object):
    def __init__(self):
        self.label = "Gaussian Contrast Stretch"
        self.description = "Performs a Gaussian contrast stretch on input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#GaussianContrastStretch' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/gaussian_contrast_stretch.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.gaussian_contrast_stretch(input, output=output, num_tones=num_tones)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HistogramEqualization(object):
    def __init__(self):
        self.label = "Histogram Equalization"
        self.description = "Performs a histogram equalization contrast enhancment on an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#HistogramEqualization' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/histogram_equalization.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.histogram_equalization(input, output=output, num_tones=num_tones)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HistogramMatching(object):
    def __init__(self):
        self.label = "Histogram Matching"
        self.description = "Alters the statistical distribution of a raster image matching it to a specified PDF. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#HistogramMatching' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/histogram_matching.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.histogram_matching(input, histo_file=histo_file, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HistogramMatchingTwoImages(object):
    def __init__(self):
        self.label = "Histogram Matching Two Images"
        self.description = "This tool alters the cumulative distribution function of a raster image to that of another image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#HistogramMatchingTwoImages' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/histogram_matching_two_images.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.histogram_matching_two_images(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class MinMaxContrastStretch(object):
    def __init__(self):
        self.label = "Min Max Contrast Stretch"
        self.description = "Performs a min-max contrast stretch on an input greytone image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#MinMaxContrastStretch' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/min_max_contrast_stretch.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.min_max_contrast_stretch(input, output=output, min_val=min_val, max_val=max_val, num_tones=num_tones)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PanchromaticSharpening(object):
    def __init__(self):
        self.label = "Panchromatic Sharpening"
        self.description = "Increases the spatial resolution of image data by combining multispectral bands with panchromatic data. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#PanchromaticSharpening' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/pan_sharpening.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.panchromatic_sharpening(red, green=green, blue=blue, composite=composite, pan=pan, output=output, method=method)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PercentageContrastStretch(object):
    def __init__(self):
        self.label = "Percentage Contrast Stretch"
        self.description = "Performs a percentage linear contrast stretch on input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#PercentageContrastStretch' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/percentage_contrast_stretch.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.percentage_contrast_stretch(input, output=output, clip=clip, tail=tail, num_tones=num_tones)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SigmoidalContrastStretch(object):
    def __init__(self):
        self.label = "Sigmoidal Contrast Stretch"
        self.description = "Performs a sigmoidal contrast stretch on input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#SigmoidalContrastStretch' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/sigmoidal_contrast_stretch.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sigmoidal_contrast_stretch(input, output=output, cutoff=cutoff, gain=gain, num_tones=num_tones)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StandardDeviationContrastStretch(object):
    def __init__(self):
        self.label = "Standard Deviation Contrast Stretch"
        self.description = "Performs a standard-deviation contrast stretch on input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/image_processing_tools.html#StandardDeviationContrastStretch' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/image_analysis/stdev_contrast_stretch.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.standard_deviation_contrast_stretch(input, output=output, stdev=stdev, num_tones=num_tones)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ClassifyOverlapPoints(object):
    def __init__(self):
        self.label = "Classify Overlap Points"
        self.description = "Classifies or filters LAS points in regions of overlapping flight lines. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#ClassifyOverlapPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/classify_overlap_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.classify_overlap_points(input, output=output, resolution=resolution, filter=filter)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ClipLidarToPolygon(object):
    def __init__(self):
        self.label = "Clip Lidar To Polygon"
        self.description = "Clips a LiDAR point cloud to a vector polygon or polygons. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#ClipLidarToPolygon' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/clip_lidar_to_polygon.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.clip_lidar_to_polygon(input, polygons=polygons, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ErasePolygonFromLidar(object):
    def __init__(self):
        self.label = "Erase Polygon From Lidar"
        self.description = "Erases (cuts out) a vector polygon or polygons from a LiDAR point cloud. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#ErasePolygonFromLidar' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/erase_polygon_from_lidar.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.erase_polygon_from_lidar(input, polygons=polygons, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FilterLidarScanAngles(object):
    def __init__(self):
        self.label = "Filter Lidar Scan Angles"
        self.description = "Removes points in a LAS file with scan angles greater than a threshold. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#FilterLidarScanAngles' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/filter_lidar_scan_angles.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.filter_lidar_scan_angles(input, output=output, threshold=threshold)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindFlightlineEdgePoints(object):
    def __init__(self):
        self.label = "Find Flightline Edge Points"
        self.description = "Identifies points along a flightline's edge in a LAS file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#FindFlightlineEdgePoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/find_flightline_edge_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_flightline_edge_points(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FlightlineOverlap(object):
    def __init__(self):
        self.label = "Flightline Overlap"
        self.description = "Reads a LiDAR (LAS) point file and outputs a raster containing the number of overlapping flight lines in each grid cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#FlightlineOverlap' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/flightline_overlap.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.flightline_overlap(input, output=output, resolution=resolution)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LasToAscii(object):
    def __init__(self):
        self.label = "Las To Ascii"
        self.description = "Converts one or more LAS files into ASCII text files. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LasToAscii' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/las_to_ascii.rs' target='_blank'>GitHub</a>."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input LiDAR Files",
            name="inputs",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True
        inputs.filter.list = ["las", "zip"]

        params = [inputs]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        inputs = parameters[0].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.las_to_ascii(inputs)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LasToMultipointShapefile(object):
    def __init__(self):
        self.label = "Las To Multipoint Shapefile"
        self.description = "Converts one or more LAS files into MultipointZ vector Shapefiles. When the input parameter is not specified, the tool grids all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LasToMultipointShapefile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/las_to_multipoint_shapefile.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.las_to_multipoint_shapefile(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LasToShapefile(object):
    def __init__(self):
        self.label = "Las To Shapefile"
        self.description = "Converts one or more LAS files into a vector Shapefile of POINT ShapeType. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LasToShapefile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/las_to_shapefile.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.las_to_shapefile(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarBlockMaximum(object):
    def __init__(self):
        self.label = "Lidar Block Maximum"
        self.description = "Creates a block-maximum raster from an input LAS file. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarBlockMaximum' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/block_maximum.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_block_maximum(input, output=output, resolution=resolution)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarBlockMinimum(object):
    def __init__(self):
        self.label = "Lidar Block Minimum"
        self.description = "Creates a block-minimum raster from an input LAS file. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarBlockMinimum' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/block_minimum.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_block_minimum(input, output=output, resolution=resolution)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarClassifySubset(object):
    def __init__(self):
        self.label = "Lidar Classify Subset"
        self.description = "Classifies the values in one LiDAR point cloud that correpond with points in a subset cloud. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarClassifySubset' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_classify_subset.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_classify_subset(base, subset=subset, output=output, subset_class=subset_class, nonsubset_class=nonsubset_class)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarColourize(object):
    def __init__(self):
        self.label = "Lidar Colourize"
        self.description = "Adds the red-green-blue colour fields of a LiDAR (LAS) file based on an input image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarColourize' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_colourize.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_colourize(in_lidar, in_image=in_image, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarConstructVectorTin(object):
    def __init__(self):
        self.label = "Lidar Construct Vector Tin"
        self.description = "Creates a vector triangular irregular network (TIN) fitted to LiDAR points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarConstructVectorTin' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_construct_vector_tin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_construct_vector_tin(input, output=output, returns=returns, exclude_cls=exclude_cls, minz=minz, maxz=maxz)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarElevationSlice(object):
    def __init__(self):
        self.label = "Lidar Elevation Slice"
        self.description = "Outputs all of the points within a LiDAR (LAS) point file that lie between a specified elevation range. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarElevationSlice' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_elevation_slice.rs' target='_blank'>GitHub</a>."
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

        cls = arcpy.Parameter(
            displayName="Retain but reclass points outside the specified elevation range?",
            name="cls",
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

        params = [input, output, minz, maxz, cls, inclassval, outclassval]

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
        cls = parameters[4].valueAsText
        inclassval = parameters[5].valueAsText
        outclassval = parameters[6].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_elevation_slice(input, output=output, minz=minz, maxz=maxz, cls=cls, inclassval=inclassval, outclassval=outclassval)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarGroundPointFilter(object):
    def __init__(self):
        self.label = "Lidar Ground Point Filter"
        self.description = "Identifies ground points within LiDAR dataset using a slope-based method. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarGroundPointFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_ground_point_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_ground_point_filter(input, output=output, radius=radius, min_neighbours=min_neighbours, slope_threshold=slope_threshold, height_threshold=height_threshold, classify=classify, slope_norm=slope_norm)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarHexBinning(object):
    def __init__(self):
        self.label = "Lidar Hex Binning"
        self.description = "Hex-bins a set of LiDAR points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarHexBinning' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_hex_bin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_hex_binning(input, output=output, width=width, orientation=orientation)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarHillshade(object):
    def __init__(self):
        self.label = "Lidar Hillshade"
        self.description = "Calculates a hillshade value for points within a LAS file and stores these data in the RGB field. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarHillshade' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_hillshade.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_hillshade(input, output=output, azimuth=azimuth, altitude=altitude, radius=radius)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarHistogram(object):
    def __init__(self):
        self.label = "Lidar Histogram"
        self.description = "Creates a histogram of LiDAR data. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarHistogram' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_histogram.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_histogram(input, output=output, parameter=parameter, clip=clip)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarIdwInterpolation(object):
    def __init__(self):
        self.label = "Lidar Idw Interpolation"
        self.description = "Interpolates LAS files using an inverse-distance weighted (IDW) scheme. When the input/output parameters are not specified, the tool interpolates all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarIdwInterpolation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_idw_interpolation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_idw_interpolation(input, output=output, parameter=parameter, returns=returns, resolution=resolution, weight=weight, radius=radius, exclude_cls=exclude_cls, minz=minz, maxz=maxz)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarInfo(object):
    def __init__(self):
        self.label = "Lidar Info"
        self.description = "Prints information about a LiDAR (LAS) dataset, including header, point return frequency, and classification data and information about the variable length records (VLRs) and geokeys. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarInfo' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_info.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_info(input, output=output, vlr=vlr, geokeys=geokeys)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarJoin(object):
    def __init__(self):
        self.label = "Lidar Join"
        self.description = "Joins multiple LiDAR (LAS) files into a single LAS file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarJoin' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_join.rs' target='_blank'>GitHub</a>."
        self.category = "LiDAR Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input LiDAR Files",
            name="inputs",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_join(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarKappaIndex(object):
    def __init__(self):
        self.label = "Lidar Kappa Index"
        self.description = "Performs a kappa index of agreement (KIA) analysis on the classifications of two LAS files. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarKappaIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_kappa.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_kappa_index(input1, input2=input2, output=output, class_accuracy=class_accuracy, resolution=resolution)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarNearestNeighbourGridding(object):
    def __init__(self):
        self.label = "Lidar Nearest Neighbour Gridding"
        self.description = "Grids LAS files using nearest-neighbour scheme. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarNearestNeighbourGridding' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_nn_gridding.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_nearest_neighbour_gridding(input, output=output, parameter=parameter, returns=returns, resolution=resolution, radius=radius, exclude_cls=exclude_cls, minz=minz, maxz=maxz)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarPointDensity(object):
    def __init__(self):
        self.label = "Lidar Point Density"
        self.description = "Calculates the spatial pattern of point density for a LiDAR data set. When the input/output parameters are not specified, the tool grids all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarPointDensity' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_point_density.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_point_density(input, output=output, returns=returns, resolution=resolution, radius=radius, exclude_cls=exclude_cls, minz=minz, maxz=maxz)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarPointStats(object):
    def __init__(self):
        self.label = "Lidar Point Stats"
        self.description = "Creates several rasters summarizing the distribution of LAS point data. When the input/output parameters are not specified, the tool works on all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarPointStats' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_point_stats.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_point_stats(input, resolution=resolution, num_points=num_points, num_pulses=num_pulses, z_range=z_range, intensity_range=intensity_range, predom_class=predom_class)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarRemoveDuplicates(object):
    def __init__(self):
        self.label = "Lidar Remove Duplicates"
        self.description = "Removes duplicate points from a LiDAR data set. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarRemoveDuplicates' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/remove_duplicates.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_remove_duplicates(input, output=output, include_z=include_z)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarRemoveOutliers(object):
    def __init__(self):
        self.label = "Lidar Remove Outliers"
        self.description = "Removes outliers (high and low points) in a LiDAR point cloud. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarRemoveOutliers' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_outliers.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_remove_outliers(input, output=output, radius=radius, elev_diff=elev_diff)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarSegmentation(object):
    def __init__(self):
        self.label = "Lidar Segmentation"
        self.description = "Segments a LiDAR point cloud based on normal vectors. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarSegmentation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_segmentation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_segmentation(input, output=output, radius=radius, norm_diff=norm_diff, maxzdiff=maxzdiff)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarSegmentationBasedFilter(object):
    def __init__(self):
        self.label = "Lidar Segmentation Based Filter"
        self.description = "Identifies ground points within LiDAR point clouds using a segmentation based approach. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarSegmentationBasedFilter' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_segmentation_based_filter.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_segmentation_based_filter(input, output=output, radius=radius, norm_diff=norm_diff, maxzdiff=maxzdiff, classify=classify)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarThin(object):
    def __init__(self):
        self.label = "Lidar Thin"
        self.description = "Thins a LiDAR point cloud, reducing point density. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarThin' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_thin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_thin(input, output=output, resolution=resolution, method=method, save_filtered=save_filtered)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarThinHighDensity(object):
    def __init__(self):
        self.label = "Lidar Thin High Density"
        self.description = "Thins points from high density areas within a LiDAR point cloud. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarThinHighDensity' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_thin_high_density.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_thin_high_density(input, output=output, resolution=resolution, density=density, save_filtered=save_filtered)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarTile(object):
    def __init__(self):
        self.label = "Lidar Tile"
        self.description = "Tiles a LiDAR LAS file into multiple LAS files. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarTile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_tile.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_tile(input, width=width, height=height, origin_x=origin_x, origin_y=origin_y, min_points=min_points)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarTileFootprint(object):
    def __init__(self):
        self.label = "Lidar Tile Footprint"
        self.description = "Creates a vector polygon of the convex hull of a LiDAR point cloud. When the input/output parameters are not specified, the tool works with all LAS files contained within the working directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarTileFootprint' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_tile_footprint.rs' target='_blank'>GitHub</a>."
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

        hull = arcpy.Parameter(
            displayName="Create Convex Hull Around Points",
            name="hull",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        hull.value = "false"

        params = [input, output, hull]

        return params

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        input = parameters[0].valueAsText
        output = parameters[1].valueAsText
        hull = parameters[2].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_tile_footprint(input, output=output, hull=hull)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarTinGridding(object):
    def __init__(self):
        self.label = "Lidar Tin Gridding"
        self.description = "Creates a raster grid based on a Delaunay triangular irregular network (TIN) fitted to LiDAR points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarTinGridding' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_tin_gridding.rs' target='_blank'>GitHub</a>."
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

        max_triangle_edge_length = arcpy.Parameter(
            displayName="Maximum Triangle Edge Length",
            name="max_triangle_edge_length",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        params = [input, output, parameter, returns, resolution, exclude_cls, minz, maxz, max_triangle_edge_length]

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
        max_triangle_edge_length = parameters[8].valueAsText
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_tin_gridding(input, output=output, parameter=parameter, returns=returns, resolution=resolution, exclude_cls=exclude_cls, minz=minz, maxz=maxz, max_triangle_edge_length=max_triangle_edge_length)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LidarTophatTransform(object):
    def __init__(self):
        self.label = "Lidar Tophat Transform"
        self.description = "Performs a white top-hat transform on a Lidar dataset; as an estimate of height above ground, this is useful for modelling the vegetation canopy. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#LidarTophatTransform' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/lidar_tophat_transform.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.lidar_tophat_transform(input, output=output, radius=radius)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NormalVectors(object):
    def __init__(self):
        self.label = "Normal Vectors"
        self.description = "Calculates normal vectors for points within a LAS file and stores these data (XYZ vector components) in the RGB field. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#NormalVectors' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/normal_vectors.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.normal_vectors(input, output=output, radius=radius)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SelectTilesByPolygon(object):
    def __init__(self):
        self.label = "Select Tiles By Polygon"
        self.description = "Copies LiDAR tiles overlapping with a polygon into an output directory. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html#SelectTilesByPolygon' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/lidar_analysis/select_tiles_by_polygon.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.select_tiles_by_polygon(indir, outdir=outdir, polygons=polygons)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class And(object):
    def __init__(self):
        self.label = "And"
        self.description = "Performs a logical AND operator on two Boolean raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#And' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/and.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.And(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Not(object):
    def __init__(self):
        self.label = "Not"
        self.description = "Performs a logical NOT operator on two Boolean raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Not' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/not.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.Not(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Or(object):
    def __init__(self):
        self.label = "Or"
        self.description = "Performs a logical OR operator on two Boolean raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Or' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/or.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.Or(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AbsoluteValue(object):
    def __init__(self):
        self.label = "Absolute Value"
        self.description = "Calculates the absolute value of every cell in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#AbsoluteValue' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/abs.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.absolute_value(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Add(object):
    def __init__(self):
        self.label = "Add"
        self.description = "Performs an addition operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Add' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/add.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.add(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Anova(object):
    def __init__(self):
        self.label = "Anova"
        self.description = "Performs an analysis of variance (ANOVA) test on a raster dataset. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Anova' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/anova.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.anova(input, features=features, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ArcCos(object):
    def __init__(self):
        self.label = "Arc Cos"
        self.description = "Returns the inverse cosine (arccos) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ArcCos' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/arccos.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.arc_cos(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ArcSin(object):
    def __init__(self):
        self.label = "Arc Sin"
        self.description = "Returns the inverse sine (arcsin) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ArcSin' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/arcsin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.arc_sin(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ArcTan(object):
    def __init__(self):
        self.label = "Arc Tan"
        self.description = "Returns the inverse tangent (arctan) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ArcTan' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/arctan.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.arc_tan(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Atan2(object):
    def __init__(self):
        self.label = "Atan2"
        self.description = "Returns the 2-argument inverse tangent (atan2). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Atan2' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/atan2.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.atan2(input_y, input_x=input_x, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AttributeCorrelation(object):
    def __init__(self):
        self.label = "Attribute Correlation"
        self.description = "Performs a correlation analysis on attribute fields from a vector database. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#AttributeCorrelation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/attribute_correlation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.attribute_correlation(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AttributeHistogram(object):
    def __init__(self):
        self.label = "Attribute Histogram"
        self.description = "Creates a histogram for the field values of a vector's attribute table. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#AttributeHistogram' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/attribute_histogram.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.attribute_histogram(input, field=field, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class AttributeScattergram(object):
    def __init__(self):
        self.label = "Attribute Scattergram"
        self.description = "Creates a scattergram for two field values of a vector's attribute table. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#AttributeScattergram' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/attribute_scattergram.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.attribute_scattergram(input, fieldx=fieldx, fieldy=fieldy, output=output, trendline=trendline)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Ceil(object):
    def __init__(self):
        self.label = "Ceil"
        self.description = "Returns the smallest (closest to negative infinity) value that is greater than or equal to the values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Ceil' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/ceil.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.ceil(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Cos(object):
    def __init__(self):
        self.label = "Cos"
        self.description = "Returns the cosine (cos) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Cos' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/cos.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cos(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Cosh(object):
    def __init__(self):
        self.label = "Cosh"
        self.description = "Returns the hyperbolic cosine (cosh) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Cosh' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/cosh.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cosh(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CrispnessIndex(object):
    def __init__(self):
        self.label = "Crispness Index"
        self.description = "Calculates the Crispness Index, which is used to quantify how crisp (or conversely how fuzzy) a probability image is. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#CrispnessIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/crispness_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.crispness_index(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CrossTabulation(object):
    def __init__(self):
        self.label = "Cross Tabulation"
        self.description = "Performs a cross-tabulation on two categorical images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#CrossTabulation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/cross_tabulation.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cross_tabulation(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class CumulativeDistribution(object):
    def __init__(self):
        self.label = "Cumulative Distribution"
        self.description = "Converts a raster image to its cumulative distribution function. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#CumulativeDistribution' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/cumulative_dist.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.cumulative_distribution(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Decrement(object):
    def __init__(self):
        self.label = "Decrement"
        self.description = "Decreases the values of each grid cell in an input raster by 1.0 (see also InPlaceSubtract). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Decrement' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/decrement.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.decrement(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Divide(object):
    def __init__(self):
        self.label = "Divide"
        self.description = "Performs a division operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Divide' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/divide.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.divide(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class EqualTo(object):
    def __init__(self):
        self.label = "Equal To"
        self.description = "Performs a equal-to comparison operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#EqualTo' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/equal_to.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.equal_to(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Exp(object):
    def __init__(self):
        self.label = "Exp"
        self.description = "Returns the exponential (base e) of values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Exp' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/exp.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.exp(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Exp2(object):
    def __init__(self):
        self.label = "Exp2"
        self.description = "Returns the exponential (base 2) of values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Exp2' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/exp2.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.exp2(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExtractRasterStatistics(object):
    def __init__(self):
        self.label = "Extract Raster Statistics"
        self.description = "Extracts descriptive statistics for a group of patches in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ExtractRasterStatistics' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/extract_statistics.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.extract_raster_statistics(input, features=features, output=output, stat=stat, out_table=out_table)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Floor(object):
    def __init__(self):
        self.label = "Floor"
        self.description = "Returns the largest (closest to positive infinity) value that is less than or equal to the values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Floor' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/floor.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.floor(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class GreaterThan(object):
    def __init__(self):
        self.label = "Greater Than"
        self.description = "Performs a greater-than comparison operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#GreaterThan' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/greater_than.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.greater_than(input1, input2=input2, output=output, incl_equals=incl_equals)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ImageAutocorrelation(object):
    def __init__(self):
        self.label = "Image Autocorrelation"
        self.description = "Performs Moran's I analysis on two or more input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ImageAutocorrelation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/image_autocorrelation.rs' target='_blank'>GitHub</a>."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.image_autocorrelation(inputs, contiguity=contiguity, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ImageCorrelation(object):
    def __init__(self):
        self.label = "Image Correlation"
        self.description = "Performs image correlation on two or more input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ImageCorrelation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/image_correlation.rs' target='_blank'>GitHub</a>."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.image_correlation(inputs, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ImageRegression(object):
    def __init__(self):
        self.label = "Image Regression"
        self.description = "Performs image regression analysis on two input images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ImageRegression' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/image_regression.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.image_regression(input1, input2=input2, output=output, out_residuals=out_residuals, standardize=standardize)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class InPlaceAdd(object):
    def __init__(self):
        self.label = "In Place Add"
        self.description = "Performs an in-place addition operation (input1 += input2). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#InPlaceAdd' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/inplace_add.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.in_place_add(input1, input2=input2)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class InPlaceDivide(object):
    def __init__(self):
        self.label = "In Place Divide"
        self.description = "Performs an in-place division operation (input1 /= input2). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#InPlaceDivide' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/inplace_divide.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.in_place_divide(input1, input2=input2)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class InPlaceMultiply(object):
    def __init__(self):
        self.label = "In Place Multiply"
        self.description = "Performs an in-place multiplication operation (input1 *= input2). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#InPlaceMultiply' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/inplace_multiply.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.in_place_multiply(input1, input2=input2)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class InPlaceSubtract(object):
    def __init__(self):
        self.label = "In Place Subtract"
        self.description = "Performs an in-place subtraction operation (input1 -= input2). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#InPlaceSubtract' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/inplace_subtract.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.in_place_subtract(input1, input2=input2)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Increment(object):
    def __init__(self):
        self.label = "Increment"
        self.description = "Increases the values of each grid cell in an input raster by 1.0. (see also InPlaceAdd). View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Increment' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/increment.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.increment(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class IntegerDivision(object):
    def __init__(self):
        self.label = "Integer Division"
        self.description = "Performs an integer division operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#IntegerDivision' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/integer_division.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.integer_division(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class IsNoData(object):
    def __init__(self):
        self.label = "Is No Data"
        self.description = "Identifies NoData valued pixels in an image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#IsNoData' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/isnodata.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.is_no_data(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class KappaIndex(object):
    def __init__(self):
        self.label = "Kappa Index"
        self.description = "Performs a kappa index of agreement (KIA) analysis on two categorical raster files. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#KappaIndex' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/kappa_index.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.kappa_index(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class KsTestForNormality(object):
    def __init__(self):
        self.label = "Ks Test For Normality"
        self.description = "Evaluates whether the values in a raster are normally distributed. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#KsTestForNormality' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/ks_normality_test.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.ks_test_for_normality(input, output=output, num_samples=num_samples)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LessThan(object):
    def __init__(self):
        self.label = "Less Than"
        self.description = "Performs a less-than comparison operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#LessThan' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/less_than.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.less_than(input1, input2=input2, output=output, incl_equals=incl_equals)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ListUniqueValues(object):
    def __init__(self):
        self.label = "List Unique Values"
        self.description = "Lists the unique values contained in a field witin a vector's attribute table. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ListUniqueValues' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/list_unique_values.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.list_unique_values(input, field=field, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Ln(object):
    def __init__(self):
        self.label = "Ln"
        self.description = "Returns the natural logarithm of values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Ln' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/ln.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.ln(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Log10(object):
    def __init__(self):
        self.label = "Log10"
        self.description = "Returns the base-10 logarithm of values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Log10' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/log10.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.log10(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Log2(object):
    def __init__(self):
        self.label = "Log2"
        self.description = "Returns the base-2 logarithm of values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Log2' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/log2.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.log2(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Max(object):
    def __init__(self):
        self.label = "Max"
        self.description = "Performs a MAX operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Max' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/max.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.max(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Min(object):
    def __init__(self):
        self.label = "Min"
        self.description = "Performs a MIN operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Min' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/min.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.min(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Modulo(object):
    def __init__(self):
        self.label = "Modulo"
        self.description = "Performs a modulo operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Modulo' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/modulo.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.modulo(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Multiply(object):
    def __init__(self):
        self.label = "Multiply"
        self.description = "Performs a multiplication operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Multiply' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/multiply.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.multiply(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Negate(object):
    def __init__(self):
        self.label = "Negate"
        self.description = "Changes the sign of values in a raster or the 0-1 values of a Boolean raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Negate' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/negate.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.negate(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class NotEqualTo(object):
    def __init__(self):
        self.label = "Not Equal To"
        self.description = "Performs a not-equal-to comparison operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#NotEqualTo' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/not_equal_to.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.not_equal_to(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Power(object):
    def __init__(self):
        self.label = "Power"
        self.description = "Raises the values in grid cells of one rasters, or a constant value, by values in another raster or constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Power' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/power.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.power(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class PrincipalComponentAnalysis(object):
    def __init__(self):
        self.label = "Principal Component Analysis"
        self.description = "Performs a principal component analysis (PCA) on a multi-spectral dataset. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#PrincipalComponentAnalysis' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/principal_component_analysis.rs' target='_blank'>GitHub</a>."
        self.category = "Math and Stats Tools"

    def getParameterInfo(self):
        inputs = arcpy.Parameter(
            displayName="Input Files",
            name="inputs",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")
        inputs.multiValue = True

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.principal_component_analysis(inputs, output=output, num_comp=num_comp, standardized=standardized)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Quantiles(object):
    def __init__(self):
        self.label = "Quantiles"
        self.description = "Transforms raster values into quantiles. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Quantiles' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/quantiles.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.quantiles(input, output=output, num_quantiles=num_quantiles)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RandomField(object):
    def __init__(self):
        self.label = "Random Field"
        self.description = "Creates an image containing random values. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#RandomField' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/random_field.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.random_field(base, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RandomSample(object):
    def __init__(self):
        self.label = "Random Sample"
        self.description = "Creates an image containing randomly located sample grid cells with unique IDs. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#RandomSample' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/random_sample.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.random_sample(base, output=output, num_samples=num_samples)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterHistogram(object):
    def __init__(self):
        self.label = "Raster Histogram"
        self.description = "Creates a histogram from raster values. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#RasterHistogram' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/raster_histogram.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_histogram(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterSummaryStats(object):
    def __init__(self):
        self.label = "Raster Summary Stats"
        self.description = "Measures a rasters min, max, average, standard deviation, num. non-nodata cells, and total. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#RasterSummaryStats' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/raster_summary_stats.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_summary_stats(input)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Reciprocal(object):
    def __init__(self):
        self.label = "Reciprocal"
        self.description = "Returns the reciprocal (i.e. 1 / z) of values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Reciprocal' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/reciprocal.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.reciprocal(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RescaleValueRange(object):
    def __init__(self):
        self.label = "Rescale Value Range"
        self.description = "Performs a min-max contrast stretch on an input greytone image. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#RescaleValueRange' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/rescale_value_range.rs' target='_blank'>GitHub</a>."
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
            displayName="Lower-Tail Clip Value",
            name="clip_min",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        clip_max = arcpy.Parameter(
            displayName="Upper-Tail Clip Value",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.rescale_value_range(input, output=output, out_min_val=out_min_val, out_max_val=out_max_val, clip_min=clip_min, clip_max=clip_max)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RootMeanSquareError(object):
    def __init__(self):
        self.label = "Root Mean Square Error"
        self.description = "Calculates the RMSE and other accuracy statistics. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#RootMeanSquareError' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/root_mean_square_error.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.root_mean_square_error(input, base=base)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Round(object):
    def __init__(self):
        self.label = "Round"
        self.description = "Rounds the values in an input raster to the nearest integer value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Round' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/round.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.round(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Sin(object):
    def __init__(self):
        self.label = "Sin"
        self.description = "Returns the sine (sin) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Sin' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/sin.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sin(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Sinh(object):
    def __init__(self):
        self.label = "Sinh"
        self.description = "Returns the hyperbolic sine (sinh) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Sinh' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/sinh.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.sinh(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Square(object):
    def __init__(self):
        self.label = "Square"
        self.description = "Squares the values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Square' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/square.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.square(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class SquareRoot(object):
    def __init__(self):
        self.label = "Square Root"
        self.description = "Returns the square root of the values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#SquareRoot' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/sqrt.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.square_root(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Subtract(object):
    def __init__(self):
        self.label = "Subtract"
        self.description = "Performs a differencing operation on two rasters or a raster and a constant value. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Subtract' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/subtract.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.subtract(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Tan(object):
    def __init__(self):
        self.label = "Tan"
        self.description = "Returns the tangent (tan) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Tan' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/tan.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.tan(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Tanh(object):
    def __init__(self):
        self.label = "Tanh"
        self.description = "Returns the hyperbolic tangent (tanh) of each values in a raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Tanh' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/tanh.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.tanh(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ToDegrees(object):
    def __init__(self):
        self.label = "To Degrees"
        self.description = "Converts a raster from radians to degrees. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ToDegrees' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/to_degrees.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.to_degrees(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ToRadians(object):
    def __init__(self):
        self.label = "To Radians"
        self.description = "Converts a raster from degrees to radians. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ToRadians' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/to_radians.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.to_radians(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TrendSurface(object):
    def __init__(self):
        self.label = "Trend Surface"
        self.description = "Estimates the trend surface of an input raster file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#TrendSurface' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/trend_surface.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.trend_surface(input, output=output, order=order)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TrendSurfaceVectorPoints(object):
    def __init__(self):
        self.label = "Trend Surface Vector Points"
        self.description = "Estimates a trend surface from vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#TrendSurfaceVectorPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/trend_surface_vector_points.rs' target='_blank'>GitHub</a>."
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
            displayName="Cell Size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Required",
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.trend_surface_vector_points(input, field=field, output=output, order=order, cell_size=cell_size)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Truncate(object):
    def __init__(self):
        self.label = "Truncate"
        self.description = "Truncates the values in a raster to the desired number of decimal places. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Truncate' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/truncate.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.truncate(input, output=output, num_decimals=num_decimals)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TurningBandsSimulation(object):
    def __init__(self):
        self.label = "Turning Bands Simulation"
        self.description = "Creates an image containing random values based on a turning-bands simulation. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#TurningBandsSimulation' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/turning_bands.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.turning_bands_simulation(base, output=output, range=range, iterations=iterations)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class Xor(object):
    def __init__(self):
        self.label = "Xor"
        self.description = "Performs a logical XOR operator on two Boolean raster images. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#Xor' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/xor.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.xor(input1, input2=input2, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ZScores(object):
    def __init__(self):
        self.label = "Z Scores"
        self.description = "Standardizes the values in an input raster by converting to z-scores. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html#ZScores' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/math_stat_analysis/zscores.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.z_scores(input, output=output)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class DistanceToOutlet(object):
    def __init__(self):
        self.label = "Distance To Outlet"
        self.description = "Calculates the distance of stream grid cells to the channel network outlet cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#DistanceToOutlet' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/dist_to_outlet.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.distance_to_outlet(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExtractStreams(object):
    def __init__(self):
        self.label = "Extract Streams"
        self.description = "Extracts stream grid cells from a flow accumulation raster. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#ExtractStreams' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/extract_streams.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.extract_streams(flow_accum, output=output, threshold=threshold, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ExtractValleys(object):
    def __init__(self):
        self.label = "Extract Valleys"
        self.description = "Identifies potential valley bottom grid cells based on local topolography alone. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#ExtractValleys' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/extract_valleys.rs' target='_blank'>GitHub</a>."
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
        variant.filter.list = ['LQ', 'JandR', 'PandD']

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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.extract_valleys(dem, output=output, variant=variant, line_thin=line_thin, filter=filter)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FarthestChannelHead(object):
    def __init__(self):
        self.label = "Farthest Channel Head"
        self.description = "Calculates the distance to the furthest upstream channel head for each stream cell. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#FarthestChannelHead' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/farthest_channel_head.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.farthest_channel_head(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class FindMainStem(object):
    def __init__(self):
        self.label = "Find Main Stem"
        self.description = "Finds the main stem, based on stream lengths, of each stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#FindMainStem' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/find_main_stem.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.find_main_stem(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HackStreamOrder(object):
    def __init__(self):
        self.label = "Hack Stream Order"
        self.description = "Assigns the Hack stream order to each tributary in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#HackStreamOrder' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/hack_order.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.hack_stream_order(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class HortonStreamOrder(object):
    def __init__(self):
        self.label = "Horton Stream Order"
        self.description = "Assigns the Horton stream order to each tributary in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#HortonStreamOrder' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/horton_order.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.horton_stream_order(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LengthOfUpstreamChannels(object):
    def __init__(self):
        self.label = "Length Of Upstream Channels"
        self.description = "Calculates the total length of channels upstream. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#LengthOfUpstreamChannels' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/total_length_channels.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.length_of_upstream_channels(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LongProfile(object):
    def __init__(self):
        self.label = "Long Profile"
        self.description = "Plots the stream longitudinal profiles for one or more rivers. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#LongProfile' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/long_profile.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.long_profile(d8_pntr, streams=streams, dem=dem, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class LongProfileFromPoints(object):
    def __init__(self):
        self.label = "Long Profile From Points"
        self.description = "Plots the longitudinal profiles from flow-paths initiating from a set of vector points. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#LongProfileFromPoints' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/long_profile_from_points.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.long_profile_from_points(d8_pntr, points=points, dem=dem, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterStreamsToVector(object):
    def __init__(self):
        self.label = "Raster Streams To Vector"
        self.description = "Converts a raster stream file into a vector file. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#RasterStreamsToVector' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/raster_streams_to_vector.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.raster_streams_to_vector(streams, d8_pntr=d8_pntr, output=output, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RasterizeStreams(object):
    def __init__(self):
        self.label = "Rasterize Streams"
        self.description = "Rasterizes vector streams based on Lindsay (2016) method. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#RasterizeStreams' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/rasterize_streams.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.rasterize_streams(streams, base=base, output=output, nodata=nodata, feature_id=feature_id)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class RemoveShortStreams(object):
    def __init__(self):
        self.label = "Remove Short Streams"
        self.description = "Removes short first-order streams from a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#RemoveShortStreams' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/remove_short_streams.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.remove_short_streams(d8_pntr, streams=streams, output=output, min_length=min_length, esri_pntr=esri_pntr)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class ShreveStreamMagnitude(object):
    def __init__(self):
        self.label = "Shreve Stream Magnitude"
        self.description = "Assigns the Shreve stream magnitude to each link in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#ShreveStreamMagnitude' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/shreve_magnitude.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.shreve_stream_magnitude(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StrahlerStreamOrder(object):
    def __init__(self):
        self.label = "Strahler Stream Order"
        self.description = "Assigns the Strahler stream order to each link in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#StrahlerStreamOrder' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/strahler_order.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.strahler_stream_order(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StreamLinkClass(object):
    def __init__(self):
        self.label = "Stream Link Class"
        self.description = "Identifies the exterior/interior links and nodes in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#StreamLinkClass' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/stream_link_class.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.stream_link_class(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StreamLinkIdentifier(object):
    def __init__(self):
        self.label = "Stream Link Identifier"
        self.description = "Assigns a unique identifier to each link in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#StreamLinkIdentifier' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/stream_link_id.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.stream_link_identifier(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StreamLinkLength(object):
    def __init__(self):
        self.label = "Stream Link Length"
        self.description = "Estimates the length of each link (or tributary) in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#StreamLinkLength' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/stream_link_length.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.stream_link_length(d8_pntr, linkid=linkid, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StreamLinkSlope(object):
    def __init__(self):
        self.label = "Stream Link Slope"
        self.description = "Estimates the average slope of each link (or tributary) in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#StreamLinkSlope' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/stream_link_slope.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.stream_link_slope(d8_pntr, linkid=linkid, dem=dem, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class StreamSlopeContinuous(object):
    def __init__(self):
        self.label = "Stream Slope Continuous"
        self.description = "Estimates the slope of each grid cell in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#StreamSlopeContinuous' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/stream_slope_continuous.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.stream_slope_continuous(d8_pntr, streams=streams, dem=dem, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TopologicalStreamOrder(object):
    def __init__(self):
        self.label = "Topological Stream Order"
        self.description = "Assigns each link in a stream network its topological order. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#TopologicalStreamOrder' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/topological_stream_order.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.topological_stream_order(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return


class TributaryIdentifier(object):
    def __init__(self):
        self.label = "Tributary Identifier"
        self.description = "Assigns a unique identifier to each tributary in a stream network. View detailed help documentation on <a href='https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html#TributaryIdentifier' target='_blank'>WhiteboxTools User Manual</a> and source code on <a href='https://github.com/jblindsay/whitebox-tools//tree/master/src/tools/stream_network_analysis/tributary_id.rs' target='_blank'>GitHub</a>."
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
        old_stdout = sys.stdout
        result = StringIO()
        sys.stdout = result
        wbt.tributary_identifier(d8_pntr, streams=streams, output=output, esri_pntr=esri_pntr, zero_background=zero_background)
        sys.stdout = old_stdout
        result_string = result.getvalue()
        messages.addMessage(result_string)
        return

#
# Private functions
#
def ConvertToFullPath(value):
    '''
    This function takes a single input which could be a layer object or a full path to a dataset source. If it is a layer object it returns the full path to the dataset
    otherwise it simply returns the data source.

    This is done as the user could have navigated to a data source via the browse button or selected the layer as it was already loaded in the map document.

    WhiteboxTools requires a full path to the dataset.
    '''
    desc = arcpy.Describe(value)
    return desc.catalogPath
/*
 * =====================================================================================
 *
 *       Filename:  FeatureSurface.cxx
 *
 *    Description:  Extract Topologically distinct surfaces form data set
 *
 *        Version:  1.0
 *        Created:  12/14/2013 03:46:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

// =======
// Headers
// =======

#include "FeatureSurface.h"

// Pipeline
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkDemandDrivenPipeline.h>

// Data
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

// Arrays
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkDataArray.h>

// General
#include <vtkSmartPointer.h>
#include <vtkMath.h>

// Algorithms
#include <vtkDataSetSurfaceFilter.h>
#include <vtkFeatureEdges.h>

// Test
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>

// ======
// Macros
// ======

vtkStandardNewMacro(FeatureSurface);
vtkCxxRevisionMacro(FeatureSurface,"$Revision 1.0$");
#define EPSILON 1e-4

// ===========
// Constructor
// ===========

FeatureSurface::FeatureSurface()
{
    // Default Member data
    this->DebugStatus = false;
    this->FeatureAngle = 30;
    this->SurfaceColorMode = SurfaceColor::POLY_COLOR;
    this->EdgeColorMode = EdgeColor::OPPOSITE_SURFACE_COLOR;
}

// ==========
// Destructor
// ==========

FeatureSurface::~FeatureSurface()
{
}

// ====================
// Accessors / Mutators
// ====================

// Set Surface Color
void FeatureSurface::SetSurfaceColorMode(int SurfaceColorInt)
{
    switch(SurfaceColorInt)
    {
        // Body Color
        case static_cast<int>(SurfaceColor::BODY_COLOR):
        {
            this->SurfaceColorMode = SurfaceColor::BODY_COLOR;
            break;
        }

        // Mono Color
        case static_cast<int>(SurfaceColor::MONO_COLOR):
        {
            this->SurfaceColorMode = SurfaceColor::MONO_COLOR;
            break;
        }

        // Poly color
        case static_cast<int>(SurfaceColor::POLY_COLOR):
        {
            this->SurfaceColorMode = SurfaceColor::POLY_COLOR;
            break;
        }

        // Undefined mode
        default:
        {
            vtkErrorMacro("Surface color is undefined.");
        }
    }
}

// Set Edge Color Mode
void FeatureSurface::SetEdgeColorMode(int EdgeColorInt)
{
    switch(EdgeColorInt)
    {
        // Body color
        case static_cast<int>(EdgeColor::BODY_COLOR):
        {
            this->EdgeColorMode = EdgeColor::BODY_COLOR;
            break;
        }

        // Mono color
        case static_cast<int>(EdgeColor::MONO_COLOR):
        {
            this->EdgeColorMode = EdgeColor::MONO_COLOR;
            break;
        }

        // Surface color
        case static_cast<int>(EdgeColor::SURFACE_COLOR):
        {
            this->EdgeColorMode = EdgeColor::SURFACE_COLOR;
            break;
        }

        // Opposite color
        case static_cast<int>(EdgeColor::OPPOSITE_SURFACE_COLOR):
        {
            this->EdgeColorMode = EdgeColor::OPPOSITE_SURFACE_COLOR;
            break;
        }

        // Undefined color
        default:
        {
            vtkErrorMacro("Edge color mode is undefined.");
        }
    }
}

// ==========
// Print Self
// ==========

void FeatureSurface::PrintSelf(ostream &os,vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

// ============================
// Fill Output Port Information
// ============================

int FeatureSurface::FillOutputPortInformation(int port,vtkInformation *info)
{
    if(port == 0)
    {
        info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkPolyData");
        return 1;
    }

    return 0;
}

// ============
// Request Data
// ============

int FeatureSurface::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
    // Input
    vtkInformation *InputInfo = inputVector[0]->GetInformationObject(0);
    vtkDataSet *InputDataSet = vtkDataSet::SafeDownCast(InputInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Output
    vtkInformation *OutputInfo = outputVector->GetInformationObject(0);

    // Extract Surface
    vtkSmartPointer<vtkDataSetSurfaceFilter> DataSetSurfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
        DataSetSurfaceFilter->SetInput(InputDataSet);
    #else
        DataSetSurfaceFilter->SetInputData(InputDataSet);
    #endif
    DataSetSurfaceFilter->PassThroughPointIdsOn();
    DataSetSurfaceFilter->Update();

    // Input DataSet Surface
    vtkPolyData *InputDataSetSurface = DataSetSurfaceFilter->GetOutput();

    // Feature Edges
    vtkSmartPointer<vtkFeatureEdges> FeatureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
    FeatureEdges->SetInputConnection(DataSetSurfaceFilter->GetOutputPort());
    FeatureEdges->BoundaryEdgesOff();
    FeatureEdges->ManifoldEdgesOff();
    FeatureEdges->NonManifoldEdgesOff();
    FeatureEdges->FeatureEdgesOn();
    FeatureEdges->SetFeatureAngle(this->FeatureAngle);
    FeatureEdges->ColoringOn();
    FeatureEdges->Update();

    // Edges
    vtkPolyData *Edges = FeatureEdges->GetOutput();

    // Check Edges
    if(Edges->GetNumberOfPoints() < 3)
    {
        // Number of edges are not enough
        vtkErrorMacro("Dataset does not appear to have edges.");
        return 0;
    }

    // Separated Edges
    vtkSmartPointer<vtkPolyData> SeparatedEdges = vtkSmartPointer<vtkPolyData>::New();
    unsigned int NumberOfSeparatedEdges = this->FindSeparatedEdges(Edges,SeparatedEdges,InputDataSet);

    // Check if Separated Edges found
    if(!NumberOfSeparatedEdges)
    {
        vtkErrorMacro("No separated edges found.");
        return 0;
    }

    // Add Original and Surface Ids to Separated Edges
    this->AddIdConversionArraysToSeparatedEdges(
            Edges,
            SeparatedEdges,
            InputDataSetSurface);

    // Test //
    // this->ConvertEdgeToSurfaceIdsTest(SeparatedEdges,InputDataSetSurface);
    // this->OriginalIdTest(InputDataSet,SeparatedEdges);
    // this->OriginalIdTest(InputDataSet,InputDataSetSurface);

    // Find Edge Surfaces
    vtkSmartPointer<vtkCellArray> FeatureSurfaceVertices = vtkSmartPointer<vtkCellArray>::New();
    this->FindSeparatedEdgeSurfaces(
            InputDataSet,
            InputDataSetSurface,
            Edges,
            SeparatedEdges,
            FeatureSurfaceVertices);  // Output

    // Colorize Suface
    vtkSmartPointer<vtkIntArray> SurfaceColors = vtkSmartPointer<vtkIntArray>::New();
    this->ColorizeSurface(
            InputDataSetSurface,
            SeparatedEdges,
            FeatureSurfaceVertices,
            SurfaceColors);  // Output

    // Create Output
    vtkSmartPointer<vtkPolyData> OutputDataSetSurface = vtkSmartPointer<vtkPolyData>::New();
    this->GenerateOutputDataSetSurface(
            InputDataSetSurface,
            SeparatedEdges,
            FeatureSurfaceVertices,
            SurfaceColors,
            OutputDataSetSurface);  // Output

    // Set Output
    OutputInfo->Set(vtkDataObject::DATA_OBJECT(),OutputDataSetSurface);

    // Write to file (for Debuggind) //

    // 1- input DataSet
    // vtkSmartPointer<vtkPolyDataWriter> SurfaceWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    // SurfaceWriter->SetInputData(InputDataSetSurface);
    // SurfaceWriter->SetFileName("/home/sia/Desktop/input_DataSetSurface.vtk");
    // SurfaceWriter->Update();
    
    // 2- output DataSet
    // vtkSmartPointer<vtkUnstructuredGridWriter> DataSetWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    // DataSetWriter->SetInputData(InputDataSet);
    // DataSetWriter->SetFileName("/home/sia/Desktop/output_DataSet.vtk");
    // DataSetWriter->Update();

    // 3- output Edge
    // vtkSmartPointer<vtkPolyDataWriter> EdgeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    // EdgeWriter->SetInputData(Edges);
    // EdgeWriter->SetFileName("/home/sia/Desktop/output_Edge.vtk");
    // EdgeWriter->Update();

    // 4- output Separated Edges
    // vtkSmartPointer<vtkPolyDataWriter> SeparatedEdgeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    // SeparatedEdgeWriter->SetInputData(SeparatedEdges);
    // SeparatedEdgeWriter->SetFileName("/home/sia/Desktop/output_SeparatedEdge.vtk");
    // SeparatedEdgeWriter->Update();

    return 1;
}

// ====================
// Find Separated Edges
// ====================

unsigned int FeatureSurface::FindSeparatedEdges(
        vtkPolyData *Edges,
        vtkPolyData *SeparatedEdges,
        vtkDataSet *InputDataSet)
{
    // Initialize number of separated edges
    unsigned int NumberOfSeparatedEdges = 0;

    // Number of Edges points
    unsigned int NumberOfEdgePoints = Edges->GetNumberOfPoints();

    if(NumberOfEdgePoints < 3)
    {
        // Number of edges points are not enough for closed polygon.
        return NumberOfSeparatedEdges;
    }

    // Declare Processed Points Array
    bool *PointsProcessed = new bool[NumberOfEdgePoints];
    
    // Initialize Processed Points array with false values
    for(unsigned int PointIterator = 0; PointIterator < NumberOfEdgePoints; PointIterator++)
    {
        PointsProcessed[PointIterator] = false;
    }

    // Edge Polygon and Array of Polygons
    vtkSmartPointer<vtkCellArray> Polygons = vtkSmartPointer<vtkCellArray>::New();

    // Search for Polygon Rings
    for(unsigned int InquiryPointId = 0; InquiryPointId < NumberOfEdgePoints; InquiryPointId++)
    {
        // Leave out processed points
        if(PointsProcessed[InquiryPointId] == true)
        {
            continue;
        }

        // Declare a polygon for new edge
        vtkSmartPointer<vtkPolygon> Polygon = vtkSmartPointer<vtkPolygon>::New();

        // Search for polygon ring
        bool ClosedPolygonFound = this->FindClosedPolygon(
                InquiryPointId,
                Edges,
                PointsProcessed,
                Polygon);

        if(ClosedPolygonFound)
        {
            // Add Polygon to array of Polygons
            Polygons->InsertNextCell(Polygon);
        }
    }

    // Free memory
    if(PointsProcessed != NULL)
    {
        delete [] PointsProcessed;
        PointsProcessed = NULL;
    }

    // Check if Polygons has at least a polygon
    if(Polygons->GetNumberOfCells() < 1)
    {
        NumberOfSeparatedEdges = 0;

        if(this->DebugStatus)
        {
            std::cout << "No separated edges found. " << std::endl;
        }
    }
    else
    {
        // PolyData of Polygons
        SeparatedEdges->SetPoints(Edges->GetPoints());
        SeparatedEdges->SetLines(Edges->GetLines());
        SeparatedEdges->SetPolys(Polygons);

        // Number of Separated Edges
        NumberOfSeparatedEdges = SeparatedEdges->GetPolys()->GetNumberOfCells();

        if(this->DebugStatus)
        {
            std::cout << "Total number of closed polygons: " << NumberOfSeparatedEdges << std::endl;
        }
    }

    return NumberOfSeparatedEdges;
}

// ===================
// Find Closed Polygon
// ===================

bool FeatureSurface::FindClosedPolygon(
        vtkIdType InquiryPointId,
        vtkPolyData *Edges,
        bool *PointsProcessed,
        vtkPolygon *Polygon)
{
    // Status of output
    bool ClosedPolygonFound = true;

    // Find connectivity of Inquiry point
    vtkSmartPointer<vtkIdList> InquiryPointConnectivities = vtkSmartPointer<vtkIdList>::New();
    this->GetPointConnectivities(InquiryPointId,Edges,InquiryPointConnectivities);

    // Check if inquiry point has at least a connectivity
    if(InquiryPointConnectivities->GetNumberOfIds() < 2)
    {
        // Inquiry point can not be on a closed polygon
        ClosedPolygonFound = false;
        return ClosedPolygonFound;
    }

    // Choose an adjacent point to the inquiry point
    vtkIdType AdjacentToInquiryPointId = InquiryPointConnectivities->GetId(0);

    // Add Inquiry point and it's adjacent to the polygon
    Polygon->GetPointIds()->InsertId(0,InquiryPointId);
    Polygon->GetPointIds()->InsertId(1,AdjacentToInquiryPointId);

    PointsProcessed[InquiryPointId] = true;
    PointsProcessed[AdjacentToInquiryPointId] = true;

    // Declare and Shift Current, Previous and Next point for marching on edge
    vtkIdType PreviousAdjacentPointId = InquiryPointId;
    vtkIdType CurrentPointId = AdjacentToInquiryPointId;
    vtkIdType NextAdjacentPointId;

    // Search for Polygon ring for current point
    bool BoundaryRepeated = false;

    // Marching on Edge for a closed polygon
    while(!BoundaryRepeated)
    {
        NextAdjacentPointId = this->FindAdjacentEdgePoint(
                CurrentPointId,
                PreviousAdjacentPointId,
                Edges);

        // Check if next adjacent point is valid
        if(NextAdjacentPointId < 0)
        {
            // PointsProcessed[NextAdjacentPointId] = true;
            ClosedPolygonFound = false;
            break;
        }

        // Check if Polygon is closed
        if(NextAdjacentPointId == InquiryPointId)
        {
            // Polygin closed.
            BoundaryRepeated = true;
            break;
        }

        // Check if next adjacent point is processed before
        if(PointsProcessed[NextAdjacentPointId] == true)
        {
            // Next adjacent point was processed before but the polygon is not closed yet.
            ClosedPolygonFound = false;
            return ClosedPolygonFound;
        }
        else
        {
            PointsProcessed[NextAdjacentPointId] = true;
        }

        // Add the next adjacent point to the polygon
        Polygon->GetPointIds()->InsertNextId(NextAdjacentPointId);

        // Update Marching points
        PreviousAdjacentPointId = CurrentPointId;
        CurrentPointId = NextAdjacentPointId;
    }

    // Check if polygon has enough points
    if(Polygon->GetPointIds()->GetNumberOfIds() < 3)
    {
        // Polygon does not have enough points.
        ClosedPolygonFound = false;
    }

    // A closed polygon found
    if(ClosedPolygonFound)
    {
        if(this->DebugStatus)
        {
            std::cout << "Closed polygon found. Number of points: " << 
                Polygon->GetPointIds()->GetNumberOfIds() << std::endl;
        }
    }

    return ClosedPolygonFound;
}

// ========================
// Find Adjacent Edge Point
// ========================

int FeatureSurface::FindAdjacentEdgePoint(
        vtkIdType CurrentPointId,
        vtkIdType PreviousAdjacentPointId,
        vtkPolyData *Edges)
{
    // Current Point Connectivities
    vtkSmartPointer<vtkIdList> CurrentPointConnectivities = vtkSmartPointer<vtkIdList>::New();
    this->GetPointConnectivities(CurrentPointId,Edges,CurrentPointConnectivities);

    // Connectivities should only contain two points on Edges
    unsigned int NumberOfConnectivities = CurrentPointConnectivities->GetNumberOfIds();
    if(NumberOfConnectivities != 2)
    {
        if(this->DebugStatus)
        {
            std::cout << "Fatal Point: " << CurrentPointId << 
                ", Number of connectivities: " << NumberOfConnectivities << std::endl;
        }
        return -1;
    }

    int NextAdjacentPointId;

    // Choose new adjacent point but not the previous adjacent point
    if(CurrentPointConnectivities->GetId(0) == PreviousAdjacentPointId)
    {
        NextAdjacentPointId = CurrentPointConnectivities->GetId(1);
    }
    else
    {
        NextAdjacentPointId = CurrentPointConnectivities->GetId(0);
    }

    return NextAdjacentPointId;
}

// ========================
// Get Point Connectivities
// ========================

void FeatureSurface::GetPointConnectivities(
        vtkIdType InquiryPointId,
        vtkDataSet *DataSet,
        vtkIdList *Connectivities)
{
    // Neighbor Cells of a point
    vtkSmartPointer<vtkIdList> NeighborCellsIdList = vtkSmartPointer<vtkIdList>::New();
    DataSet->GetPointCells(InquiryPointId,NeighborCellsIdList);
    vtkIdType NumberOfNeighborCells = NeighborCellsIdList->GetNumberOfIds();

    // Initialize number of connectivities
    vtkIdType NumberOfConnectivities = 0;

    // Loop over neighbor cells
    for(unsigned int NeighborCellIterator = 0;
        NeighborCellIterator < static_cast<unsigned int>(NumberOfNeighborCells);
        NeighborCellIterator++)
    {
        // Select a cell from neighbor cells
        vtkIdType SelectedCellId = NeighborCellsIdList->GetId(NeighborCellIterator);

        // Get points of selected cell
        vtkSmartPointer<vtkIdList> PointsOfSelectedCellList = vtkSmartPointer<vtkIdList>::New();
        DataSet->GetCellPoints(SelectedCellId,PointsOfSelectedCellList);
        vtkIdType NumberOfCellPoints = PointsOfSelectedCellList->GetNumberOfIds();

        // Loop over points in selected neighbor cell
        for(unsigned int PointIterator = 0;
            PointIterator < static_cast<unsigned int>(NumberOfCellPoints);
            PointIterator++)
        {
            // Select a point from points of selected neighbor cell
            vtkIdType CheckPointId = PointsOfSelectedCellList->GetId(PointIterator);
            
            // Check if the selected point is not the original inquiry point
            if(CheckPointId != InquiryPointId)
            {
                bool RepeatedBefore = false;

                // Check if the selected point is one of the previous connectivity points or not
                if(NumberOfConnectivities != 0)
                {
                    // Loop over previous connected points
                    for(vtkIdType ConnectivityIterator = 0;
                        ConnectivityIterator < NumberOfConnectivities;
                        ConnectivityIterator++)
                    {
                        // Get a previous connectivity point
                        vtkIdType PreviousConnectedPoint = Connectivities->GetId(ConnectivityIterator);

                        // Check if selected point is not the previous connected point
                        if(CheckPointId == PreviousConnectedPoint)
                        {
                            RepeatedBefore = true;
                        }
                    }
                }
                
                // Add new checked points if it is not connected before
                if(RepeatedBefore == false)
                {
                    Connectivities->InsertNextId(CheckPointId);
                    NumberOfConnectivities++;
                }
            }
        }
    }
}

// ============================
// Find Separated Edge Surfaces
// ============================

void FeatureSurface::FindSeparatedEdgeSurfaces(
        vtkDataSet *InputDataSet,
        vtkPolyData *InputDataSetSurface,
        vtkPolyData *Edges,
        vtkPolyData *SeparatedEdges,
        vtkCellArray *FeatureSurfaceVertices)
{
    // Polygon Arrays
    vtkSmartPointer<vtkCellArray> ClosedEdgesArray = SeparatedEdges->GetPolys();
    
    // Number of Closed Edges
    unsigned int NumberOfClosedEdges = ClosedEdgesArray->GetNumberOfCells();

    // Initialize output
    FeatureSurfaceVertices->InitTraversal();

    // Loop over Separated Edges
    ClosedEdgesArray->InitTraversal();
    for(unsigned int ClosedEdgeIterator = 0; ClosedEdgeIterator < NumberOfClosedEdges; ClosedEdgeIterator++)
    {
        // Get Closed Edge as IdList
        vtkSmartPointer<vtkIdList> ClosedEdge = vtkSmartPointer<vtkIdList>::New();
        ClosedEdgesArray->GetNextCell(ClosedEdge);

        // Find Geometric Center of Closed Edge
        double GeometricCenterOfClosedEdge[3];
        this->FindGeometricCenterOfClosedEdge(
                InputDataSet,
                SeparatedEdges,
                ClosedEdge,
                GeometricCenterOfClosedEdge);  // Output

        // Find a candidate point on Edge surface
        vtkIdType ACandidatePointOnEdgeSurfaceSurfaceId;
        this->FindACandidatePointOnEdgeSurface(
                InputDataSet,
                InputDataSetSurface,
                Edges,
                SeparatedEdges,
                ClosedEdge,
                GeometricCenterOfClosedEdge,
                ACandidatePointOnEdgeSurfaceSurfaceId);   // Output

        // Find portion of surface that encircled with the closed edge
        vtkSmartPointer<vtkIdList> ClosedEdgeSurface = vtkSmartPointer<vtkIdList>::New();

        this->FindClosedEdgeSurface(
                InputDataSetSurface,
                ClosedEdge,
                SeparatedEdges,
                ACandidatePointOnEdgeSurfaceSurfaceId,
                ClosedEdgeSurface);   // Output

        // Debug //
        if(this->DebugStatus)
        {
            std::cout << "Closed Edge: " << ClosedEdgeIterator+1 << ",\t Surface points found: " <<
                ClosedEdgeSurface->GetNumberOfIds() << std::endl;
        }

        // Add ClosedEdgeSurface to a vertex
        FeatureSurfaceVertices->InsertNextCell(ClosedEdgeSurface);
    }
}

// ====================================
// Find Geometric Center of Closed Edge
// ====================================

void FeatureSurface::FindGeometricCenterOfClosedEdge(
        vtkDataSet *InputDataSet,
        vtkPolyData *SeparatedEdges,
        vtkIdList *ClosedEdge,
        double GeometricCenter[3])
{
    // Original Edge Point Ids
    vtkSmartPointer<vtkIdTypeArray> OriginalEdgePointIdsArray = 
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalPointIds"));
    if(!OriginalEdgePointIdsArray)
    {
        vtkErrorMacro("Original Edge point ids array is NULL.");
    }

    // Initialize Geometric Center
    for(unsigned int Dimension = 0; Dimension < 3; Dimension ++)
    {
        GeometricCenter[Dimension] = 0;
    }

    // Number of Edge Points
    unsigned int NumberOfEdgePoints = ClosedEdge->GetNumberOfIds();

    // Averaging point coordinates
    for(unsigned int PointIterator = 0; PointIterator < NumberOfEdgePoints; PointIterator++)
    {
        // Point Id on Edge dataset
        vtkIdType PointId = ClosedEdge->GetId(PointIterator);

        // Point Id on original dataset
        vtkIdType OriginalPointId = OriginalEdgePointIdsArray->GetTuple1(PointId);

        // Point Coordinates
        double PointCoordinates[3];
        InputDataSet->GetPoint(OriginalPointId,PointCoordinates);

        // Summing coordinates for averaging
        for(unsigned int Dimension = 0; Dimension < 3; Dimension++)
        {
            GeometricCenter[Dimension] += PointCoordinates[Dimension];
        }
    }

    // Averaging coordinates
    for(unsigned int Dimension = 0; Dimension < 3; Dimension++)
    {
        GeometricCenter[Dimension] /= double(NumberOfEdgePoints);
    }
}

// ======================================
// Find A Candidate Point On Edge Surface
// ======================================

void FeatureSurface::FindACandidatePointOnEdgeSurface(
        vtkDataSet *InputDataSet,
        vtkPolyData *InputDataSetSurface,
        vtkPolyData *Edges,
        vtkPolyData *SeparatedEdges,
        vtkIdList *ClosedEdge,
        double GeometricCenterOfClosedEdge[3],
        vtkIdType &ACandidatePointOnEdgeSurfaceSurfaceId)
{
    // Check polygon
    if(ClosedEdge == NULL)
    {
        vtkErrorMacro("Closed Edge is NULL");
    }
    else if(ClosedEdge->GetNumberOfIds() < 3)
    {
        vtkErrorMacro("ClosedEdge does not have enough points.");
    }

    // Original Edge Point Ids
    vtkSmartPointer<vtkIdTypeArray> OriginalEdgePointIdsArray = 
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalPointIds"));
    if(!OriginalEdgePointIdsArray)
    {
        vtkErrorMacro("Original edge point ids array is NULL.");
    }

    // Pick a point on edge
    vtkIdType EdgePointId = ClosedEdge->GetId(0);
    vtkIdType EdgePointOriginalId = OriginalEdgePointIdsArray->GetTuple1(EdgePointId);
    double EdgePointCoordinates[3];
    InputDataSet->GetPoint(EdgePointOriginalId,EdgePointCoordinates);

    // Select an adjacent point of inquiry point on the edge, used for finding Ray angle to Polygon surface
    vtkSmartPointer<vtkIdList> EdgePointConnectivitiesOnEdge = vtkSmartPointer<vtkIdList>::New();
    this->GetPointConnectivities(EdgePointId,Edges,EdgePointConnectivitiesOnEdge);
    vtkIdType AdjacentOfEdgePointId = EdgePointConnectivitiesOnEdge->GetId(0);
    vtkIdType AdjacentOfEdgePointOriginalId = OriginalEdgePointIdsArray->GetTuple1(AdjacentOfEdgePointId);
    double AdjacentOfEdgePointCoordinates[3];
    InputDataSet->GetPoint(AdjacentOfEdgePointOriginalId,AdjacentOfEdgePointCoordinates);

    // Neighbor points on edge point on the surface, but not on the edge (excluding two adjacent points)
    vtkSmartPointer<vtkIdList> EdgePointConnectivitiesOnSurfaceExcludingEdge =
        vtkSmartPointer<vtkIdList>::New();
    this->GetEdgePointConnectivitiesOnSurfaceExcludingEdge(
            EdgePointId,
            EdgePointConnectivitiesOnEdge,
            InputDataSetSurface,
            SeparatedEdges,
            EdgePointConnectivitiesOnSurfaceExcludingEdge);

    // convert Surface Ids to Original Ids of Connectivities
    vtkSmartPointer<vtkIdList> EdgePointOriginalConnectivitiesOnSurfaceExcludingEdge =
        vtkSmartPointer<vtkIdList>::New();
    this->ConvertSurfaceIdListToOriginalIdList(
            EdgePointConnectivitiesOnSurfaceExcludingEdge,
            InputDataSetSurface,
            EdgePointOriginalConnectivitiesOnSurfaceExcludingEdge);

    // Number of Neighbor points
    unsigned int NumberOfNeighborPoints = EdgePointOriginalConnectivitiesOnSurfaceExcludingEdge->GetNumberOfIds();

    if(NumberOfNeighborPoints < 2)
    {
        vtkErrorMacro("Not enough surface neighbor points (excluding edges) for the edge point.");
    }

    // Ray Angle of neighbor point toward edge point w.r.t surface
    std::vector<double> RayAngleProjectionVector(NumberOfNeighborPoints);

    // Loop over Neighbor points
    for(unsigned int NeighborPointIterator = 0;
        NeighborPointIterator < NumberOfNeighborPoints;
        NeighborPointIterator++)
    {
        // Get Neighbor point Id
        vtkIdType NeighborOfEdgePointOriginalId = 
            EdgePointOriginalConnectivitiesOnSurfaceExcludingEdge->GetId(NeighborPointIterator);

        // Get Neighbor of edge point coordinates
        double NeighborOfEdgePointCoordinates[3];
        InputDataSet->GetPoint(NeighborOfEdgePointOriginalId,NeighborOfEdgePointCoordinates);

        // Find Ray Angle of edge point toward neighbor point w.r.t surface
        RayAngleProjectionVector[NeighborPointIterator] = GetRayAngleProjectionOnSurface(
                EdgePointCoordinates,
                AdjacentOfEdgePointCoordinates,
                GeometricCenterOfClosedEdge,
                NeighborOfEdgePointCoordinates);
    }

    // Candidate point is a point with maximum of projection of  Ray Angle w.r.t surface
    size_t MaxRayAngleProjectionIndex = GetMaxIndexOfVector(RayAngleProjectionVector);
    ACandidatePointOnEdgeSurfaceSurfaceId = 
        EdgePointConnectivitiesOnSurfaceExcludingEdge->GetId(MaxRayAngleProjectionIndex);
}

// =======================================================
// Get Edge Point Connectivities On Surface Excluding Edge
// =======================================================

bool FeatureSurface::GetEdgePointConnectivitiesOnSurfaceExcludingEdge(
        vtkIdType EdgePointId,
        vtkIdList *EdgePointConnectivitiesOnEdge,
        vtkPolyData *InputDataSetSurface,
        vtkPolyData *SeparatedEdges,
        vtkIdList *EdgePointConnectivitiesOnSurfaceExcludingEdge)
{
    // Two adjacent points of the edge point
    vtkIdType LeftAdjacentPointId = EdgePointConnectivitiesOnEdge->GetId(0);
    vtkIdType RightAdjacentPointId = EdgePointConnectivitiesOnEdge->GetId(1);

    // Converting Ids to their surface Id
    vtkSmartPointer<vtkIdTypeArray> EdgePointsSurfaceIds =
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalSurfaceIds"));

    if(!EdgePointsSurfaceIds)
    {
        vtkErrorMacro("EdgePointsSurfaceIds is NULL.");
    }

    vtkIdType EdgePointSurfaceId = EdgePointsSurfaceIds->GetTuple1(EdgePointId);
    vtkIdType LeftAdjacentPointSurfaceId = EdgePointsSurfaceIds->GetTuple1(LeftAdjacentPointId);
    vtkIdType RightAdjacentPointSurfaceId = EdgePointsSurfaceIds->GetTuple1(RightAdjacentPointId);

    // Get Connectivities of edge point on surface with surface Ids
    vtkSmartPointer<vtkIdList> EdgePointConnectivitiesOnSurface =
        vtkSmartPointer<vtkIdList>::New();
    this->GetPointConnectivities(EdgePointSurfaceId,InputDataSetSurface,EdgePointConnectivitiesOnSurface);

    // Exclude Adjacent points from connectivities
    for(unsigned int NeighborPointIterator = 0; 
        NeighborPointIterator < static_cast<unsigned int>(EdgePointConnectivitiesOnSurface->GetNumberOfIds());
        NeighborPointIterator++)
    {
        // Get Surface Id of current neighbor point
        vtkIdType NeighborPointSurfaceId = EdgePointConnectivitiesOnSurface->GetId(NeighborPointIterator);

        // Compare with Adjacent points
        if(NeighborPointSurfaceId != LeftAdjacentPointSurfaceId &&
           NeighborPointSurfaceId != RightAdjacentPointSurfaceId)
        {
            EdgePointConnectivitiesOnSurfaceExcludingEdge->InsertNextId(NeighborPointSurfaceId);
        }
    }

    // Id List need to have more than one Id
    if(EdgePointConnectivitiesOnSurfaceExcludingEdge->GetNumberOfIds() > 1)
    {
        // Enough points for looking for candidate point
        return true;
    }
    else
    {
        // Not enough points for looking for candidate point
        return false;
    }
}

// =========================================
// Convert Surface IdList to Original IdList
// =========================================

void FeatureSurface::ConvertSurfaceIdListToOriginalIdList(
        vtkIdList *SurfaceIdList,
        vtkPolyData *InputDataSetSurface,
        vtkIdList *OriginalIdList)
{
    // Initialize Output
    OriginalIdList->Initialize();

    // Surface Points Original Ids
    vtkSmartPointer<vtkIdTypeArray> SurfacePointsOriginalIds = 
        vtkIdTypeArray::SafeDownCast(InputDataSetSurface->GetPointData()->GetArray("vtkOriginalPointIds"));

    if(!SurfacePointsOriginalIds)
    {
        vtkErrorMacro("SurfacePointsOriginalIds is NULL.");
    }

    // Iterate over Ids
    for(unsigned int SurfaceIdIterator = 0;
        SurfaceIdIterator < static_cast<unsigned int>(SurfaceIdList->GetNumberOfIds());
        SurfaceIdIterator++)
    {
        // Get a Surface Id
        vtkIdType SurfaceId = SurfaceIdList->GetId(SurfaceIdIterator);

        // Find Original Id
        vtkIdType OriginalId = SurfacePointsOriginalIds->GetTuple1(SurfaceId);

        // Insert into Original IdList
        OriginalIdList->InsertNextId(OriginalId);
    }
}

// ===================================
// Get Ray Angle Projection On Surface
// ===================================

// Description:
// Computes projection of a ray angle onto the Closed edge surface (Polygon face)
// Ray is defined as a vector that starts from a selected point on edge and
// ends to a neighbor point of the edge point. The neighbor point is on the
// surface, but it is not on edge itself.
// The Ray angle is the angle between the ray and the polygon surface (closed edge).

double FeatureSurface::GetRayAngleProjectionOnSurface(
        double EdgePointCoordinates[3],
        double AdjacentOfEdgePointCoordinates[3],
        double GeometricCenterOfClosedEdge[3],
        double NeighborOfEdgePointCoordinates[3])
{
    // Define points as vectors from Edge point
    double AdjacentVector[3];
    double GeometricCenterVector[3];
    double NeighborVector[3];

    // Subtract from edge point to get vectors
    vtkMath::Subtract(AdjacentOfEdgePointCoordinates,EdgePointCoordinates,AdjacentVector);
    vtkMath::Subtract(GeometricCenterOfClosedEdge,EdgePointCoordinates,GeometricCenterVector);
    vtkMath::Subtract(NeighborOfEdgePointCoordinates,EdgePointCoordinates,NeighborVector);

    // Normalize to unit vectors
    vtkMath::Normalize(AdjacentVector);
    vtkMath::Normalize(GeometricCenterVector);
    vtkMath::Normalize(NeighborVector);

    // Normal vector of the surface (normal to adjacent and center vector)
    double SurfaceNormalVector[3];
    vtkMath::Cross(AdjacentVector,GeometricCenterVector,SurfaceNormalVector);

    // Compute component of NeighrVector on SurfaceNormalVector
    double NeighborToSurfaceNormal = vtkMath::Dot(NeighborVector,SurfaceNormalVector);
    vtkMath::MultiplyScalar(SurfaceNormalVector,NeighborToSurfaceNormal);

    // Normal component of NeighborVector w.r.t SurfaceNormalVector
    double NeighborOnSurfaceVector[3];
    vtkMath::Subtract(NeighborVector,SurfaceNormalVector,NeighborOnSurfaceVector);
    vtkMath::Normalize(NeighborOnSurfaceVector);

    // Project Angle between NeighborVector and NeighborOnSurfaceVector, onto the surface
    double RayAngleProjectionOnSurface = vtkMath::Dot(NeighborVector,NeighborOnSurfaceVector);

    return RayAngleProjectionOnSurface;
}

// =======================
// Get Max Index Of Vector
// =======================

// Description:
// Returns the Index of maximum value of a vector.
// If vecotr has multiple max values, it returns the first occurentce of max.

size_t GetMaxIndexOfVector(std::vector<double> const &Vector)
{
    // Initialize outout
    size_t MaxIndex = 0;
    double MaxValue = Vector[MaxIndex];

    // Iterate for smaller values
    for(size_t VectorIterator = 0; VectorIterator < Vector.size(); VectorIterator++)
    {
        if(Vector[VectorIterator] > MaxValue)
        {
            MaxValue = Vector[VectorIterator];
            MaxIndex = VectorIterator;
        }
    }

    return MaxIndex;
}

// =========================
//  Find Closed Edge Surface
//  ========================

void FeatureSurface::FindClosedEdgeSurface(
        vtkPolyData *InputDataSetSurface,
        vtkIdList *ClosedEdge,
        vtkPolyData *SeparatedEdges,
        vtkIdType ACandidatePointOnEdgeSurfaceSurfaceId,
        vtkIdList *ClosedEdgeSurface)   // Output
{
    // Number of DataSet Surface Points
    unsigned int NumberOfDataSetSurfacePoints = InputDataSetSurface->GetNumberOfPoints();

    // Track proceess of points
    bool *PointsProcessed = new bool[NumberOfDataSetSurfacePoints];
    for(unsigned int SurfacePointIterator = 0;
        SurfacePointIterator < NumberOfDataSetSurfacePoints;
        SurfacePointIterator++)
    {
        PointsProcessed[SurfacePointIterator] = false;
    }

    // Edge Point Ids on Surface
    vtkSmartPointer<vtkIdTypeArray> EdgePointsSurfaceIds =
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalSurfaceIds"));

    if(!EdgePointsSurfaceIds)
    {
        vtkErrorMacro("EdgePointSurfaceIds is NULL.");
    }

    // Mark all edge points as already-processed points
    unsigned int NumberOfClosedEdgePoints = ClosedEdge->GetNumberOfIds();
    for(unsigned int ClosedEdgePointIterator = 0;
        ClosedEdgePointIterator < NumberOfClosedEdgePoints;
        ClosedEdgePointIterator++)
    {
        // Closed Edge Point Id
        vtkIdType ClosedEdgePointId = ClosedEdge->GetId(ClosedEdgePointIterator);

        // Closed EdgePoint Id on surface
        vtkIdType ClosedEdgePointIdOnSurface = EdgePointsSurfaceIds->GetTuple1(ClosedEdgePointId);

        // Mark the edge point as processed
        PointsProcessed[ClosedEdgePointIdOnSurface] = true;
    }

    // Declare buffer points, those that to be marched latar
    vtkSmartPointer<vtkIdList> SurfacePointIdsBuffer = vtkSmartPointer<vtkIdList>::New();

    // Add the candidate point into the buffer
    SurfacePointIdsBuffer->InsertNextId(ACandidatePointOnEdgeSurfaceSurfaceId);

    // March points, starting from the candidate point
    while(SurfacePointIdsBuffer->GetNumberOfIds())
    {
        // Pick first point in buffer as marching point
        vtkIdType MarchingPointId = SurfacePointIdsBuffer->GetId(0);

        // Add Marching point to output
        if(PointsProcessed[MarchingPointId] == false)
        {
            // Add marching point to output
            ClosedEdgeSurface->InsertNextId(MarchingPointId);
            PointsProcessed[MarchingPointId] = true;

            // Remove current marching point from buffer
            SurfacePointIdsBuffer->DeleteId(MarchingPointId);
        }

        // Get Connectivities of marching point
        vtkSmartPointer<vtkIdList> MarchingPointConnectivities = vtkSmartPointer<vtkIdList>::New();
        this->GetPointConnectivities(MarchingPointId,InputDataSetSurface,MarchingPointConnectivities);

        // Add Connected points into buffer
        for(unsigned int NeighborPointsIterator = 0;
            NeighborPointsIterator < static_cast<unsigned int>(MarchingPointConnectivities->GetNumberOfIds());
            NeighborPointsIterator++)
        {
            // Get Neighbor point Id
            vtkIdType NeighborPointId = MarchingPointConnectivities->GetId(NeighborPointsIterator);

            // Avoid adding processed points
            if(PointsProcessed[NeighborPointId] == false)
            {
                // Add to buffer
                SurfacePointIdsBuffer->InsertNextId(NeighborPointId);
            }
        }
    }
}

// ============
// Id Converter
// ============

int FeatureSurface::IdConverter(
        vtkDataSet *SourceDataSet,
        vtkDataSet *DestinationDataSet,
        vtkIdType SourcePointId)
{
    // Get Source Point coordinates
    double *SourcePointCoordinates = SourceDataSet->GetPoint(SourcePointId);

    // Look for a point in destination dataset with the same source point coordinate
    vtkIdType DestinationPointId = DestinationDataSet->FindPoint(SourcePointCoordinates);

    // Check validity of Destination Id
    if(DestinationPointId < 0)
    {
        vtkErrorMacro("Source point in destination dataset can not be found.");
        return -1;
    }

    // Get Destination Point coordinates
    double *DestinationPointCoordinates = DestinationDataSet->GetPoint(DestinationPointId);

    // Compare point coordinates with a small threshold
    double PointDifference[3];
    vtkMath::Subtract(SourcePointCoordinates,DestinationPointCoordinates,PointDifference);
    double Error = vtkMath::Norm(PointDifference);
    if(Error < EPSILON)
    {
        vtkErrorMacro("Destination point is not close to source point.");
        return -1;
    }
    
    return DestinationPointId;
}

// =======================
// Find Original Point Ids
// =======================

void FeatureSurface::FindOriginalPointIds(
        vtkPolyData *ExtractedDataSet,
        vtkDataSet *OriginalDataSet,
        vtkIdTypeArray *OriginalCellIds,
        vtkIdTypeArray *OriginalPointIds)
{
    // Get Line Cells in extracted dataset
    vtkSmartPointer<vtkCellArray> CellsArray = ExtractedDataSet->GetLines();

    // Number of Line Cell in Extracted DataSet
    unsigned int NumberOfCells = CellsArray->GetNumberOfCells();
    if(NumberOfCells < 1)
    {
        vtkErrorMacro("No line cell to process.");
    }

    // Loop over cells in Extracted dataset
    vtkSmartPointer<vtkIdList> CellPointsIdList = vtkSmartPointer<vtkIdList>::New();

    vtkIdType CellId = 0;
    for(CellsArray->InitTraversal(); CellsArray->GetNextCell(CellPointsIdList); CellId++)
    {
        // Pick a cell
        // vtkSmartPointer<vtkIdList> CellPointsIdList = vtkSmartPointer<vtkIdList>::New();
        // CellsArray->GetNextCell(CellPointsIdList);

        // Find an original cell in original dataset
        // vtkIdType OriginalCellId = OriginalCellIds->GetTuple1(CellId);

        // Loop over points of the selected cell
        unsigned int NumberOfPointsInCell = CellPointsIdList->GetNumberOfIds();
        for(unsigned int PointIterator = 0; PointIterator < NumberOfPointsInCell; PointIterator++)
        {
            // TODO
        }
        
    }

}

// ===============================
// Convert Edge Ids To Surface Ids
// ===============================

// Description:
// Coverting Ids on Edge polydata onto Ids on Surface polydata.
// Based on matching original point ids of both edges and surfaces.
// Original point ids belong to Input dataset.

bool FeatureSurface::ConvertEdgeIdsToSurfaceIds(
        vtkPolyData *Edges,
        vtkPolyData *SeparatedEdges,
        vtkPolyData *InputDataSetSurface,
        vtkIdTypeArray *EdgePointsSurfaceIds)
{
    // Status of output
    bool SurfaceIdFound = false;

    // Edge Original Point Ids Array
    vtkSmartPointer<vtkIdTypeArray> EdgePointsOriginalIdsArray = 
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalPointIds"));

    if(!EdgePointsOriginalIdsArray)
    {
        vtkErrorMacro("Edge points original Ids array is NULL.");
    }

    // Surface Original Point Ids Array
    vtkSmartPointer<vtkIdTypeArray> SurfacePointsOriginalIdsArray = 
        vtkIdTypeArray::SafeDownCast(InputDataSetSurface->GetPointData()->GetArray("vtkOriginalPointIds"));

    if(!SurfacePointsOriginalIdsArray)
    {
        vtkErrorMacro("Surface points original Ids array is NULL.");
    }

    // Initialize Output of this function
    EdgePointsSurfaceIds->SetNumberOfTuples(Edges->GetNumberOfPoints());
    EdgePointsSurfaceIds->SetNumberOfComponents(1);
    EdgePointsSurfaceIds->SetName("vtkOriginalSurfaceIds");

    for(unsigned int IdIterator = 0;
        IdIterator < static_cast<unsigned int>(EdgePointsSurfaceIds->GetNumberOfTuples());
        IdIterator++)
    {
        EdgePointsSurfaceIds->SetTuple1(IdIterator,-1);
    }

    // Iterate over closed edges
    vtkSmartPointer<vtkCellArray> ClosedEdgesArray = SeparatedEdges->GetPolys();
    vtkSmartPointer<vtkIdList> ClosedEdge = vtkSmartPointer<vtkIdList>::New();

    for(ClosedEdgesArray->InitTraversal(); ClosedEdgesArray->GetNextCell(ClosedEdge);)
    {
        // Pick a point to start with
        vtkIdType StartPointId = ClosedEdge->GetId(0);
        vtkIdType StartPointOriginalId = EdgePointsOriginalIdsArray->GetTuple1(StartPointId);
        vtkIdType StartPointIdOnSurface;  // To be found in the next for loop

        // Find Id of start point on the surface, using global id search
        for(unsigned int SurfacePointIterator = 0;
            SurfacePointIterator < static_cast<unsigned int>(InputDataSetSurface->GetNumberOfPoints());
            SurfacePointIterator++)
        {
            // Get Original Surface point Id
            vtkIdType SurfacePointOriginalId = SurfacePointsOriginalIdsArray->GetTuple1(SurfacePointIterator);

            // Compare Original Ids of Edge and Surface
            if(StartPointOriginalId == SurfacePointOriginalId)
            {
                StartPointIdOnSurface = SurfacePointIterator;
                EdgePointsSurfaceIds->SetTuple1(StartPointId,StartPointIdOnSurface);
                SurfaceIdFound = true;
                break;
            }
        }

        // Check if Surface Id not found
        if(SurfaceIdFound == false)
        {
            vtkErrorMacro("No surface point found with the same original Id as the edge original Id.");
            return SurfaceIdFound;
        }

        // Pick an adjacent point to the starting point
        vtkSmartPointer<vtkIdList> StartPointConnectivitiesOnEdge = vtkSmartPointer<vtkIdList>::New();
        this->GetPointConnectivities(StartPointId,Edges,StartPointConnectivitiesOnEdge);
        vtkIdType AdjacentStartPointId = StartPointConnectivitiesOnEdge->GetId(0);

        // Marching on the edge
        vtkIdType PreviousEdgePointId = StartPointId;
        vtkIdType PreviousEdgePointIdOnSurface = StartPointIdOnSurface;
        vtkIdType CurrentEdgePointId = AdjacentStartPointId;
        vtkIdType CurrentEdgePointIdOnSurface;  // To be found in the next while loop

        // Iterate over Points
        while(CurrentEdgePointId != StartPointId)
        {
            // Get Surface connectivities of previous point
            vtkSmartPointer<vtkIdList> PreviousPointConnectivitiesOnSurface = vtkSmartPointer<vtkIdList>::New();
            this->GetPointConnectivities(
                    PreviousEdgePointIdOnSurface,
                    InputDataSetSurface,
                    PreviousPointConnectivitiesOnSurface);

            // Iterate over previous point connectivities
            SurfaceIdFound = false;
            for(unsigned int NeighborPointIterator = 0;
                NeighborPointIterator < static_cast<unsigned int>(PreviousPointConnectivitiesOnSurface->GetNumberOfIds());
                NeighborPointIterator++)
            {
                // Get Neighbor point Id on surface
                vtkIdType NeighborPointSurfaceId = PreviousPointConnectivitiesOnSurface->GetId(NeighborPointIterator);
                vtkIdType NeighborPointOriginalId = SurfacePointsOriginalIdsArray->GetTuple1(NeighborPointSurfaceId);

                // Current point original Id
                vtkIdType CurrentEdgePointOriginalId = EdgePointsOriginalIdsArray->GetTuple1(CurrentEdgePointId);

                // Check if a neighbor of previous point is the same of current point
                if(CurrentEdgePointOriginalId == NeighborPointOriginalId)
                {
                    CurrentEdgePointIdOnSurface = NeighborPointSurfaceId;
                    EdgePointsSurfaceIds->SetTuple1(CurrentEdgePointId,CurrentEdgePointIdOnSurface);
                    SurfaceIdFound = true;
                    break;
                }
            }

            // Check if surface point no found
            if(SurfaceIdFound == false)
            {
                vtkErrorMacro("Could not match original Ids of current marching point with surface neighbor points.");
                return SurfaceIdFound;
            }

            // Update Marching points
            vtkIdType NextEdgePointId = this->FindAdjacentEdgePoint(CurrentEdgePointId,PreviousEdgePointId,Edges);
            PreviousEdgePointId = CurrentEdgePointId;
            PreviousEdgePointIdOnSurface = CurrentEdgePointIdOnSurface;
            CurrentEdgePointId = NextEdgePointId;
        }
    }

    return SurfaceIdFound;
}

// ===========================================
// Add Id Conversion Arrays To Separated Edges
// ===========================================

// Description:
// This method adds the following two vtkIdTypeArray to the PointData of SeparatedEdges
// 1- OriginalPointIds: Equivalent Ids of edge points in original InputDataSet:
// 2- EdgePointsSurfaceIds: Equivalent Ids of edge points in InputDataSetSurface

void FeatureSurface::AddIdConversionArraysToSeparatedEdges(
        vtkPolyData *Edges,
        vtkPolyData *SeparatedEdges,
        vtkPolyData *InputDataSetSurface)
{
    // 1- Adding Original Point Ids //

    // Adding Original points Ids from DataSet to Edge
    vtkSmartPointer<vtkIdTypeArray> OriginalPointIds = 
        vtkIdTypeArray::SafeDownCast(Edges->GetPointData()->GetArray("vtkOriginalPointIds"));
    if(!OriginalPointIds)
    {
        vtkErrorMacro("Edges polydata does not have vtkOriginalPointIds field.");
    }
    SeparatedEdges->GetPointData()->AddArray(OriginalPointIds);

    // 2- Adding Surface Point Ids //

    // Convert Edge Ids to Surface Ids
    vtkSmartPointer<vtkIdTypeArray> EdgePointsSurfaceIds = vtkSmartPointer<vtkIdTypeArray>::New();
    bool SurfaceIdFound = this->ConvertEdgeIdsToSurfaceIds(
            Edges,
            SeparatedEdges,
            InputDataSetSurface,
            EdgePointsSurfaceIds);

    // Check success of finding Ids
    if(!SurfaceIdFound)
    {
        vtkErrorMacro("Error in converting edge Ids to surface Ids occured.");
    }

    // Check Array
    if(!EdgePointsSurfaceIds)
    {
        vtkErrorMacro("EdgePointSurfaceIds is NULL.");
    }
    else if(EdgePointsSurfaceIds->GetNumberOfTuples() != SeparatedEdges->GetNumberOfPoints())
    {
        vtkErrorMacro("Number of points does not match.");
    }

    // Add array to SeparatedEdge
    SeparatedEdges->GetPointData()->AddArray(EdgePointsSurfaceIds);
}

// ================
// Original Id Test
// ================

bool FeatureSurface::OriginalIdTest(
        vtkDataSet *OriginalDataSet,
        vtkPolyData *ExtractedDataSet)
{
    // Original Ids array
    vtkSmartPointer<vtkIdTypeArray> OriginalPointIdsArray =
        vtkIdTypeArray::SafeDownCast(ExtractedDataSet->GetPointData()->GetArray("vtkOriginalPointIds"));
    if(!OriginalPointIdsArray)
    {
        vtkErrorMacro("OriginalPointIdsArray is NULL.");
    }

    // Number Of Points in Extracted PolyData
    unsigned int NumberOfPoints = ExtractedDataSet->GetNumberOfPoints();

    // Loop over points
    bool IdsMatched = true;
    for(unsigned int PointId = 0; PointId < NumberOfPoints; PointId++)
    {
        // Extracted Point Coordinates
        double ExtractedPointCoordinates[3];
        ExtractedDataSet->GetPoint(PointId,ExtractedPointCoordinates);

        // Original Point Id
        vtkIdType OriginalPointId = OriginalPointIdsArray->GetTuple1(PointId);

        // Get Original Point Coordinates
        double OriginalPointCoordinates[3];
        OriginalDataSet->GetPoint(OriginalPointId,OriginalPointCoordinates);

        // Difference of points
        double Difference[3];
        vtkMath::Subtract(ExtractedPointCoordinates,OriginalPointCoordinates,Difference);
        double Error = vtkMath::Norm(Difference);
        
        if(Error > 1e-8)
        {
            IdsMatched = false;
            vtkErrorMacro("Testing original point Ids: some Ids do not match.");
        }
    }

    if(IdsMatched)
    {
        std::cout << "All original Ids matched successfully." << std::endl;
    }

    return IdsMatched;
}

// ================================
// Convert Edge To Surface Ids Test
// ================================

bool FeatureSurface::ConvertEdgeToSurfaceIdsTest(
        vtkPolyData *SeparatedEdges,
        vtkPolyData *InputDataSetSurface)
{
    // Status of output Test
    bool IdsMatched = true;

    // Edge Point Original Ids
    vtkSmartPointer<vtkIdTypeArray> EdgePointsOriginalIds =
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalPointIds"));

    if(!EdgePointsOriginalIds)
    {
        vtkErrorMacro("EdgePointOriginalIds is NULL.");
    }

    // Edge Points Surface Ids
    vtkSmartPointer<vtkIdTypeArray> EdgePointsSurfaceIds = 
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalSurfaceIds"));

    if(!EdgePointsSurfaceIds)
    {
        vtkErrorMacro("EdgePointsSurfaceIds is NULL.");
    }

    // Surface Point Original Ids
    vtkSmartPointer<vtkIdTypeArray> SurfacePointsOriginalIds =
        vtkIdTypeArray::SafeDownCast(InputDataSetSurface->GetPointData()->GetArray("vtkOriginalPointIds"));

    if(!SurfacePointsOriginalIds)
    {
        vtkErrorMacro("SurfacePointsOriginalIds is NULL.");
    }

    // Iterate over Separated edge points
    for(unsigned int EdgePointIterator = 0;
        EdgePointIterator < static_cast<unsigned int>(SeparatedEdges->GetNumberOfPoints());
        EdgePointIterator++)
    {
        // Get Edge Point Original Id
        vtkIdType EdgePointOriginalId = EdgePointsOriginalIds->GetTuple1(EdgePointIterator);

        // Get Surface Id of edge point
        vtkIdType EdgePointSurfaceId = EdgePointsSurfaceIds->GetTuple1(EdgePointIterator);

        // Skip Non Closed-Edge points with negative surface Id
        if(EdgePointSurfaceId < 0)
        {
            continue;
        }

        // Get Surface point Original Id
        vtkIdType SurfacePointOriginalId = SurfacePointsOriginalIds->GetTuple1(EdgePointSurfaceId);

        // Compare Surface and Edge original Ids
        if(EdgePointOriginalId != SurfacePointOriginalId)
        {
            std::cout << "Edge and Surface original Ids not match. edge point: " << EdgePointIterator << std::endl;
            IdsMatched = false;
        }
    }

    return IdsMatched;
}

// ================
// Colorize Surface
// ================

// Description:
// Adds a color to each surface points as following
// 1- Non- featured surface (body of surface) are zero
// 2- All points of each featres surface are the same color, the color are integers bigger than zero.
//    For example, All points inside the first closed edge colored 1, the second surface colored 2, etc.
// 3- All edge points are colored the same color of their surface, but in negative integer.
//    For example, the first closed edge is -1, the second polygon is -2, etc.

void FeatureSurface::ColorizeSurface(
        vtkPolyData *InputDataSetSurface,
        vtkPolyData *SeparatedEdges,
        vtkCellArray *FeatureSurfaceVertices,
        vtkIntArray *SurfaceColors)
{
    // Initialize output
    SurfaceColors->SetNumberOfTuples(InputDataSetSurface->GetNumberOfPoints());
    SurfaceColors->SetNumberOfComponents(1);
    SurfaceColors->SetName("SurfaceColors");

    // 1- Add zero background color to all points //
    
    // Define body color
    int BodyColorValue = static_cast<int>(0);

    // Iterate over all points
    for(unsigned int PointIterator = 0;
        PointIterator < static_cast<unsigned int>(InputDataSetSurface->GetNumberOfPoints());
        PointIterator++)
    {
        SurfaceColors->SetTuple1(PointIterator,BodyColorValue);
    }

    // 2- Colorize Edges //

    // Edge point Ids on Surface
    vtkSmartPointer<vtkIdTypeArray> EdgePointsSurfaceIds =
        vtkIdTypeArray::SafeDownCast(SeparatedEdges->GetPointData()->GetArray("vtkOriginalSurfaceIds"));

    if(!EdgePointsSurfaceIds)
    {
        vtkErrorMacro("EdgePointsSurfaceIDs is NULL.");
    }

    // Array of edges
    vtkSmartPointer<vtkCellArray> ClosedEdgesArray = SeparatedEdges->GetPolys();
    ClosedEdgesArray->InitTraversal();

    // Iterate over edge polygons
    for(unsigned int ClosedEdgeIterator = 0;
        ClosedEdgeIterator < static_cast<unsigned int>(ClosedEdgesArray->GetNumberOfCells());
        ClosedEdgeIterator++)
    {
        // Get a closed edge polygon
        vtkSmartPointer<vtkIdList> ClosedEdge = vtkSmartPointer<vtkIdList>::New();
        ClosedEdgesArray->GetNextCell(ClosedEdge);

        // Color of Edge group
        int EdgeColorValue = -1;
        switch(this->EdgeColorMode)
        {
            // Body color
            case EdgeColor::BODY_COLOR:
            {
                EdgeColorValue = BodyColorValue;
                break;
            }

            // Mono color
            case EdgeColor::MONO_COLOR:
            {
                EdgeColorValue = static_cast<int>(-1);
                break;
            }

            // Surface color
            case EdgeColor::SURFACE_COLOR:
            {
                EdgeColorValue = static_cast<int>(ClosedEdgeIterator)+1;
                break;
            }

            // Opposite surface color
            case EdgeColor::OPPOSITE_SURFACE_COLOR:
            {
                EdgeColorValue = -static_cast<int>(ClosedEdgeIterator)-1;
                break;
            }

            // Undefined mode
            default:
            {
                vtkErrorMacro("Edge color mode is undefined.");
            }
        }

        // Iterate over points on each edge polygon
        for(unsigned int EdgePointIterator = 0;
            EdgePointIterator < static_cast<unsigned int>(ClosedEdge->GetNumberOfIds());
            EdgePointIterator++)
        {
            // Get edge point Id
            vtkIdType EdgePointId = ClosedEdge->GetId(EdgePointIterator);

            // Get edge point Id on Surface
            vtkIdType EdgePointSurfaceId = EdgePointsSurfaceIds->GetTuple1(EdgePointId);

            // Final check if the point was not colored before
            if(SurfaceColors->GetTuple1(EdgePointSurfaceId) != 0)
            {
                vtkErrorMacro("Edge point was belonging to two or more edge groups.");
            }

            // Colorize edge point
            SurfaceColors->SetTuple1(EdgePointSurfaceId,EdgeColorValue);
        }
    }

    // 3- Colorize Feature Surfaces //

    // Iterate over feature surfaces
    FeatureSurfaceVertices->InitTraversal();

    for(unsigned int FeatureSurfaceIterator = 0;
        FeatureSurfaceIterator < static_cast<unsigned int>(FeatureSurfaceVertices->GetNumberOfCells());
        FeatureSurfaceIterator++)
    {
        // Define color of surface
        int SurfaceColorValue = 1;

        switch(this->SurfaceColorMode)
        {
            // Body color
            case SurfaceColor::BODY_COLOR:
            {
                SurfaceColorValue = BodyColorValue;
                break;
            }

            // Mono color
            case SurfaceColor::MONO_COLOR:
            {
                SurfaceColorValue = static_cast<int>(1);
                break;
            }

            // Poly Color
            case SurfaceColor::POLY_COLOR:
            {
                SurfaceColorValue = static_cast<int>(FeatureSurfaceIterator)+1;
                break;
            }

            // Undefined mode
            default:
            {
                vtkErrorMacro("Surface color mode is undefined.");
            }
        }

        // Get a vertex cell
        vtkSmartPointer<vtkIdList> ClosedEdgeSurface = vtkSmartPointer<vtkIdList>::New();
        FeatureSurfaceVertices->GetNextCell(ClosedEdgeSurface);

        // Iterate over points
        for(unsigned int SurfacePointIterator = 0;
            SurfacePointIterator < static_cast<unsigned int>(ClosedEdgeSurface->GetNumberOfIds());
            SurfacePointIterator++)
        {
            // Get a surface point
            vtkIdType SurfacePointId = ClosedEdgeSurface->GetId(SurfacePointIterator);

            // Final check if the point was not colored before
            if(SurfaceColors->GetTuple1(SurfacePointId) != 0)
            {
                vtkErrorMacro("Surface point was belonging to two or more edge groups.");
            }

            // Colorize surface point
            SurfaceColors->SetTuple1(SurfacePointId,SurfaceColorValue);
        }
    }
}

// ===============================
// Generate Output DataSet Surface
// ===============================

void FeatureSurface::GenerateOutputDataSetSurface(
        vtkPolyData *InputDataSetSurface,
        vtkPolyData *SeparatedEdges,
        vtkCellArray *FeatureSurfaceVertices,
        vtkIntArray *SurfaceColors,
        vtkPolyData *OutputDataSetSurface)
{
    // 1- Copy Struture of PolyData from Input to Output //
    OutputDataSetSurface->CopyStructure(InputDataSetSurface);

    // 2- Copy Original Point Ids array //
    vtkSmartPointer<vtkIdTypeArray> OriginalPointIdsArray = 
        vtkIdTypeArray::SafeDownCast(InputDataSetSurface->GetPointData()->GetArray("vtkOriginalPointIds"));

    if(!OriginalPointIdsArray)
    {
        vtkErrorMacro("OriginalPointIdsArray is NULL.");
    }

    OutputDataSetSurface->GetPointData()->AddArray(OriginalPointIdsArray);

    // 3- Add Surface Colors array //
    if(!SurfaceColors)
    {
        vtkErrorMacro("SurfaceColors array is NULL.");
    }
    else if(SurfaceColors->GetNumberOfTuples() != InputDataSetSurface->GetNumberOfPoints())
    {
        vtkErrorMacro("SurfaceColors array does not have correct number of points.");
    }

    OutputDataSetSurface->GetPointData()->SetScalars(SurfaceColors);

    // 4- Add Feature Surface vertices //
    // OutputDataSetSurface->SetVerts(FeatureSurfaceVertices);

    // 5- Add Separated Edges into Polygons //
    vtkSmartPointer<vtkCellArray> ClosedEdgesArray = SeparatedEdges->GetPolys();
    // OutputDataSetSurface->SetPolys(ClosedEdgesArray); // CHANGED
    // OutputDataSetSurface->SetVerts(ClosedEdgesArray);
}

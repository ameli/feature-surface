/*
 * =====================================================================================
 *
 *       Filename:  FeatureSurface.h
 *
 *    Description:  Extract Topologically distinct surfaces form data set
 *
 *        Version:  1.0
 *        Created:  12/14/2013 03:47:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Siavash Ameli
 *   Organization:  University of California, Berkeley
 *
 * =====================================================================================
 */

#ifndef __FeatureSurface_h
#define __FeatureSurface_h

// ======
// Macros
// ======

// Debugging
#define HERE std::cout << "\033[0;44m" << "DEBUG: " << __FILE__ << " at " << __LINE__ << "\033[0m" << std::endl;

// =====
// Types
// =====

// Surface color
class SurfaceColor
{
    public:
        enum SurfaceColorModeType
        {
            BODY_COLOR = 0,
            MONO_COLOR,
            POLY_COLOR,
            NumberOfSurfaceColorModes
        };
};

// Edge color
class EdgeColor
{
    public:
        enum EdgeColorModeType
        {
            BODY_COLOR = 0,
            MONO_COLOR,
            SURFACE_COLOR,
            OPPOSITE_SURFACE_COLOR,
            NumberOfEdgeColorModes
        };
};

// =====================
// Foreward Declarations
// =====================

#include <vtkDataSetAlgorithm.h>
#include <vector>

class vtkPolygon;
class vtkIdList;
class vtkIdTypeArray;
class vtkIntArray;
class vtkCellArray;

// =====================
// Feature Surface Class
// =====================

class FeatureSurface: public vtkDataSetAlgorithm
{
    public:
        static FeatureSurface *New();
        vtkTypeRevisionMacro(FeatureSurface,vtkDataSetAlgorithm);
        virtual void PrintSelf(ostream &os,vtkIndent indent);

        // Accessors / Mutators
        vtkGetMacro(DebugStatus,bool);
        vtkSetMacro(DebugStatus,bool);
        void DebugOn() { this->DebugStatus = true; };
        void DebugOff() { this->DebugStatus = false; };

        vtkGetMacro(FeatureAngle,double);
        vtkSetMacro(FeatureAngle,double);

        vtkGetMacro(SurfaceColorMode,SurfaceColor::SurfaceColorModeType);
        void SetSurfaceColorMode(int SurfaceColorInt);

        vtkGetMacro(EdgeColorMode,EdgeColor::EdgeColorModeType);
        void SetEdgeColorMode(int EdgeColorInt);

    protected:
        FeatureSurface();
        virtual ~FeatureSurface();

        // Pipeline Executives
        virtual int FillOutputPortInformation(
                int port,
                vtkInformation *info);
        
        virtual int RequestData(
                vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector);

        // Methods
        unsigned int FindSeparatedEdges(
                vtkPolyData *Edges,
                vtkPolyData *SeparatedEdges,
                vtkDataSet *InputDataSet);

        bool FindClosedPolygon(
                vtkIdType InquiryPointId,
                vtkPolyData *Edges,
                bool *PointsProcessed,
                vtkPolygon *Polygon);

        int FindAdjacentEdgePoint(
                vtkIdType CurrentPointId,
                vtkIdType PreviousAdjacentPointId,
                vtkPolyData *Edges);

        void GetPointConnectivities(
                vtkIdType InquiryPointId,
                vtkDataSet *DataSet,
                vtkIdList *Connectivities);

        void FindSeparatedEdgeSurfaces(
                vtkDataSet *InputDataSet,
                vtkPolyData *InputDataSetSurface,
                vtkPolyData *Edges,
                vtkPolyData *SeparatedEdges,
                vtkCellArray *FeatureSurfaceVertices);

        void FindGeometricCenterOfClosedEdge(
                vtkDataSet *InputDataSet,
                vtkPolyData *SeparatedEdges,
                vtkIdList *ClosedEdge,
                double GeometricCenter[3]);

        void FindACandidatePointOnEdgeSurface(
                vtkDataSet *InputDataSet,
                vtkPolyData *InputDataSetSurface,
                vtkPolyData *Edges,
                vtkPolyData *SeparatedEdges,
                vtkIdList *ClosedEdge,
                double GeometricCenterOfClosedEdge[3],
                vtkIdType &ACandidatePointOnEdgeSurfaceSurfacelId);

        bool GetEdgePointConnectivitiesOnSurfaceExcludingEdge(
                vtkIdType EdgePointId,
                vtkIdList *EdgePointConnectivitiesOnEdge,
                vtkPolyData *InputDataSetSurface,
                vtkPolyData *SeparatedEdges,
                vtkIdList *EdgePointConnectivitiesOnSurfaceExcludingEdge);

        void ConvertSurfaceIdListToOriginalIdList(
                vtkIdList *SurfaceIdList,
                vtkPolyData *InputDataSetSurface,
                vtkIdList *OriginalIdList);

        double GetRayAngleProjectionOnSurface(
                double EdgePointCoordinates[3],
                double AdjacentOfEdgePointCoordinates[3],
                double GeometricCenterOfClosedPoint[3],
                double NeighborOfEdgePointCoordinates[3]);

        void FindClosedEdgeSurface(
                vtkPolyData *InputDataSetSurface,
                vtkIdList *ClosedEdge,
                vtkPolyData *SeparatedEdges,
                vtkIdType ACandidatePointOnEdgeSurfaceId,
                vtkIdList *ClosedEdgeSurface);

        int IdConverter(
                vtkDataSet *SourceDataSet,
                vtkDataSet *DestinationDataSet,
                vtkIdType SourceId);

        void FindOriginalPointIds(
                vtkPolyData *ExtractedDataSet,
                vtkDataSet *OriginalDataSet,
                vtkIdTypeArray *OriginalCellIds,
                vtkIdTypeArray *OriginalPointIdsy);

        bool ConvertEdgeIdsToSurfaceIds(
                vtkPolyData *Edge,
                vtkPolyData *SeparatedEdges,
                vtkPolyData *InputDataSetSurface,
                vtkIdTypeArray *EdgePointsSurfaceIds);

        void AddIdConversionArraysToSeparatedEdges(
                vtkPolyData *Edges,
                vtkPolyData *SeparatedEdges,
                vtkPolyData *InputDataSetSurface);

        bool OriginalIdTest(
                vtkDataSet *OriginalDataSet,
                vtkPolyData *ExtractedDataSet);

        bool ConvertEdgeToSurfaceIdsTest(
                vtkPolyData *SeparatedEdges,
                vtkPolyData *InputDataSetSurface);

        void ColorizeSurface(
                vtkPolyData *InputDataSetSurface,
                vtkPolyData *SeparatedEdges,
                vtkCellArray *FeatureSurfaceVertices,
                vtkIntArray *SurfaceColors);

        void GenerateOutputDataSetSurface(
                vtkPolyData *InputDataSetSurface,
                vtkPolyData *SeparatedEdges,
                vtkCellArray *FeatureSurfaceVertices,
                vtkIntArray *SurfaceColors,
                vtkPolyData *OutputDataSetSurface);

        // Member Data
        bool DebugStatus;
        double FeatureAngle;
        SurfaceColor::SurfaceColorModeType SurfaceColorMode;
        EdgeColor::EdgeColorModeType EdgeColorMode;

    private:
        FeatureSurface(const FeatureSurface &);   // Not implemented.
        void operator=(const FeatureSurface &);   // Not implemented.
};

// ==========
// Prototypes
// ==========

size_t GetMaxIndexOfVector(std::vector<double> const &Vector);

#endif

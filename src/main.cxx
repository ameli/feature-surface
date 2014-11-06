/*
 * =====================================================================================
 *
 *       Filename:  TestFeatureSurface.cxx
 *
 *    Description:  Test for FeatureSurface
 *
 *        Version:  1.0
 *        Created:  12/14/2013 04:15:04 PM
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

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include "FeatureSurface.h"
#include <vtkPolyDataWriter.h>

// ====
// Main
// ====

int main(int argc, char *argv[])
{
    // Inuput filenames
    char InputFilename[256];
    if(argc > 1)
    {
        strcpy(InputFilename,argv[1]);
    }
    else
    {
        strcpy(InputFilename,"/home/sia/data/vtkfiles/restart_vel.2800.vtk");
    }

    // Output filename
    char OutputFilename[256];
    if(argc > 1)
    {
        strcpy(OutputFilename,argv[2]);
    }
    else
    {
        strcpy(OutputFilename,"/home/sia/Desktop/output_DataSetSurface.vtk");
    }

    // Read Data
    vtkSmartPointer<vtkUnstructuredGridReader> Reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    // vtkSmartPointer<vtkXMLUnstructuredGridReader> Reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    Reader->SetFileName(InputFilename);
    Reader->Update();

    // Feature Surface
    vtkSmartPointer<FeatureSurface> FeatureSurfaceFilter = vtkSmartPointer<FeatureSurface>::New();
    FeatureSurfaceFilter->SetInputConnection(Reader->GetOutputPort());
    FeatureSurfaceFilter->DebugOn();
    // FeatureSurfaceFilter->SetSurfaceColorMode(1);
    FeatureSurfaceFilter->Update();

    // Write Data
    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    Writer->SetInputConnection(FeatureSurfaceFilter->GetOutputPort());
    Writer->SetFileName(OutputFilename);
    Writer->Update();
    
    return EXIT_SUCCESS;
}

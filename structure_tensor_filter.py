# -*- coding: utf-8 -*-
"""Custom paraview filter to map line orientation to cells."""
import numpy as np
import vtk
from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid


@smproxy.filter(label="Structure Tensor")
@smproperty.input(name="Cells", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class StructureTensorFilter(VTKPythonAlgorithmBase):
    """StructureTensorFilter.

    Attributes
    ----------
    _radius : float
        Search radius

    """

    def __init__(self):
        """Set up the filter."""
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            nOutputPorts=1,
            outputType="vtkUnstructuredGrid",
        )
        self._radius = 1.0

    def FillInputPortInformation(self, port, info):
        """Require unstructured grids as input types."""
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    @smproperty.doublevector(name="Radius", default_values=[1.0])
    @smdomain.doublerange()
    def SetRadius(self, radius):
        """Create an input field for the search radius and save it."""
        if self._radius != radius:
            self._radius = radius
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        """Process the request submitted after applying the filter."""
        # Access and wrap input source
        source_input = vtkUnstructuredGrid.GetData(inInfo[0])
        source = dsa.WrapDataObject(source_input)
        N = source.GetNumberOfCells()

        # Create a vtkKdTreePointLocator for fast computation of neighborhood
        point_locator = vtk.vtkKdTreePointLocator()
        point_locator.SetDataSet(source_input)
        point_locator.BuildLocator()

        # Set up result arrays
        n = np.zeros(N)
        A = np.nan * np.ones((N, 3, 3))

        # Iterate over cells
        for i in range(N):
            # Find neighbors
            idList = vtk.vtkIdList()
            pos = np.array([0.0, 0.0, 0.0])
            cell = source.GetCell(i)
            cell.ComputeBoundingSphere(pos)
            point_locator.FindPointsWithinRadius(self._radius, pos, idList)

            # Evaluate neighbors
            N0 = idList.GetNumberOfIds()
            ids = [idList.GetId(j) for j in range(N0)]
            points = source.Points[ids]
            dist = points - pos
            dist_norm = np.linalg.norm(dist, axis=1)
            dist = dist / dist_norm
            w = np.ones_like(dist_norm)

            # Assign number of neighbors
            n[i] = N0

            # Assign structure tensors
            if N0 > 0:
                A[i, :, :] = np.einsum("k, ki, kj->ij", w, dist, dist) / np.sum(w)

        # Create output
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))
        output.ShallowCopy(source_input)
        output.CellData.append(n, "N")
        output.CellData.append(A, "Structure Tensor")

        return 1

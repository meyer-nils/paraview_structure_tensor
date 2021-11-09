[![LICENSE](https://black.readthedocs.io/en/stable/_static/license.svg)](https://raw.github.com/nilsmeyerkit/paraview_structure_tensor/master/LICENSE)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Paraview Structure Tensors
This repository defines a ParaView filter that builds a structure tensor from nodes around a cell center within a sphere. This structure tensors might help to identify the geometric configuration. If the sphere is completely filled, the result tends an isotropic state. If the geometry is plate-like, the result tends toward a transversal isotropic state. If the geometry is almost one-dimensional, it tends toward a unidirectional state.

The plugin requires Paraview 5.8 or higher.

Load the plugin to Paraview via 'Tools' -> 'Manage Plugins...' -> 'Load New'.

# Use
Load a .vtk file with cells to ParaView. Execute 'Filters' -> 'Structure Tensor' and select the cells. Click 'Apply' to start the plugin and 'Reload Python Module' to update the filter to changes in the code.

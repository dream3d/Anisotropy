Steiner Compact {#steinercompact}
======

## Group (Subgroup) ##
Unsupported (Anisotropy)

## Description ##
This **Filter** estimates the morphological anisotropy of the grains on the basis of their profiles in the section planes. Normal to the section planes is selected to be perpendicular to one of the sample sites. Typically, anisotropy is analysed in the plane of mapping (with normal parallel to Z-direction). 

Steiner compact is estimated by counting number of intersections of random test lines in specified directions with the grain boundaries (rose of intersections). This **Filter** creates a .vtk file which can be visualized in Paraview. Steiner compacts for different phases are distinguished by different Z-coordinate in the output file.

The user can set the (maximum) number of sites of the Steiner compact. This corresponds to the number of directions for which the hits of the grain boundaries by random test lines are counted. Note that the actual number of sites can be lower since some of the sites can have zero length.


## Parameters ##
| Name | Type | Description |
|------|------| ----------- |
| Section Plane | int | Section plane where the Steiner compact is found. |
| Number Of Sites | int | Number of sites of the Steiner compact. |
| Output File | File Path | The output .vtk file path where the Steiner compact is drawn. |

## Required Geometry ##
Image 

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Cell Attribute Array** | FeatureIds | int32_t | (1) | Specifies to which **Feature** each **Cell** belongs |
| **Cell Attribute Array** | Phases | int32_t | (1) | Specifies to which **Ensemble** each **Cell** belongs. |

## Created Objects ##
None

## References ##

Journal articles on _Steiner Compact_ that are useful:

+ _Analysis of planar anisotropy by means of the Steiner compact_. J. Appl. Prob. 26:490-502. Rataj, J. and Saxl, I. (1989).
+ _Planar anisotropy revisited_. Kybernetika 36(2): 149-164. Benes, V. and Gokhale, A.M. (2000). 
+ _Non-parametric estimation of the directional distribution of stationary line and fibre processes_. Adv. in Appl. Prob. 33(1): 6-24. Kiderlen, M. (2001).

## License & Copyright ##

Please see the description file distributed with this **Plugin**

## DREAM.3D Mailing Lists ##

If you need more help with a **Filter**, please consider asking your question on the [DREAM.3D Users Google group!](https://groups.google.com/forum/?hl=en#!forum/dream3d-users)



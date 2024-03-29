Adaptive Alignment (Misorientation) 
======

## Group (Subgroup) ##

Unsupported (Anisotropy)

## Description ##

This **Filter** attempts to align consecutive 'sections' perpendicular to the Z-direction of the sample by determining the position that results in the minimum amount of misorientation between **Cells** directly above-below each other.

This initial alignment is followed by a correction algorithm which adapts the shifts to a complementary criterion having one of the following forms:

1. Images of the mapped area taken from each slice by scanning electron microscope (EDAX).
2. User-defined shifts between the first and the last slice of the stack.

The initial alignment algorithm of this **Filter** is as follows:

1. Calculate the misorientation between each **Cell** in a section and the **Cell** directly above it in the next section  
2. Count the number of **Cell** pairs that have a misorientation above the user defined tolerance and store that as the misalignment value for that position
3. Repeat steps 1 and 2 for each position when shifting the second slice (relative to the first) from three (3) **Cells** to the left to three (3) **Cells** to the right, as well as from three (3) **Cells** up to three (3) **Cells** down. *Note that this creates a 7x7 grid*
4. Determine the position in the 7x7 grid that has the lowest misalignment value. (It will be the position with the fewest different **Cell** pairs)
5. Repeat steps 1-4 with the center of each (new) 7x7 grid at the best position from the last 7x7 grid until the best position in the current/new 7x7 grid is the same as the last 7x7 grid
6) Repeat steps 1-5 for each pair of neighboring sections

**Note that this is similar to a downhill simplex and can get caught in a local minimum!**

The correction alignment algorithm of this **Filter** attempts to improve the complementary fit as follows:

1. Start with the shifts obtained by the initial algorithm, i.e., define the current shifts as those obtained from the initial algorithm for each pair of consecutive cross sections.
2. For each pair of consecutive cross sections, consider all the shifts checked within the alignment algorithm (7x7 grids), different from the current shift, which improve the complementary fit. From these shifts, find the candidate shift as that with the lowest misalignment value. 
3. From the candidate shifts identified in step 2 (one candidate shift for each pair of consecutive sections), select the one with minimum difference between the misalignment value of the candidate shift and the misalignment value of the current shift. Change the current state by this candidate shift and recompute the complementary fit.
4. Repeat steps 2 - 3 until no improvement in the complementary fit is attained in step 2.

If the user elects to use a mask array, the **Cells** flagged as *false* in the mask array will not be considered during the alignment process.  

The user can choose to write the determined shift to an output file by enabling *Write Alignment Shifts File* and providing a file path.  


## Parameters ##

| Name | Type | Description |
|------|------| ----------- |
| Misorientation Tolerance | float | Tolerance used to decide if **Cells** above/below one another should be considered to be _the same_. The value selected should be similar to the tolerance one would use to define **Features** (i.e., 2-10 degrees) |
| Write Alignment Shift File | bool | Whether to write the shifts applied to each section to a file |
| Alignment File | File Path | The output file path where the user would like the shifts applied to the section to be written. Only needed if *Write Alignment Shifts File* is checked |
| Global Correction: SEM Images | bool | Whether to use SEM images for adaptive alignment. |
| Global Correction: Own Shifts | bool | Whether to set the shifts for adaptive alignment manually. |
| Total Shift In X-Direction (Microns) | float | Shift in X-direction between the first and the last slice of the stack in microns. Only needed if *Global Correction: Own Shifts* is checked. |
| Total Shift In Y-Direction (Microns) | float | Shift in Y-direction between the first and the last slice of the stack in microns. Only needed if *Global Correction: Own Shifts* is checked. |
| Use Mask Array | bool | Whether to remove some **Cells** from consideration in the alignment process |
 
## Required Geometry ##

Image 

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Image Data Array** | ImageData | uint8_t | any component size | Specifies image data containing the SEM images to be used for global correction of the shifts. |
| **Cell Attribute Array** | Quats | float | (4) | Specifies the orientation of the **Cell** in quaternion representation. |
| **Cell Attribute Array** | Phases | int32_t | (1) | Specifies to which **Ensemble** each **Cell** belongs. |
| **Cell Attribute Array** | Mask | bool | (1) | Specifies if the **Cell** is to be counted in the algorithm. Only required if *Use Mask Array* is checked. |
| **Ensemble Attribute Array** | CrystalStructures | uint32_t | (1) | Enumeration representing the crystal structure for each **Ensemble**. |

## Created Objects ##

None

## Example Pipelines ##

+ 02_Adaptive Alignment - Misorientation - Zero Shifts

## License & Copyright ##

Please see the description file distributed with this **Plugin**

## DREAM.3D Mailing Lists ##

If you need more help with a **Filter**, please consider asking your question on the [DREAM.3D Users Google group!](https://groups.google.com/forum/?hl=en#!forum/dream3d-users)



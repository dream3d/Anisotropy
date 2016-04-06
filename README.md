# Anisotropy
Anisotropy is a DREAM.3D plugin which contains filters that help to analyse an anisotropy in morphology of grains and correctly align 2D slices which build up a 3D data set of anisotropic material. 

Anisotropy is analysed via Steiner compact which is estimated by counting the intersections of grain boundaries with random lines in systematically selected directions. As a result, a polygon (in 2D) or polyhedra (in 3D) is obtained which approximates the average grain shape.

For correct alignment, a complementary information is utilized, coming from either lower-resolution SEM images (provided by EDAX microscope) where the mapped area is marked by a green rectangle, or user-defined shifts in X- and Y- direction between the first and last slice of the stack. If the user-defined shifts are set to zero, no systematic misalignment between the slices is expected and the slices are aligned accordingly.

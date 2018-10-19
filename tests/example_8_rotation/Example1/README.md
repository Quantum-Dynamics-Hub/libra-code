# Instructions 

  1.  Just run the script ```python test_rotate_frag.py ```

  2.  Load any of the produced output .xyz files to the VMD to visualize

  3.  Optionally, create a movie:

     * load all the files into VMD: in Tcl console use: ```mol addfile <filename>```
     * in the VMD use: Extensions -> Visualization -> Movie Maker 
       Select: Renderer -> Snapshot; Movie Settings -> Trajectory; uncheck "Delete files"
       Use the trajectory timestep of 4 
     * using the ffmpeg.ext as:
     * ```ffmpeg.exe -i untitled.%05d.bmp -vframes 400 -crf 25 demo.avi``` will generate a movie


# Explanations 

  Here, we have a system composed of C60 (first 60 atoms) and graphene sheet (the rest of atoms)

  1. We can group arbitrary sets of atoms to form a rigid body for the following 
     rotations and translations of them. This is done my the instructions:

     ``` syst0.GROUP_ATOMS(range(1,61), 1)```

    Note that the numbers used here are assumed to start from 1 (not from 0 as conventional in programming).
    This is because, the numbers are regarded as more intuitive atom or group ID numbers, so starting from 1

  2. We can now perform some manipulations with such fragments:

    *  **test_rotations()** function shows how to rotate the fragments

       The main function utilized here is: ```syst0.ROTATE_FRAGMENT(dphi, u, frag_id,  pivot_pt)```
       Here:

       ** dphi - a small angle by which the fragment is rotated at every function call (the effect 
        is cumulative)
       ** u - is a vector defining the axis around which the rotation occurs
       ** frag_id - is an ID of the fragment that we manipulate. Again, the numbering starts from 1 here.
       In fact, this is an identifier that you assign to a fragment when you group atoms
       ** pivot_pt - the vector defining the center of rotation


    *  **test_translations()** function shows how to translate the fragments

       The main function utilized here is: ```syst0.TRANSLATE_FRAGMENT(dr, tdir, frag_id)```
       Here:

       ** dr - a small displacement by which the fragment is moved at every function call (the effect 
        is cumulative)
       ** tdir - is a vector defining the direction of the fragment's translation
       ** frag_id - is an ID of the fragment that we manipulate. Again, the numbering starts from 1 here.
       In fact, this is an identifier that you assign to a fragment when you group atoms


    *  **test_mixed()** function shows how to rotate and translate the fragments at the same time
        in this case, we can rotate the fragments around the displaced center of mass, for instance


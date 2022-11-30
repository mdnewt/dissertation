Iterative Shape Averaging – Help Guides

Average Creation Step 1: ISA.M

To create the average atlas successfully, the following steps must be
adhered to:
 
1) Create a main directory and save ISA_1.m to this directory. Then,
   create the following sub-directories within this directory:
 
  '0 - Pre-training Images' - OPTIONAL - In cases in which the user
                              desires to resample the raw training images
                              to a lower resolution prior to drawing the
                              training masks, the raw training images can
                              be placed here. The code *WHAT* can then be
                              used to create training images from this
                              list of pre-training images.
 
  '1 - Training Images'     - This will contain .MATs of all image stacks
                              to be used to create the average atlas. If
                              the tissue of interest has a laterality (ie
                              right and left femora), these should be
                              mirrored as necessary to produce an
                              identical image set. The convention is that
                              lefts should be mirrored to look like
                              rights.
  
  '2 - Training Masks'      - This will contain SUBFOLDERS for each
                              tissue included in the average atlas. Each
                              subfolder will then contain .MATs of the
                              corresponding tissue for each training
                              image.
 
  '3 - Preprocess'          - Preprocessing is performed when necessary
                              to resample the image stack and/or masks to
                              the resolution at which they will be
                              averaged.
 
  '4 - Unknown Samples'     - This will contain .MATs of the unknown
                              images upon which the average atlas will be
                              co-registered.
 
  '5 - Registrations'       - Registrations of the average atlas onto
                              unknown samples will be saved here.
 
2) Create a .mat save file called 'vox.mat' with 4 variables saved to it,
   and place this file in the main directory. These are the 4 variables:
 
      vox_train - The voxel size of the training image stacks and masks.
                  The code assumes that all tissue masks are stored at 
                  the same resolution.
 
      vox_avg   - The average atlas creation is actually performed at 
                  this resolution, and therefore this is the resolution 
                  at which registrations should be performed later on. If
                  there is any discrepancy between vox_stack or vox_mask
                  and this resolution, then the stack / masks will be
                  resampled to vox_avg during preprocessing.
 
      vox_reg   - The voxel size at which co-registration of the average
                  atlas onto the unknown sample will be performed. If
                  needed, the unknown 
 
      vox_unk   - The voxel size of the unknown samples. If
                  co-registration is performed at a different resolution,
                  the resultant atlas will be resampled to match the
                  original resolution of the unknown image.
 
3) Populate folder 1. This can be done directly by saving the raw image 
   stacks to folder 1. Alternately, the raw images can be saved to 
   folder 0, and WHAT.m can be used to resample these and save them to 
   folder 1.
 
4) Run ISA_1.M

Average Atlas Creation Step 2: ISA_2.M

ISA_2 script takes the average image created in ISA_1.M and applies it to
all tissue masks contained in the atlas, resulting in the creation of a
multi-tissue average atlas.
 
To create the average atlas successfully, the following steps must be
adhered to:
 
1) Follow the steps outlined in ISA_1.M and run ISA_1.M. This will create
   an average image from the training dataset and will store the
   deformations in a .MAT file.
 
2) Populate folder '2 - Training Masks'. If there are N training images
   and T tissues of interest, this folder should contain T folders, each
   of which contains N .MAT files containing the corresponding tissue
   masks.
 
3) Run ISA_2.M.

Average Atlas Optional Pre=Processing Step: ISA_0.M

 

1. set up HV_config file 
   example file is HV4_config.m in workspace\Results\HV4
   line 22 - line 26, set simulationDir and abaqus_command, make sure abaqus_command can run abaqus with user-subroutine
   All files for HV4 simulationDir can be found at workspace\HaoTemp\HV4_optimization.zip for reference.
  
2. In the segmentation folder, make sure you have a file named 
   optimization_config.m, in general no need to change this one. One example file can be found in workspace\Results\HV4
   
3. Various files to be prepared before run the simualtion
   * BiVentricleVolume.mat, see workspace\Results\HV4. Always set RV_endo_correction to be 0. The other values shall from the image derived geometry. We will need to discuss what shall be used for early_diastole and end_diastole volumes. The end-diastolic volume shall be measured from images since we dont reconstruct the geometry
   * BiVenMeshRVLV_SO.inp see workspace\Results\HV4\early_diastole\abaqusInput. You can load it into Abaqus, there are several sets to be defined, mainly the LV_ENDO, RV_ENDO, in order to apply pressure
   * fibTotal.inp, see workspace\Results\HV4\early_diastole\abaqusInput, its format is fibre_dir, sheet_dir 
   * LV_RV_assignment.mat, see workspace\Results\HV4\early_diastole\abaqusInput. node_assign: 1 for LV and 2 for RV. Similarly for elem_assign. It is used to calculate the cavity volume
   
3. TestAbaqusRun_PreOptimization  (NEED TO UPDATE PRESSURES)

4. Optimization_Ca_Cb

5. Optimization_refine_Ca_Cb

6. Optimization_af_bf

7. Optimization_a_afs_Ca_RV

There are some warnings from Matlab saying a file can not be opened. I have not find out where the error is, but seems not affecting the optimization procedure.
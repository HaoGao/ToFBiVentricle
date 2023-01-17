Step 1: Prepare all dicom images
        example: C:\Users\hg67u\Dropbox\work\BiVentricularReconstruction\MRIData

Step 2: Prepare the config file 
        HV25_Config.m as in C:\Users\hg67u\Dropbox\work\BiVentricularReconstruction\Results\HV25
	   
step 3: LVWM_DicomSampleSelection
       this will pop up a dialog to load the config file
	   
step 4: LVWM_SASegManual
        choose the right phase and segment the boundary from base to apex
        if it is the first time, then set fresh segmentation to be 1
        endo and epi separetely, but LV and RV are not separated for now
		
Step 5: LVWM_LASegManual
        choose the right phase and segment the long axis boundary 
        similar as the SA segmentation

Step 6: SAAdjustmentBYLASeries
		Adjust short-axis boundary according to long-axis boundary 
		
Step 7: optional to segment valve tract: LVWM_LVOTSegManual
        example HV5

Step 8: LVWM_DicomShowTogether
        output all boundaries in text file for generating geometry
        need to resegment the LV basal cavity to define the centre

Step 9: Geometry Generation (SolidWorks)
        export using step format

Step 10: Mesh generation (ICEM)
         example: \Results\HV1\early_diastole\solidworks\icem
         Choose size 3 for the mesh, in a range around 200k
	    10.1 first output file from ICEM: abaqus_earlyDiastole
        10.2 AbaqusInputAdjustedNodeNumber, after reading into Matlab code by choosing option 1
        10.3 Abaqus to define the right surfaces: 
    		Node Sets:    NODE_BASE, NODE_EPI, NODE_LV_ENDO, NODE_RV_ENDO (including both free wall and septum)
            Surface sets: SURF_BASE, SURF_EPI, SURF_LV_ENDO, SURF_RV_ENDO (including both free wall and septum), 
                              SURF_RV_SEPTUM_ENDO (only the rv septum surface)
        10.4 output an input file named as AbaqusInputSets.inp

	 
Step 11: meshConversion: A very important step 
         5 steps will run, each for different purpose
         vtk mesh for deformetric
	 libmesh mesh for fibre generation

Step 11: Diastolic pressure estimation 
         LVWM_AnnuRing_vel
         LVWM_MVInflow

Step 12: Strain estimation 
         LVWM_StrainAnalysis_splineDeform : needs to make sure only run on the right slices 

Step 13: LVWM_StrainAnalysis_summarize
         need to prepare a file deformRes.dat, example in HV1\bSpline

Step 14: LVWM_AHADefinition
         basically it will segment one slice in middel and one in apex
		 checked: pass
		 
Step 15: LVDivision_17Segments_AHA_Apply
         This will require the mesh file in a proper format
		 
Step 16: Fibre generation
	 source cshrc (using the 2017 version IBAMR, euclid-10)
         /xlwork6/hgao/BiVentricleModelling
         Upload the two files: Libmesh.dat and .node 
         Update mesh.in
         Run fibreGen-opt
         note: the angles are not outputted yet, only fibre direction, which will be updated later 24/09/2021

Step 17: AbaqusMeshGeneration_HO_Standard
         copy a main input file (BiVenMainStandard_SO)
         copy the fortran file (myuanisohyper_inv_sommer.for)

Step 17(a) run standard Abaqus simulation, need to use the lenova laptop, trying with School Desktop, seems ok.

Step 17(b): Test the Abaqus file for optimization procedure
		TestAbaqusRun_preOptimization
		Before to do any computation, needs to set up 
		(1) setup optimization_config.m at Results\HV1
		(2) prepare BiVenMainStandard_SO.inp, which is a fixed template for diastolic filling currently
		(3) prepare the strain energy function fortran file myuanisohyper_inv_sommer.for
		(4) copy BiVenMeshRVLV_SO.inp and fibTotal.inp

Optimization procedure 
Step 18: TestAbaqusRun_preOptimization
         this will output the strain and volumes for the initial sets of paramters defined in optimization_config.m
		 There is a correction for RV cavity volume calculation due to some wall volume maybe mis-included, thus the difference from the approximaiton - solidworks reported value, will be used for correction. 

Step 19: Optimization_Ca_Cb
Seems there is no need to run Ca Cb sweeping, may just start from some value and ask Matlab to determine the best one.
a) Optimization_refine_Ca_Cb
b) Optimization_af_bf
c) Optimization_a_afs_Ca_RV


		
To run the conda
Python is in /home/staff2/hgao/anaconda3/bin
	./conda env list
	conda activate deformetrica
	cd /xlwork6/hgao/BiVentricleModelling/deformetrica/HVHao
	run.sh
	

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		

Step 1: Prepare all dicom images
        example: C:\Users\hg67u\Dropbox\work\BiVentricularReconstruction\MRIData

Step 2: Prepare the config file 
        HV25_Config.m as in C:\Users\hg67u\Dropbox\work\BiVentricularReconstruction\Results\HV25
	   
step 3: LVWM_DicomSampleSelection
       this will pop up a dialog to load the config file
	   
step 4: LVWM_SASegManual
        choose the right phase and segment the boundary from base to apex
		
Step 5: LVWM_LASegManual
        choose the right phase and segment the long axis boundary 

Step 6: SAAdjustmentBYLASeries
		Adjust short-axis boundary according to long-axis boundary 
		
Step 7: optional to segment valve tract: LVWM_LVOTSegManual
        example HV5

Step 8: LVWM_DicomShowTogether
        output all boundaries in text file for generating geometry

Step 9: Geometry Generation (SolidWorks)

Step 10: Mesh generation (ICEM)
         example: \Results\HV1\early_diastole\solidworks\icem		 

Step 11: meshConversion: A very important step 
         5 steps will run, each for different purpose
        vtk mesh for deformetric
		libmesh mesh for fibre generation

Step 11: Diastolic pressure estimation 
         LVWM_AnnuRing_vel
         LVWM_MVInflow

Step 12: Strain estimation 
         LVWM_StrainAnalysis_splineDeform : needs to make sure only run on the right slices. Note that the first slice should be the one for constructing ventricular geometry at early-diastole. 

Step 13: LVWM_StrainAnalysis_summarize
         need to prepare a file deformRes.dat, example in HV1\bSpline

Step 14: LVWM_AHADefinition
         basically it will segment one slice in middel and one in apex
		 checked: pass
		 
Step 15: LVDivision_17Segments_AHA_Apply
         This will require the mesh file in a proper format
		 
Step 16: Fibre generation
         /xlwork6/hgao/BiVentricleModelling

Step 17: Test the Abaqus file
		TestAbaqusRun_preOptimization
		Before to do any computation, needs to set up 
		(1) setup optimization_config.m at Results\HV1
		(2) prepare BiVenMainStandard_SO.inp, which is a fixed template for diastolic filling currently
		(3) prepare the strain energy function fortran file myuanisohyper_inv_sommer.for
		(4) copy BiVenMeshRVLV_SO.inp and fibTotal.inp

Optimization procedure 

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		

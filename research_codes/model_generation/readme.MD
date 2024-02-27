
* Step 1: set up the config file and prepare dicom images

* Step 2: LVWM_DicomSampleSelection

	load dicom images into the Matlab workspace and saved for later usage

* Step 3: LVWM_SASegManual

   manually segment short axial images
   - case 1: LV, RV, EPI, for the bi-ventrilce region
   - case 2: RV, PA, EPI, not used anymore
   - case 3: RV, RV EPI, PA, PA EPI, not used
   - case 4: PA, PA EPI, for the RV outflow tract

* Step 4: LVWM_LASegManual

* Step 5: SAAdjustmentByLASeries

* Step 6: LVWM_DicomShowTogether
          LVWM_DicomShowTogether_v11
		  output: data points can be used for geometry generation
		        SA_LV
				SA_RV
				SA_EPI

* Step 7: Geometry Generation (details)
         surface generation: LVENDO, RVENDO and EPI
		 solid geometry generation through workbench: HV 
		 mesh generation to Abaqus
		 Abaqus set definition

* FibreConstruction (details)
         MakeFibreDir_v8_explicit 
		 	 
* Passive Parameter Inference (details)
      in template files: 
      BiVenMainStandard_SO_EXP: 
	  BiVenMainStandard_SO_EXP - PASSIVE: only has the passive step, using explicit solver
	  BiVenMainStandard_SO_EXP - ACTIVE:

		 
  

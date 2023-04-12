
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

* Step 7: Geometry Generation

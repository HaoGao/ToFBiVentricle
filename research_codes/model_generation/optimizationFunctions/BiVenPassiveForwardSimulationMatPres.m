function ResSim = BiVenPassiveForwardSimulationMatPres(optimize_opt, A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                              Ca_RV, press)

workingDir = pwd();
abaqusSimulationDir = optimize_opt.abaqusSimulationDir;
subroutine_file = optimize_opt.subroutine_file;
cpu_number = optimize_opt.cpu_number;
abaqus_inputfile = optimize_opt.abaqus_inputfile;
pythonOriginalFilesDir = optimize_opt.pythonOriginalFilesDir;
pythonfilename = optimize_opt.pythonfilename;
lineNoForODBName = optimize_opt.lineNoForODBName;
lineNoForDisName = optimize_opt.lineNoForDisName;
abaqus_dis_out_filename = optimize_opt.abaqus_dis_out_filename;
LVendoNodes = optimize_opt.LVendoNodes;
LV_EPI_RVSeptumENDO = optimize_opt.LV_EPI_RVSeptumENDO;
BiVentricleVolume = optimize_opt.BiVentricleVolume;
                                          
%% regenerate updated abaqus input file
abaqusInputFileUpdate_MatModel(optimize_opt, A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                        Ca_RV, press);
[status, result] = RunAbaqusJobFromMatlab(abaqusSimulationDir, abaqus_inputfile, ...
                      cpu_number, subroutine_file);
SuccessB = AbaqusRunningSuccessOrNot(abaqusSimulationDir, abaqus_inputfile);

if (SuccessB == 1)
    %%% now extract the displacement fields
    DislacementExtractByPython(abaqusSimulationDir,abaqus_inputfile,...
                            pythonOriginalFilesDir,pythonfilename, ...
                            lineNoForODBName, lineNoForDisName,...
                            abaqus_dis_out_filename)

    %% volume calculation
    ResVol = biVenCavityCalculation(optimize_opt);
    %% correct the RV cavity volume, make sure RV_endo_correction is defined in the test run
    ResVol.RV_vol_update  = ResVol.RV_vol_update - BiVentricleVolume.RV_endo_correction;
    ResVol.RV_vol_ori     = ResVol.RV_vol_ori - BiVentricleVolume.RV_endo_correction;

    %% strain calculation
    using_fsn = 0; %% whether to project the strain along myofibre direction or not
    fiberStrain_crl = strainCalculationFromNodalDisplacements(optimize_opt, using_fsn);

    strain_segs= segRegionsStrainSummarization(optimize_opt, fiberStrain_crl, LVendoNodes);
    strain_segs_LVEpi= segRegionsStrainSummarization(optimize_opt, fiberStrain_crl, LV_EPI_RVSeptumENDO);
    
    %% output to the log file 
    cd(abaqusSimulationDir);
    fid_log = fopen(optimize_opt.logfile_name,'a');
    cd(workingDir);
    
    timestr =  datestr(clock());
    fprintf(fid_log, 'running a forward simulation\n');
    fprintf(fid_log, 'Step running on %s\n', timestr);
    fprintf(fid_log, 'abaqus running success: %d\n', SuccessB);
    fprintf(fid_log, 'parameters used for LV: %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f, \t%f mmHg\n', A, B, Af, Bf, An, Bn, Afs, Bfs, press);
    fprintf(fid_log, 'parameters used for RV: %f\n',Ca_RV);
%     fprintf(fid_log, 'LV volume: %f(target: %f)\n', strainComparison.LVVolumeAba, strainComparison.LVVolumeMRI);
%     %fprintf(fid_log, 'LV strain difference squared: %f\n', sum(feval_total.^2));
%     fprintf(fid_log, 'LV strain difference squared: %f using %d segments (ave: %f vs %f) \n', sum(feval_total.^2)^0.5, length(strainComparison.strainMRITotal), ...
%                                    mean(strainComparison.strainMRITotal),mean(strainComparison.strainAbaTotal) );
    fclose(fid_log);


    
%     figure; hold on;
%     plot(strain_segs(7:30), 'r-'); hold on;
%     plot(strain_segs_LVEpi(7:30), 'b-');
%     plot(  (strain_segs(7:30) + strain_segs_LVEpi(7:30))/2, 'b--');
    ResSim.ResVol = ResVol;
    ResSim.strain_segs = strain_segs;
    ResSim.strain_segs_LVEpi= strain_segs_LVEpi;
    ResSim.SuccessB = SuccessB;
    ResSim.using_fsn = using_fsn;
    ResSim.ParaMat.A = A;
    ResSim.ParaMat.B = B;
    ResSim.ParaMat.Af = Af;
    ResSim.ParaMat.Bf = Bf;
    ResSim.ParaMat.An = An;
    ResSim.ParaMat.Bn = Bn;
    ResSim.ParaMat.Afs = Afs;
    ResSim.ParaMat.Bfs = Bfs;
    ResSim.ParaMat.press = press;
    ResSim.ParaMat.Ca_RV = Ca_RV;
else
    ResSim.ResVol = [];
    ResSim.strain_segs = [];
    ResSim.strain_segs_LVEpi= [];
    ResSim.SuccessB = 0;
    ResSim.using_fsn = [];
    ResSim.ParaMat.A = A;
    ResSim.ParaMat.B = B;
    ResSim.ParaMat.Af = Af;
    ResSim.ParaMat.Bf = Bf;
    ResSim.ParaMat.An = An;
    ResSim.ParaMat.Bn = Bn;
    ResSim.ParaMat.Afs = Afs;
    ResSim.ParaMat.Bfs = Bfs;
    ResSim.ParaMat.press = press;
    ResSim.ParaMat.Ca_RV = Ca_RV;

end
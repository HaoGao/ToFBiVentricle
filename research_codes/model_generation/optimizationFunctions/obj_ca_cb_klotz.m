function feval_total_fmincon = obj_ca_cb_klotz(x)

workingDir = pwd();
global optimize_opt;
opt_log_filename = optimize_opt.logfile_name ;
abaqusDir = optimize_opt.abaqusSimulationDir;
weighting_varifold=optimize_opt.weighting_varifold;

% mpara = options.mpara;


%% that is for Ca and Cb
A = optimize_opt.mpara.A_opt*x(1);
B = optimize_opt.mpara.B_opt*x(2);
Af = optimize_opt.mpara.Af_opt*x(1);
Bf = optimize_opt.mpara.Bf_opt*x(2);
An = optimize_opt.mpara.An_opt*x(1);
Bn = optimize_opt.mpara.Bn_opt*x(2);
Afs = optimize_opt.mpara.Afs_opt*x(1);
Bfs = optimize_opt.mpara.Bfs_opt*x(2);
Ca_RV = optimize_opt.mpara.Ca_RV_opt*x(3);
press = optimize_opt.mpara.press0;

%%this is the constrain 
if A < 0.1
    A= 0.1;
end
if Af < 0.1
    Af = 0.1;
end
if An < 0.1
    An = 0.1;
end
if Afs < 0.1
    Afs = 0.1;
end

cd(abaqusDir);
fid_log = fopen(opt_log_filename ,'a');
cd(workingDir);
timestr =  datestr(clock());
fprintf(fid_log, '\n');
fprintf(fid_log, 'Step 2 optimization for Ca Cb based on Klotz curve on %s\n', timestr);
fprintf(fid_log, 'one iteration begins in obj_ca_cb_klotz: %s\n',timestr);


%% calling one forward simulation with provided 8 unknown parameters 

ResSim = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                                            A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                           Ca_RV, press);
fprintf(fid_log, 'finish one simulation with %f mmHg\n\n', press);
                                            
%% the LV volume: ResSim.ResVol.LV_vol_update
%% the RV volume: Resim.ResVol.RV_vol_update
%% from measurement: optimize_opt.BiVentricleVolume.LV_end_diastole
%% optimize_opt.BiVentricleVolume.RV_end_diastole
LVVolumeAba = ResSim.ResVol.LV_vol_update;
RVVolumeAba = ResSim.ResVol.RV_vol_update;
LVEDVMRI = optimize_opt.BiVentricleVolume.LV_end_diastole;
RVEDVMRI = optimize_opt.BiVentricleVolume.RV_end_diastole;
%% the objective function
feval_total = (LVVolumeAba-LVEDVMRI)/LVEDVMRI; %relative volume difference
feval_total_fmincon = sum(feval_total.^2);

%% need to add the right ventricle part 
feval_total_RV = (RVVolumeAba-RVEDVMRI)/RVEDVMRI;
feval_total_RV_fmincon = sum(feval_total_RV.^2);

%%%call the klotz function, the original code takes so long time to compute
%%%the klotz curve
%[SumError, press_sim, EDV_sim, EDV_sim_norm, EDP_norm_klotzPrediction]= klotz_error_function(mpara_t);
%feval_total_fmincon = sum(SumError.^2)+feval_total_fmincon;

%% now we will only try to calculate the LEVDP=20mmHg, not all
V0 = optimize_opt.BiVentricleVolume.LV_early_diastole;
Ppred = optimize_opt.LVEDP_High;
Vpred = klotz_predictive_volume(press, LVEDVMRI, V0, Ppred);
ResSim_high = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                                            A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                           Ca_RV, Ppred);
fprintf(fid_log, 'finish one simulation with %f mmHg with volume: %f (%f klotz) \n\n', ...
                    Ppred, ResSim_high.ResVol.LV_vol_update,Vpred);

V_high = ResSim_high.ResVol.LV_vol_update;
SumError = (Vpred - V_high)/Vpred;


%% here we can try to use the varifold_distance to define the distance 
abaqusInput = optimize_opt.abaqusInput;
abaqusInput_end_diastole = optimize_opt.abaqusInput_end_diastole;
cd(abaqusDir);
dis = load(optimize_opt.abaqus_dis_out_filename); 
cd(workingDir);
[vari_dis, vari_dis_LV_endo] = varifold_distance_abaqus(abaqusInput, abaqusInput_end_diastole, dis);
vari_dis_max = optimize_opt.vari_dis_max;
feval_vari_fmincon = 3*weighting_varifold*(vari_dis.dis/(vari_dis_max.dis))^2;


%V0 = strainComparison.LVVolumeOri;
%V_ed = strainComparison.LVVolumeAba;
%V_ed = options.LVEDVMRI; %using the end-diastolic volume for prodicting what the value should be for LVEDP30
%V_30 = ( (V_ed - V0)*EDV_high  + EDV_norm*V0)/EDV_norm;%%predicted using klotz curve
%V_30_klotz = LVVolumeAba_30;
%SumError = (V_30_klotz-V_30)/V_30; 

feval_total_fmincon = sum(SumError.^2)+feval_total_fmincon + feval_total_RV_fmincon + feval_vari_fmincon;

cd(abaqusDir);
fid_loss = fopen('loss_function.data', 'a');
cd(workingDir);
fprintf(fid_loss, '%f\t %f\t %f\t %f\n', feval_total_fmincon,feval_total_fmincon,feval_total_RV_fmincon, feval_vari_fmincon);                 

fprintf(fid_log, 'abaqus running success for step 2 klotz: %d\n', ResSim.SuccessB);
fprintf(fid_log, 'x updated:          %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n', x(1), x(2), x(1), x(2), x(1), x(2), x(1), x(2), x(3));
fprintf(fid_log, 'parameters updated: %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f, \t%f\n', A, B, Af, Bf, An, Bn, Afs, Bfs, Ca_RV);
fprintf(fid_log, 'LV volume: %f(target: %f)\n', LVVolumeAba, LVEDVMRI);
fprintf(fid_log, 'RV volume: %f(target: %f)\n', RVVolumeAba, RVEDVMRI);
fprintf(fid_log, 'Difference (total, LV volume, klotz, RV volume): %f\t %f\t %f\t %f\n', sqrt(feval_total_fmincon),  ...
                                       sqrt(sum(feval_total.^2)), sqrt(sum(SumError.^2)), sqrt(feval_total_RV_fmincon));
fprintf(fid_log, 'EDV_high_sim: %f,\t EDV_high_predicted using klotz curve: %f with pressure %f\n', V_high, Vpred, Ppred);

fprintf(fid_log, 'vari_dis: %f (max: %f)\n', vari_dis.dis, vari_dis_max.dis );
fprintf(fid_log, 'one iteration ends\n');
fprintf(fid_log, '\n');
fclose(fid_log);

%assert(SuccessB==1);
%assert(SuccessB_30==1);
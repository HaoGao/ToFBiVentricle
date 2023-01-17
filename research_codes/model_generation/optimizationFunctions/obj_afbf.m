function feval_total_fmincon = obj_afbf(x)

workingDir = pwd();
global optimize_opt;
opt_log_filename = optimize_opt.logfile_name ;
abaqusDir = optimize_opt.abaqusSimulationDir;

strainData = optimize_opt.strainData;
strainDataFlag = optimize_opt.strainDataFlag;

%%only Af, Bf, As, Bs will be optimized
A = optimize_opt.mpara.A_opt;
B = optimize_opt.mpara.B_opt;
Af = optimize_opt.mpara.Af_opt*x(1);
Bf = optimize_opt.mpara.Bf_opt*x(2);
An = optimize_opt.mpara.An_opt;
Bn = optimize_opt.mpara.Bn_opt;
Afs = optimize_opt.mpara.Afs_opt;
Bfs = optimize_opt.mpara.Bfs_opt;
Ca_RV = optimize_opt.mpara.Ca_RV_opt;
press = optimize_opt.mpara.press0;

%%this is the constraint, need to think about what to do with this
if An >= Af/2
    An = Af/2;
end

if Bn>=Bf
    Bn = Bf;
end

cd(abaqusDir);
fid_log = fopen(opt_log_filename ,'a');
cd(workingDir);
timestr =  datestr(clock());
fprintf(fid_log, '\n');
fprintf(fid_log, 'Step 3 optimization for af and bf %s\n', timestr);
fprintf(fid_log, 'one iteration begins in LVPassive_model_optimization_fmincon_afbf: %s\n',timestr);

ResSim = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                                            A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                           Ca_RV, press);
fprintf(fid_log, 'finish one simulation with %f mmHg\n\n', press);

LVVolumeAba = ResSim.ResVol.LV_vol_update;
RVVolumeAba = ResSim.ResVol.RV_vol_update;
LVEDVMRI = optimize_opt.BiVentricleVolume.LV_end_diastole;
RVEDVMRI = optimize_opt.BiVentricleVolume.RV_end_diastole;


%% the objectie funciton for cicumferential strain 
NStrains = length(strainData);
diff_strain_Endo = [];
strainDataMRI = [];
strainDataAba = [];
for i = 1 : length(NStrains)
    if ~isnan(ResSim.strain_segs(i)) && strainDataFlag(i)>0.1
            diff_strain_Endo = [diff_strain_Endo ResSim.strain_segs(i) - strainData(i)];
            strainDataMRI = [strainDataMRI strainData(i)];
            strainDataAba = [strainDataAba ResSim.strain_segs(i)];
    end
end
 
%% the objective function for volume
feval_total = [diff_strain_Endo, (LVVolumeAba-LVEDVMRI)/LVEDVMRI]; %relative volume difference
feval_total_fmincon = sum(feval_total.^2); 

SuccessB = ResSim.SuccessB;
fprintf(fid_log, 'abaqus running success for updating af and bf: %d\n', SuccessB);
fprintf(fid_log, 'x updated:          %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n', 1.0, 1.0, x(1), x(2), 1.0, 1.0, 1.0, 1.0);
fprintf(fid_log, 'parameters updated: %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n', A, B, Af, Bf, An, Bn, Afs, Bfs);
fprintf(fid_log, 'LV volume: %f(target: %f)\n', LVVolumeAba, LVEDVMRI);
fprintf(fid_log, 'RV volume: %f(target: %f)\n', RVVolumeAba, RVEDVMRI);
fprintf(fid_log, 'strain: %f (target: %f)\n', mean(strainDataAba),   ...
                                              mean(strainDataMRI) );
fprintf(fid_log, 'Difference (total): %f\n', feval_total_fmincon);

fprintf(fid_log, 'iteration ends\n');
fprintf(fid_log, '\n');
fclose(fid_log);












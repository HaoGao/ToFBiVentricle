function feval_total_fmincon = obj_a_b_Ca_RV(x)

workingDir = pwd();
global optimize_opt;
opt_log_filename = optimize_opt.logfile_name ;
abaqusDir = optimize_opt.abaqusSimulationDir;

%%only Af, Bf, As, Bs will be optimized
A = optimize_opt.mpara.A_opt*x(1);
B = optimize_opt.mpara.B_opt*x(2);
Af = optimize_opt.mpara.Af_opt;
Bf = optimize_opt.mpara.Bf_opt;
An = optimize_opt.mpara.An_opt;
Bn = optimize_opt.mpara.Bn_opt;
Afs = optimize_opt.mpara.Afs_opt*x(1);
Bfs = optimize_opt.mpara.Bfs_opt*x(2);
Ca_RV = optimize_opt.mpara.Ca_RV_opt*x(3);
press = optimize_opt.mpara.press0;

%% constrain
if A < 0.02
    A = 0.02;
end
if Afs < 0.02
    Afs = 0.02;
end

cd(abaqusDir);
fid_log = fopen(opt_log_filename ,'a');
cd(workingDir);
timestr =  datestr(clock());
fprintf(fid_log, '\n');
fprintf(fid_log, 'Step 4 optimization for a and b for a final tune %s\n', timestr);
fprintf(fid_log, 'one iteration begins at : %s\n',timestr);
fclose(fid_log);

ResSim = BiVenPassiveForwardSimulationMatPres(optimize_opt, ...
                                            A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                           Ca_RV, press);
cd(abaqusDir);
fid_log = fopen(opt_log_filename ,'a');
cd(workingDir);                                       
fprintf(fid_log, 'finish one simulation with %f mmHg\n\n', press);
fclose(fid_log);

LVVolumeAba = ResSim.ResVol.LV_vol_update;
RVVolumeAba = ResSim.ResVol.RV_vol_update;
LVEDVMRI = optimize_opt.BiVentricleVolume.LV_end_diastole;
RVEDVMRI = optimize_opt.BiVentricleVolume.RV_end_diastole;

feval_total = [(LVVolumeAba-LVEDVMRI)/LVEDVMRI (RVVolumeAba-RVEDVMRI)/RVEDVMRI]; %relative volume difference
feval_total_fmincon = sum(feval_total.^2)*1000.0; 


SuccessB = ResSim.SuccessB;
cd(abaqusDir);
fid_log = fopen(opt_log_filename ,'a');
cd(workingDir);
fprintf(fid_log, 'abaqus running success for updating a, b, afs, and bfs: %d\n', SuccessB);
fprintf(fid_log, 'x updated:          %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n', x(1), x(2), 1, 1, 1.0, x(1), x(2), x(3));
fprintf(fid_log, 'parameters updated: %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f, \t%f\n', A, B, Af, Bf, An, Bn, Afs, Bfs, Ca_RV);
fprintf(fid_log, 'LV volume: %f(target: %f)\n', LVVolumeAba, LVEDVMRI);
fprintf(fid_log, 'RV volume: %f(target: %f)\n', RVVolumeAba, RVEDVMRI);
fprintf(fid_log, 'Relative Difference (total): %f\n', feval_total_fmincon);

fprintf(fid_log, 'iteration ends\n');
fprintf(fid_log, '\n');
fclose(fid_log);



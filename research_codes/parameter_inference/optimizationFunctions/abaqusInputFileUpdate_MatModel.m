function abaqusInputFileUpdate_MatModel(optimize_opt, A, B, Af, Bf, An, Bn, Afs, Bfs, ...
                                        Ca_RV, press)
 %%here we will assume right ventricular is propotial to LV, they should be
 %%more less the same
 
 abaqus_dir = optimize_opt.abaqusSimulationDir;
 abaqus_inputfile = optimize_opt.abaqus_inputfile;
 abaqus_inputfile_original = optimize_opt.abaqus_inputfile_original;
 materialParam_startLine_LV = optimize_opt.materialParam_startLine_LV;
 materialParam_startLine_RV = optimize_opt.materialParam_startLine_RV;
 pressureParam_startLine = optimize_opt.pressureParam_startLine;
 
 workingDir = pwd();
 
 cd(abaqus_dir);

 %cmd_delete_abaqus_odbfile = sprintf('%s.odb', abaqus_inputfile);
 %if exist(cmd_delete_abaqus_odbfile, 'file')
 %   dos(cmd_delete_abaqus_odbfile);
 %end
 
 %abaqus_inputfile = sprintf('%s.inp', abaqus_inputfile);
 %abaqus_inputfile_original = sprintf('%s.inp', abaqus_inputfile_original);
 fid_abaqus_inputfile = fopen(abaqus_inputfile, 'w');
 fid_abaqus_original = fopen(abaqus_inputfile_original, 'r');
 
 cd(workingDir);
 
 lineIndex = 1;
 tline = fgetl(fid_abaqus_original);
 fprintf(fid_abaqus_inputfile,'%s\r\n',tline);
 
 
 while ~feof(fid_abaqus_original)
    lineIndex = lineIndex + 1;
    tline = fgetl(fid_abaqus_original);
    
    if lineIndex == materialParam_startLine_LV
        tline = sprintf('A=%f*UNITCHANGEP',A);
    end
    if lineIndex == materialParam_startLine_LV+1
        tline = sprintf('B=%f',B);
    end
    if lineIndex == materialParam_startLine_LV+2
        tline = sprintf('Af=%f*UNITCHANGEP',Af);
    end
    if lineIndex == materialParam_startLine_LV+3
        tline = sprintf('Bf=%f',Bf);
    end
    if lineIndex == materialParam_startLine_LV+4
        tline = sprintf('An=%f*UNITCHANGEP',An);
    end
    if lineIndex == materialParam_startLine_LV+5
        tline = sprintf('Bn=%f',Bn);
    end
    if lineIndex == materialParam_startLine_LV+6
        tline = sprintf('Afs=%f*UNITCHANGEP',Afs);
    end
    if lineIndex == materialParam_startLine_LV+7
        tline = sprintf('Bfs=%f',Bfs);
    end
    
    %% now for right ventricle 
    if lineIndex == materialParam_startLine_RV
        tline = sprintf('Ca_RV = %f',Ca_RV);
    end
    if lineIndex == materialParam_startLine_RV + 1
        tline = sprintf('aR=%f*Ca_RV*UNITCHANGEP', A);
    end
    if lineIndex == materialParam_startLine_RV + 2
        tline = sprintf('bR=%f', B);
    end
    if lineIndex == materialParam_startLine_RV + 3
        tline = sprintf('afR=%f*Ca_RV*UNITCHANGEP', Af);
    end
    if lineIndex == materialParam_startLine_RV + 4
        tline = sprintf('bfR=%f', Bf);
    end
    if lineIndex == materialParam_startLine_RV + 5
        tline = sprintf('anR=%f*Ca_RV*UNITCHANGEP', An);
    end
    if lineIndex == materialParam_startLine_RV + 6
        tline = sprintf('bnR=%f', Bn);
    end
    if lineIndex == materialParam_startLine_RV + 7
        tline = sprintf('afsR=%f*Ca_RV*UNITCHANGEP', Afs);
    end
    if lineIndex == materialParam_startLine_RV + 8
        tline = sprintf('bfsR=%f', Bfs);
    end
    
    %% updating pressure
    if lineIndex == pressureParam_startLine 
       tline = sprintf(' P_mmHg = %f', press); 
    end
    
    
    fprintf(fid_abaqus_inputfile,'%s\r\n',tline);
 end

 fclose(fid_abaqus_inputfile);
 fclose(fid_abaqus_original);
 
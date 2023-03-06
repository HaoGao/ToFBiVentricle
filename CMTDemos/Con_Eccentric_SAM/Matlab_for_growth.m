%%Initially prepared by Hao Gao, modified by DB
clear all; close all; clc;



%read fibre data
fs_di=load('fibresheet.txt');
fs_dir=fs_di(:,2:7); 
%read node data from update file
node=load('node.txt');
thick=load('Thickness_Element.txt');
%load('SDV_C1');
SDV_stretch=load('Stretch_normal.txt');
SDV_stress=load('Stress_Normal.txt');

if exist('CVOLV.txt', 'file')==2
  delete('CVOLV.txt');
end
if exist('PCALV.txt', 'file')==2
  delete('PCALV.txt');
end

Make_InitialFile;

stateSum=ones(size(fs_dir,1),4); %should be volume fraction and total
stateSum(:,1)=0.274;stateSum(:,2)=0.7;stateSum(:,3)=0.026; %should be volume fraction and total

stateNew=ones(size(fs_dir,1),4); %should be growth amount at this cycle
stateOld=ones(size(fs_dir,1),4); %should be growth amount at this cycle
gv=ones(size(fs_dir,1),4);

filecopy={'FibreInp' 'FibreD' 'FGPI_sum' 'FGPI_C' 'FGPI_G' 'FGPI_M' 'Fr_C' 'Fr_G' 'Fr_M'...
'Udata' 'Ndata' 'Stretch' 'Stress' 'GrowthRatio'};

cycle_BIG=8

    pressure=(8+cycle_BIG)*133.322/1000000;
    stenosis=10.0;
    error_v=110;
    m_stretch=mean(SDV_stretch(:,3));
    m_stress=mean(SDV_stress(:,2));
%%    
for cycle=1:1
    %cases={'RBM_updata'  'Coarse_fibresheetDir.inp'};
    ['cycle=' num2str(cycle)]
    
%%***********************************************************************************************************
%%***********************************************************************************************************

    
	%update input file
	'update input file'
        fidr = fopen('Normal.inp','rt');
        filename='RBM_updata';
        fidw = fopen([filename '.inp'],'wt');   
        
        line=0;
		lineupdate=0;lineupdate2=0;
        while feof(fidr) == 0  %from first line to the last one 
            line=line+1;
            s = fgetl(fidr);  
            if line>=11 & line <= 26020
				lineupdate=lineupdate+1;
                fprintf(fidw,'\t%i,\t%14.10f,\t%14.10f,\t%14.10f\n',node(lineupdate,1),node(lineupdate,2),node(lineupdate,3),node(lineupdate,4));
            elseif  line>=161621 & line<=294662
                lineupdate2=lineupdate2+1;
                fprintf(fidw,'%i,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n',lineupdate2,fs_dir(lineupdate2,1),fs_dir(lineupdate2,2),fs_dir(lineupdate2,3), ...
                    fs_dir(lineupdate2,4),fs_dir(lineupdate2,5),fs_dir(lineupdate2,6));
            elseif  line==295011
                fprintf(fidw, [num2str(stenosis) ',0.\n']);
            elseif  line==295140
                fprintf(fidw, ['RP-LA, 8, 8, ' num2str(pressure) '\n']);
            elseif  line==295143
                fprintf(fidw, ['RP-CAV, 8, 8, ' num2str(pressure) '\n']);
            elseif  line==295221
                fprintf(fidw, ['RP-LA, 8, 8, ' num2str(pressure) '\n']);
            elseif  line==295426
                fprintf(fidw, ['RP-CAV, 8, 8, ' num2str(pressure) '\n']);
            else
                fprintf(fidw,'%s\n',s);
            end
            
        end
        
        fclose(fidr);
        fclose(fidw);
        
%%***********************************************************************************************************
%%***********************************************************************************************************

    %run abaqus
	'run abaqus'
    abaqus_inputfile=filename;
    command = sprintf('/maths/DassaultSystemes/SIMULIA/Commands/abq2018 job=%s user=newtest_Residual cpus=64 interactive double=both ask_delete=OFF',abaqus_inputfile);  
    [status,result] = system(command,'-echo');
    
    system('killall pre');
    system('killall standard');
    system('killall explicit');
    system('killall mpirun');
     
   
    copyfile('RBM_updata.inp', [abaqus_inputfile num2str(cycle) '.inp']);
    copyfile('RBM_updata.odb', [abaqus_inputfile num2str(cycle) '.odb']);

%%***********************************************************************************************************
%%***********************************************************************************************************

	%run python code to read node date after computation    
	'run python code to read node date after computation'  
    system(['/maths/DassaultSystemes/SIMULIA/Commands/abq2018 ' 'script=readNode_outSDV.py']); %Windows system? and output node.txt	


    ComputeFFGFR;
    %ComputeFFGFR_NoResidual;
    for fl=1:length(filecopy)
        copyfile([filecopy{fl} '.txt'], [filecopy{fl} num2str(cycle) '.txt']);
    end
    
        
    
   filename_workspace = ['cycle' num2str(cycle) '.mat'];
   save(filename_workspace)
        
end







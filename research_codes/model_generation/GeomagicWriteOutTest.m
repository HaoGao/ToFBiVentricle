%% write out for geomatric
%% SA
figure
hold on
endo_lv_data=[];endo_rv_data=[];epi_c_data=[];npo=50;

for imIndex = 1 : size(DataSegSARotated, 2)
            
    endo_lv = DataSegSARotated(imIndex).endo_lvReal;
    endo_rv = DataSegSARotated(imIndex).endo_rvReal;
    epi_c = DataSegSARotated(imIndex).epi_cReal;
    endo_lv_int=[];endo_rv_int=[];epi_c_int=[];
    if ~isempty(endo_lv)
        
        points=[];
        points = fnplt(cscvn(endo_lv));
        if length(points)>npo
            dp=round(length(points)/npo);
        else
            dp=1;
        end
        in=0;
        for i=1:dp:length(points)
            in=in+1;
            endo_lv_int(:,in)=points(:,i);
        end

        points=[];
        points = fnplt(cscvn(endo_rv));
        if length(points)>npo
            dp=round(length(points)/npo);
        else
            dp=1;
        end
        in=0;
        for i=1:dp:length(points)
            in=in+1;
            endo_rv_int(:,in)=points(:,i);
        end
        
        points=[];
        points = fnplt(cscvn(epi_c));
        if length(points)>npo+30
            dp=round(length(points)/(npo+30));
        else
            dp=1;
        end
        in=0;
        for i=1:dp:length(points)
            in=in+1;
            epi_c_int(:,in)=points(:,i);
        end
        
        endo_lv_data=[endo_lv_data; endo_lv_int'];
        endo_rv_data=[endo_rv_data; endo_rv_int'];
        epi_c_data=[epi_c_data; epi_c_int'];
        
        plot3(endo_lv_int(1,:),endo_lv_int(2,:),endo_lv_int(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2)
        plot3(endo_rv_int(1,:),endo_rv_int(2,:),endo_rv_int(3,:),'LineStyle', '-', 'Color', ...
                            'y', 'LineWidth',2)
        plot3(epi_c_int(1,:),epi_c_int(2,:),epi_c_int(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)
    end
            
end
    
if Bwritten
                    %%%output the bc boundaries
   cd(phase_resultDir);
   filename = sprintf('SA_LV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(endo_lv_data,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',endo_lv_data(i,1),endo_lv_data(i,2),endo_lv_data(i,3));
   end
   fclose(fidSA);
   
   filename = sprintf('SA_RV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(endo_rv_data,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',endo_rv_data(i,1),endo_rv_data(i,2),endo_rv_data(i,3));
   end
   fclose(fidSA);
   
   filename = sprintf('SA_EPI.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(epi_c_data,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',epi_c_data(i,1),epi_c_data(i,2),epi_c_data(i,3));
   end
   fclose(fidSA);
   
end


%%%%% LA
figure
hold on
endo_lv_data=[];endo_rv_data=[];epi_c_data=[];npo=50;

for imIndex = 1 : size(DataSegLARotated, 2)
            
    endo_lv = DataSegLARotated(imIndex).endo_lvReal;
    endo_rv = DataSegLARotated(imIndex).endo_rvReal;
    epi_c = DataSegLARotated(imIndex).epi_cReal;
    endo_lv_int=[];endo_rv_int=[];epi_c_int=[];
    if ~isempty(endo_lv)
        
        points=[];
        points = fnplt(cscvn(endo_lv));
        if length(points)>npo
            dp=round(length(points)/npo);
        else
            dp=1;
        end
        in=0;
        for i=1:dp:length(points)
            in=in+1;
            endo_lv_int(:,in)=points(:,i);
        end

        points=[];
        points = fnplt(cscvn(endo_rv));
        if length(points)>npo
            dp=round(length(points)/npo);
        else
            dp=1;
        end
        in=0;
        for i=1:dp:length(points)
            in=in+1;
            endo_rv_int(:,in)=points(:,i);
        end
        
        points=[];
        points = fnplt(cscvn(epi_c));
        if length(points)>npo+30
            dp=round(length(points)/(npo+30));
        else
            dp=1;
        end
        in=0;
        for i=1:dp:length(points)
            in=in+1;
            epi_c_int(:,in)=points(:,i);
        end
        
        
        endo_lv_data=[endo_lv_data; endo_lv_int'];
        endo_rv_data=[endo_rv_data; endo_rv_int'];
        epi_c_data=[epi_c_data; epi_c_int'];
        
        plot3(endo_lv_int(1,:),endo_lv_int(2,:),endo_lv_int(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2)
        plot3(endo_rv_int(1,:),endo_rv_int(2,:),endo_rv_int(3,:),'LineStyle', '-', 'Color', ...
                            'y', 'LineWidth',2)
        plot3(epi_c_int(1,:),epi_c_int(2,:),epi_c_int(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)
    end
            
end
    
if Bwritten
                    %%%output the bc boundaries
   cd(phase_resultDir);
   filename = sprintf('LA_LV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(endo_lv_data,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',endo_lv_data(i,1),endo_lv_data(i,2),endo_lv_data(i,3));
   end
   fclose(fidSA);
   
   filename = sprintf('LA_RV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(endo_rv_data,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',endo_rv_data(i,1),endo_rv_data(i,2),endo_rv_data(i,3));
   end
   fclose(fidSA);
   
   filename = sprintf('LA_EPI.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(epi_c_data,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',epi_c_data(i,1),epi_c_data(i,2),epi_c_data(i,3));
   end
   fclose(fidSA);
   
end



function GeomagicInputForClosedSplineGeneration()

%% SA
figure
hold on
endo_lv_data=[];endo_rv_data=[];epi_c_data=[];
theta_new = linspace(-pi, pi-0.01, 50);

for imIndex = 1 : size(DataSegSARotated, 2)
            
    endo_lv = DataSegSARotated(imIndex).endo_lvReal;
    endo_rv = DataSegSARotated(imIndex).endo_rvReal;
    epi_c = DataSegSARotated(imIndex).epi_cReal;
    endo_lv_int=[];endo_rv_int=[];epi_c_int=[];
    if ~isempty(endo_lv)
        %%%%% lv
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(endo_lv));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
  
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(length(theta_new))=center(3);
        endo_lv_int=[x; y; z_new];
        endo_lv_data=[endo_lv_data; endo_lv_int];
        
                
        plot3(endo_lv_int(1,:),endo_lv_int(2,:),endo_lv_int(3,:),'LineStyle', '-', 'Color', ...
                            'b', 'LineWidth',2)

 
        %%%%% rv
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(endo_rv));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(length(theta_new))=center(3);
        endo_rv_int=[x; y; z_new];
        endo_rv_data=[endo_rv_data; endo_rv_int];
        
        plot3(endo_rv_int(1,:),endo_rv_int(2,:),endo_rv_int(3,:),'LineStyle', '-', 'Color', ...
                            'y', 'LineWidth',2)

        
        %%%%%% epi
        points=[];theta=[];rho=[];z=[];
        points = fnplt(cscvn(epi_c));
        %[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
        center=mean(points,2);
        for i=1:size(points,2)
            vec=points(1:2,i)-center(1:2);
            theta(i)=atan2(vec(2),vec(1));
            rho(i)=norm(vec);
        end
   
        [theta_uni,m1,n1] = unique(theta,'stable');
        rho_uni=rho(m1);
        
        r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');
        for i=1:length(theta_new)
            x(i)=r_new(i)*cos(theta_new(i))+center(1);
            y(i)=r_new(i)*sin(theta_new(i))+center(2);
        end
        
        z_new(length(theta_new))=center(3);
        epi_c_int=[x; y; z_new];
        epi_c_data=[epi_c_data; epi_c_int];
        
        plot3(epi_c_int(1,:),epi_c_int(2,:),epi_c_int(3,:),'LineStyle', '-', 'Color', ...
                            'r', 'LineWidth',2)

    end
    
    
        
end  


%%%%%%interpolation along z axis
%%%% LV
min_z=min(DataSegLARotated(2).endo_lvReal(3,:));
max_z=max(DataSegLARotated(2).endo_lvReal(3,:));
z_int= linspace(min_z, max_z, 20);
LV_new=[];
for i=1:50
    new_data=[];
    x=[];y=[];z=[];
    x=endo_lv_data(1:3:size(endo_lv_data),i);
    y=endo_lv_data(2:3:size(endo_lv_data),i);
    z=endo_lv_data(3:3:size(endo_lv_data),i);
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; zz];
    LV_new=[LV_new new_data];
end

%%%% RV
min_z=min(DataSegLARotated(2).endo_rvReal(3,:));
max_z=max(DataSegLARotated(2).endo_rvReal(3,:));
z_int= linspace(min_z, max_z, 20);
RV_new=[];
for i=1:50
    new_data=[];
    x=[];y=[];z=[];
    x=endo_rv_data(1:3:size(endo_rv_data),i);
    y=endo_rv_data(2:3:size(endo_rv_data),i);
    z=endo_rv_data(3:3:size(endo_rv_data),i);
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; zz];
    RV_new=[RV_new new_data];
end

%%%% EPI
min_z=min(DataSegLARotated(2).epi_cReal(3,:));
max_z=max(DataSegLARotated(2).epi_cReal(3,:));
z_int= linspace(min_z, max_z, 20);
EPI_new=[];
for i=1:50
    new_data=[];
    x=[];y=[];z=[];
    x=epi_c_data(1:3:size(epi_c_data),i);
    y=epi_c_data(2:3:size(epi_c_data),i);
    z=epi_c_data(3:3:size(epi_c_data),i);
    
    xx = spline(z,x,z_int);
    yy = spline(z,y,z_int);
    
    new_data=[xx; yy; zz];
    EPI_new=[EPI_new new_data];
end


    
    



 


        
if Bwritten
                    %%%output the bc boundaries
   cd(phase_resultDir);
   filename = sprintf('SA_LV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(LV_new,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',LV_new(i,1),LV_new(i,2),LV_new(i,3));
   end
   fclose(fidSA);
   
   filename = sprintf('SA_RV.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(RV_new,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',RV_new(i,1),RV_new(i,2),RV_new(i,3));
   end
   fclose(fidSA);
   
   filename = sprintf('SA_EPI.txt');
   fidSA = fopen(filename,'w');
   for i=1:size(EPI_new,1)
       fprintf(fidSA,'%14.10f \t%14.10f \t%14.10f\n',EPI_new(i,1),EPI_new(i,2),EPI_new(i,3));
   end
   fclose(fidSA);
   
end

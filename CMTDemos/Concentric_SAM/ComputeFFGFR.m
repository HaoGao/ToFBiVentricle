%clear all; close all;
'post processing'
%%***********************************************************************************************************
%read data after python
t1=readtable('FGPI_G.txt', 'HeaderLines', 1);
fgiv_G=table2array(t1);
t2=readtable('Fr_G.txt', 'HeaderLines', 1);
VFRS_G=table2array(t2);
t3=readtable('FGPI_M.txt', 'HeaderLines', 1);
fgiv_M=table2array(t3);
t4=readtable('Fr_M.txt', 'HeaderLines', 1);
VFRS_M=table2array(t4);
t5=readtable('FGPI_C.txt', 'HeaderLines', 1);
fgiv_C=table2array(t5);
t6=readtable('Fr_C.txt', 'HeaderLines', 1);
VFRS_C=table2array(t6);

%VFRS=load('');
%fs_di=load('fibresheet.txt');
    
node_d=load('Ndata.txt');
element_d=load('Edata.txt');
node_dxdydz_d=load('Udata.txt');
SDV1=load('Stretch.txt');
SDV2=load('Stress.txt');
    
% re order
idx=[];
[~,idx] = sort(node_d(:,1)); % sort just the first column
node_d = node_d(idx,:);

idx=[];
[~,idx] = sort(node_dxdydz_d(:,1)); % sort just the first column
node_dxdydz_d = node_dxdydz_d(idx,:);

idx=[];
[~,idx] = sort(SDV1(:,1)); % sort just the first column
SDV1 = SDV1(idx,:);

idx=[];
[~,idx] = sort(SDV2(:,1)); % sort just the first column
SDV2 = SDV2(idx,:);
%%***********************************************************************************************************
%update node set
node=node_d+node_dxdydz_d;
node(:,1)=node_d(:,1);

%%***********************************************************************************************************
%parameters for growth ground matrix
p_lamd_crg=1.1; p_pre_crg=0.012;
pf_vmaxg=1.4; pf_taog=1.0; pf_gamag=2.0;
ps_vmaxg=2.0; ps_taog=3.2; ps_gamag=2.0;

%parameters for growth myofibre
p_lamd_crm=m_stretch; p_pre_crm=m_stress;
pf_vmaxm=1.4; pf_taom=0.4; pf_gamam=1.0;
ps_vmaxm=2.0; ps_taom=0.2; ps_gamam=1.0;

%parameters for growth collagen fibre
p_lamd_crc=1.1; p_pre_crc=0.012;
pf_vmaxc=1.4; pf_taoc=1.0; pf_gamac=2.0;
ps_vmaxc=2.0; ps_taoc=3.2; ps_gamac=2.0;

pf_dt=1;ps_dt=1;
%%***********************************************************************************************************
%Compute total F
%read the message in the input file
%fs_dir=fs_di(:,2:7); 
%stateNew=ones(size(fs_dir,1),3); %should be growth amount at this cycle
%1: ground matrix 2: myofibre 3: collagen fibre


for i=1:size(fs_dir,1)
     xyztet=[];
     for j=1:4
            xyztet(1,j)=node_d(element_d(i,j+1),1+1);
            xyztet(2,j)=node_d(element_d(i,j+1),2+1);
            xyztet(3,j)=node_d(element_d(i,j+1),3+1);
        
            dxdydz(1,j)=node_dxdydz_d(element_d(i,j+1),1+1);
            dxdydz(2,j)=node_dxdydz_d(element_d(i,j+1),2+1);
            dxdydz(3,j)=node_dxdydz_d(element_d(i,j+1),3+1);
      end
    
        [abc, Vcol]=IsoTet4ShapeFunDer(xyztet);
    
	%the deformation gradient tensor
        F=[];
        F=eye(3)+(dxdydz*abc)/6/Vcol; % defined in global cardisen coordinate
	% determinant of F
		F=F*(det(F))^(-1/3);
        
        f0=fs_dir(i,1:3)';
        s0=fs_dir(i,4:6)';
        n0=cross(f0,s0);
        n0=n0/norm(n0);
        
        f1=F*f0;
        s1=F*s0;
        
        f2=f1/norm(f1);
        s2=s1/norm(s1);
        
        n2=cross(f2,s2);
        n2=n2/norm(n2);
        s2=cross(n2,f2);
        s2=s2/norm(s2);
        
  %%***********************************************************************************************************
  % remodeling fibre direction
        fs_dir(i,1:3)=f2;fs_dir(i,4:6)=s2;
        
  %%***********************************************************************************************************
  % compute Fg,Fr1 at current coordinate
        Fg1_G=[fgiv_G(i,1) fgiv_G(i,4) fgiv_G(i,9); fgiv_G(i,7) fgiv_G(i,2) fgiv_G(i,5);fgiv_G(i,6) fgiv_G(i,8) fgiv_G(i,3)];
        Fg1_M=[fgiv_M(i,1) fgiv_M(i,4) fgiv_M(i,9); fgiv_M(i,7) fgiv_M(i,2) fgiv_M(i,5);fgiv_M(i,6) fgiv_M(i,8) fgiv_M(i,3)];
		Fg1_C=[fgiv_C(i,1) fgiv_C(i,4) fgiv_C(i,9); fgiv_C(i,7) fgiv_C(i,2) fgiv_C(i,5);fgiv_C(i,6) fgiv_C(i,8) fgiv_C(i,3)];
		R1=[f0 s0 n0]; R2=[f2 s2 n2];
        FL=R1'*F*R1; % rotate to local coordinate
        vr_G=[VFRS_G(i,1) VFRS_G(i,4) VFRS_G(i,9); VFRS_G(i,7) VFRS_G(i,2) VFRS_G(i,5); VFRS_G(i,6) VFRS_G(i,8) VFRS_G(i,3)];
		vr_M=[VFRS_M(i,1) VFRS_M(i,4) VFRS_M(i,9); VFRS_M(i,7) VFRS_M(i,2) VFRS_M(i,5); VFRS_M(i,6) VFRS_M(i,8) VFRS_M(i,3)];
		vr_C=[VFRS_C(i,1) VFRS_C(i,4) VFRS_C(i,9); VFRS_C(i,7) VFRS_C(i,2) VFRS_C(i,5); VFRS_C(i,6) VFRS_C(i,8) VFRS_C(i,3)];
        %vr=[1 0 0; 0 1 0; 0 0 1];
        Fr1_G=FL*vr_G*Fg1_G;
		Fr1_M=FL*vr_M*Fg1_M;
		Fr1_C=FL*vr_C*Fg1_C;
        
  %%***********************************************************************************************************      
  % compute growth along fibre by maximum stretch SDV1(:,2)

  %1: groumd matrix
        %p_lamd=SDV1(i,2);
        %pf_v=stateOld(i,1);
        %if p_lamd >= p_lamd_crg	
        %    if pf_v>=pf_vmaxg
        %        stateNew(i,1)=stateOld(i,1);
        %    else
        %        for k=1:50			
	%				p_zhi=(pf_vmaxg-pf_v)/(pf_vmaxg-1.0);
	%				p_k=1.0/pf_taog*p_zhi^pf_gamag;
	%				p_phi=p_lamd-p_lamd_crg;
	%			
	%				p_kphi=p_phi*pf_gamag*p_k/(pf_v-pf_vmaxg);
	%				%p_phik=-1.0d0*p_lamd/pf_v/pf_v*p_k;
	%			
	%				p_R=pf_v-stateOld(i,1)-pf_dt*p_k*p_phi;
	%				p_K=1.0d0-(p_kphi)*pf_dt;
	%			
	%				p_RK=p_R/p_K;
	%				pf_v=pf_v-p_RK;
    %
	%				if abs(p_RK) <= 1.0E-9
	%					stateNew(i,1)=real(pf_v);
	%					break;
    %                end 
    %                pf_v=real(pf_v);
    %                
    %            end 	
    %            stateNew(i,1)=real(pf_v);
    %
    %        end
    %    else
    %        stateNew(i,1)=stateOld(i,1);
    %    end
        %%% no growth in ground matrix
        stateNew(i,1)=1.0;
        


  %2: myofibre
  %      p_lamd=SDV1(i,3);
  %      pf_v=stateOld(i,2);
   %     if p_lamd >= p_lamd_crm	
    %        if pf_v>=pf_vmaxm
     %           stateNew(i,2)=stateOld(i,2);
      %      else
       %         for k=1:100			
		%			p_zhi=(pf_vmaxm-pf_v)/(pf_vmaxm-1.0);
		%			p_k=1.0/pf_taom*p_zhi^pf_gamam;
		%			p_phi=p_lamd-p_lamd_crm;
		%		
		%			p_kphi=p_phi*pf_gamam*p_k/(pf_v-pf_vmaxm);
		%			%p_phik=-1.0d0*p_lamd/pf_v/pf_v*p_k;
		%		
		%			p_R=pf_v-stateOld(i,2)-pf_dt*p_k*p_phi;
		%			p_K=1.0d0-(p_kphi)*pf_dt;
		%		
		%			p_RK=p_R/p_K;
		%			pf_v=pf_v-p_RK;
%
%					if abs(p_RK) <= 1.0E-9
%						stateNew(i,2)=real(pf_v);
%						break;
 %                   end 
%                    pf_v=real(pf_v);
 %                   
  %              end 	
  %              stateNew(i,2)=real(pf_v);
%
 %           end
  %      else
   %         stateNew(i,2)=stateOld(i,2);
    %    end
   %%% no growth in myofibre
        stateNew(i,2)=1.0;
       
        %3: collagen fibre
   %     p_lamd=SDV1(i,4);
   %     pf_v=stateOld(i,3);
   %     if p_lamd >= p_lamd_crc	
   %         if pf_v>=pf_vmaxc
   %             stateNew(i,3)=stateOld(i,3);
   %         else
   %             for k=1:50			
	%				p_zhi=(pf_vmaxc-pf_v)/(pf_vmaxc-1.0);
	%				p_k=1.0/pf_taoc*p_zhi^pf_gamac;
	%				p_phi=p_lamd-p_lamd_crc;
	%			
	%				p_kphi=p_phi*pf_gamac*p_k/(pf_v-pf_vmaxc);
	%				%p_phik=-1.0d0*p_lamd/pf_v/pf_v*p_k;
	%			
	%				p_R=pf_v-stateOld(i,3)-pf_dt*p_k*p_phi;
	%				p_K=1.0d0-(p_kphi)*pf_dt;
	%			
	%				p_RK=p_R/p_K;
	%				pf_v=pf_v-p_RK;
    %
	%				if abs(p_RK) <= 1.0E-9
	%					stateNew(i,3)=real(pf_v);
	%					break;
    %                end 
    %                pf_v=real(pf_v);
    %                
    %            end 	
    %            stateNew(i,3)=real(pf_v);
    %
    %        end
    %    else
    %        stateNew(i,3)=stateOld(i,3);
    %    end
    %%% no growth in collagen fibre
        stateNew(i,3)=1.0;
    %%% same growth in collagen fibre as myocycte
        %stateNew(i,3)=stateNew(i,2);
        
 % compute growth cross fibre by maximum stress SDV1(:,2)+     
  
	    
        p_pre=SDV2(i,2);
        ps_v=stateOld(i,4);
        if  p_pre>= p_pre_crm
            if ps_v>=ps_vmaxm
                stateNew(i,4)=ps_vmaxm;
            else
                for k=1:100			
					p_zhi=abs(ps_vmaxm-ps_v)/(ps_vmaxm-1.0);
					p_k=1.0/ps_taom*p_zhi^ps_gamam;
					p_phi=abs(p_pre-p_pre_crm);
				
					p_kphi=p_phi*ps_gamam*p_k/(ps_v-ps_vmaxm);
					%p_phik=-1.0d0*p_lamd/pf_v/pf_v*p_k;
				
					p_R=ps_v-stateOld(i,4)-ps_dt*p_k*p_phi;
					p_K=1.0d0-(p_kphi)*ps_dt;
				
					p_RK=p_R/p_K;
					ps_v=ps_v-p_RK;

					if abs(p_RK) <= 1.0E-9
						stateNew(i,4)=real(ps_v);
						break;
                    end 
                    ps_v=real(ps_v);
                    
                end 	
                stateNew(i,4)=real(ps_v);

            end
        else
            stateNew(i,4)=stateOld(i,4);
        end 

  %%***********************************************************************************************************
        
  %2: myofibre
  % compute Fr then rotation Fr from last coorditnate to next coordinate     
        Fr2_M=R2'*R1*Fr1_M*R1'*R2;
  %VR decomposition toward new residual tensor
        C=Fr2_M'*Fr2_M;
        [A,B]=eig(C);
        D=sqrtm(B);
        U=A*D*A';
        R=Fr2_M*inv(U);
        Vr2_M=Fr2_M*inv(R);
        VI=inv(Vr2_M);      
  %equivalent growth tensor
        FgM=[];
        vf_M=VI(:,1)/norm(VI(:,1));
        vs_M=VI(:,2)/norm(VI(:,2));
        vn_M=VI(:,3)/norm(VI(:,3));	
		% no need to be orthogonal
        vn_M=cross(vf_M,vs_M);
        vs_M=cross(vn_M,vf_M);
        vn_M=vn_M/norm(vn_M);
        vs_M=vs_M/norm(vs_M);
        
        
        gv(i,2)=stateNew(i,2)/stateOld(i,2);  gv(i,4)=stateNew(i,4)/stateOld(i,4);
        stateOld(i,2)=stateNew(i,2);  stateOld(i,4)=stateNew(i,4);
        gml=gv(i,4);
        %vf=[1 0 0]';vs=[0 1 0]';
		
        for ia=1:3
            for ib=1:3
                FgM(ia,ib)=gv(i,2)*vf_M(ia)*vf_M(ib)+gv(i,4)*vs_M(ia)*vs_M(ib)+vn_M(ia)*vn_M(ib);
            end
        end
        %FgM=FgM+eye(3);       
        FgMI=inv(FgM); 
        
        
        
  %%***********************************************************************************************************
  %1: ground matrix
  % compute Fr then rotation Fr from last coorditnate to next coordinate     
        Fr2_G=R2'*R1*Fr1_G*R1'*R2;
  %VR decomposition toward new residual tensor
        C=Fr2_G'*Fr2_G;
        [A,B]=eig(C);
        D=sqrtm(B);
        U=A*D*A';
        R=Fr2_G*inv(U);
        Vr2_G=Fr2_G*inv(R);
        VI=inv(Vr2_G);      
  %equivalent growth tensor
        FgG=[];
        vf_G=VI(:,1)/norm(VI(:,1));
        vs_G=VI(:,2)/norm(VI(:,2));
        vn_G=VI(:,3)/norm(VI(:,3));
        vn_G=cross(vf_G,vs_G);
        vs_G=cross(vn_G,vf_G);
        vn_G=vn_G/norm(vn_G);
        vs_G=vs_G/norm(vs_G);
        
		gv(i,1)=stateNew(i,1)/stateOld(i,1);
        stateOld(i,1)=stateNew(i,1);
        %vf=[1 0 0]';vs=[0 1 0]';
        for ia=1:3
            for ib=1:3
                FgG(ia,ib)=1/sqrt(gml)*vf_M(ia)*vf_M(ib)+gml*vs_M(ia)*vs_M(ib)+1/sqrt(gml)*vn_M(ia)*vn_M(ib);
            end
        end
        %FgG=FgG+eye(3);       
        FgGI=inv(FgG); 
        

        
  %3: collagen fibre
  % compute Fr then rotation Fr from last coorditnate to next coordinate     
        Fr2_C=R2'*R1*Fr1_C*R1'*R2;
  %VR decomposition toward new residual tensor
        C=Fr2_C'*Fr2_C;
        [A,B]=eig(C);
        D=sqrtm(B);
        U=A*D*A';
        R=Fr2_C*inv(U);
        Vr2_C=Fr2_C*inv(R);
        VI=inv(Vr2_C);      
  %equivalent growth tensor
        FgC=[];
        vf_C=VI(:,1)/norm(VI(:,1));
        vs_C=VI(:,2)/norm(VI(:,2));
        vn_C=VI(:,3)/norm(VI(:,3));	
        vn_C=cross(vf_C,vs_C);
        vs_C=cross(vn_C,vf_C);
        vn_C=vn_C/norm(vn_C);
        vs_C=vs_C/norm(vs_C);
        
        
        gv(i,3)=stateNew(i,3)/stateOld(i,3);
        stateOld(i,3)=stateNew(i,3);
        %vf=[1 0 0]';vs=[0 1 0]';
        %vf=[1 0 0]';vs=[0 1 0]';
        for ia=1:3
            for ib=1:3
                FgC(ia,ib)=vf_M(ia)*vf_M(ib)+gml*vs_M(ia)*vs_M(ib)+vn_M(ia)*vn_M(ib);
                %FgC(ia,ib)=vf_C(ia)*vf_C(ib)+vs_C(ia)*vs_C(ib)+vn_C(ia)*vn_C(ib);
            end
        end
        %FgC=FgC+eye(3);
        %lambda_g=gv(i,2);
        %%%%%%%%%%%%%%%%%%%%%%%%% same growth as myocyte
        %FgC(1,1)=gv(i,3);FgC(2,2)=1; FgC(3,3)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  elastic
        %FgC(1,1)=1;FgC(2,2)=1; FgC(3,3)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   plastic
        %FgC(1,1)=lambda_g;FgC(2,2)=1/sqrt(lambda_g); FgC(3,3)=1/sqrt(lambda_g);
        FgCI=inv(FgC); 
        
        FrwriteG(i,:,:)=Vr2_G;
        FrwriteM(i,:,:)=Vr2_M;
		FrwriteC(i,:,:)=Vr2_C;
        
        FgwriteG(i,:,:)=FgGI;
        FgwriteM(i,:,:)=FgMI;
        FgwriteC(i,:,:)=FgCI;
        
        Ffwrite(i,1,:)=vf_M;
        Ffwrite(i,2,:)=vf_C;
        Ffwrite(i,3,:)=vs_C;
		Ffwrite(i,4,:)=vn_C;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newmass=stateSum(i,1)*gv(i,1)+stateSum(i,2)*gv(i,2)*gv(i,4)+stateSum(i,3)*gv(i,3)*gv(i,4);
        stateSum(i,1)=stateSum(i,1)*gv(i,1)/newmass;
        stateSum(i,2)=stateSum(i,2)*gv(i,4)*gv(i,2)/newmass;
        stateSum(i,3)=stateSum(i,3)*gv(i,4)*gv(i,3)/newmass;
        stateSum(i,4)=newmass;
  end

    
growthvaluegv(cycle,:)=mean(gv);
growthvaluetotal(cycle,:)=mean(stateNew);      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    id1=0; id2=0; id3=0; id4=0; id5=0; id6=0; id7=0; id8=0; id9=0;id10=0;
    ic1= zeros(1,4); ic2= zeros(1,4); ic3= zeros(1,4); ic4= zeros(1,4); ic5= zeros(1,4); ic6=0; ic7=0; ic8=0; ic9=0;ic10=0;
    ir1= zeros(1,4); ir2= zeros(1,4); ir3= zeros(1,4); ir4= zeros(1,4); ir5= zeros(1,4); ir6=0; ir7=0; ir8=0; ir9=0;ir10=0;
    il1= zeros(1,4); il2= zeros(1,4); il3= zeros(1,4); il4= zeros(1,4); il5= zeros(1,4); il6=0; il7=0; il8=0; il9=0;il10=0;
    iG1= zeros(1,4); iG2= zeros(1,4); iG3= zeros(1,4); iG4= zeros(1,4); iG5= zeros(1,4); il6=0; il7=0; il8=0; il9=0;il10=0;
%%%%base
for i=1:size(fs_dir,1)
    
    tic=thick(i,2);
    if tic>0.0 & tic<=0.2
        id1=id1+1;
        ic1=ic1+stateNew(i,:);
        ir1=ir1+gv(i,:);
        il1=il1+SDV1(i,:);
        iG1=iG1+stateSum(i,:);
    elseif tic>0.2 & tic<=0.4
        id2=id2+1;
        ic2=ic2+stateNew(i,:);
        ir2=ir2+gv(i,:);
        il2=il2+SDV1(i,:);
        iG2=iG2+stateSum(i,:);
    elseif tic>0.4 & tic<=0.6
        id3=id3+1;
        ic3=ic3+stateNew(i,:);
        ir3=ir3+gv(i,:);
        il3=il3+SDV1(i,:);
        iG3=iG3+stateSum(i,:);
    elseif tic>0.6 & tic<=0.8
        id4=id4+1;
        ic4=ic4+stateNew(i,:);
        ir4=ir4+gv(i,:);
        il4=il4+SDV1(i,:);
        iG4=iG4+stateSum(i,:);
    else
        id5=id5+1;
        ic5=ic5+stateNew(i,:);
        ir5=ir5+gv(i,:);
        il5=il5+SDV1(i,:);
        iG5=iG5+stateSum(i,:);
    end
    
end

growlayerTO(cycle,1,:)=ic1/id1;
growlayerTO(cycle,2,:)=ic2/id2;
growlayerTO(cycle,3,:)=ic3/id3;
growlayerTO(cycle,4,:)=ic4/id4;
growlayerTO(cycle,5,:)=ic5/id5;

growlayerSI(cycle,1,:)=ir1/id1;
growlayerSI(cycle,2,:)=ir2/id2;
growlayerSI(cycle,3,:)=ir3/id3;
growlayerSI(cycle,4,:)=ir4/id4;
growlayerSI(cycle,5,:)=ir5/id5;

lamdlayerSI(cycle,1,:)=il1/id1;
lamdlayerSI(cycle,2,:)=il2/id2;
lamdlayerSI(cycle,3,:)=il3/id3;
lamdlayerSI(cycle,4,:)=il4/id4;
lamdlayerSI(cycle,5,:)=il5/id5;

GMlayerSI(cycle,1,:)=iG1/id1;
GMlayerSI(cycle,2,:)=iG2/id2;
GMlayerSI(cycle,3,:)=iG3/id3;
GMlayerSI(cycle,4,:)=iG4/id4;
GMlayerSI(cycle,5,:)=iG5/id5;


   % fid1 = fopen('fibresheet.txt','w');
   % for i = 1 : size(fs_dir,1)
   %     fprintf(fid1, '%i,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
   %         i, fs_dir(i,1), fs_dir(i,2), fs_dir(i,3), fs_dir(i,4), fs_dir(i,5),fs_dir(i,6));
   % end
   % fclose(fid1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fid2 = fopen('Fr_G.txt','w');
    fprintf(fid2, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid2, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            FrwriteG(i,1,1), FrwriteG(i,2,2), FrwriteG(i,3,3), FrwriteG(i,1,2), FrwriteG(i,2,3), FrwriteG(i,3,1),FrwriteG(i,2,1),FrwriteG(i,3,2),FrwriteG(i,1,3));
    end
    fclose(fid2);
	
	fid6 = fopen('Fr_M.txt','w');
    fprintf(fid6, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid6, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            FrwriteM(i,1,1), FrwriteM(i,2,2), FrwriteM(i,3,3), FrwriteM(i,1,2), FrwriteM(i,2,3), FrwriteM(i,3,1),FrwriteM(i,2,1),FrwriteM(i,3,2),FrwriteM(i,1,3));
    end
    fclose(fid6);
	
	fid7 = fopen('Fr_C.txt','w');
    fprintf(fid7, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid7, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            FrwriteC(i,1,1), FrwriteC(i,2,2), FrwriteC(i,3,3), FrwriteC(i,1,2), FrwriteC(i,2,3), FrwriteC(i,3,1),FrwriteC(i,2,1),FrwriteC(i,3,2),FrwriteC(i,1,3));
    end
    fclose(fid7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fid3 = fopen('FGPI_G.txt','w');
    fprintf(fid3, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid3, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            FgwriteG(i,1,1), FgwriteG(i,2,2), FgwriteG(i,3,3), FgwriteG(i,1,2), FgwriteG(i,2,3), FgwriteG(i,3,1),FgwriteG(i,2,1),FgwriteG(i,3,2),FgwriteG(i,1,3));
    end
    fclose(fid3);
    
        
    fid4 = fopen('FGPI_M.txt','w');
    fprintf(fid4, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid4, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            FgwriteM(i,1,1), FgwriteM(i,2,2), FgwriteM(i,3,3), FgwriteM(i,1,2), FgwriteM(i,2,3), FgwriteM(i,3,1),FgwriteM(i,2,1),FgwriteM(i,3,2),FgwriteM(i,1,3));
    end
    fclose(fid4);
    
    
    fid5 = fopen('FGPI_C.txt','w');
    fprintf(fid5, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid5, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            FgwriteC(i,1,1), FgwriteC(i,2,2), FgwriteC(i,3,3), FgwriteC(i,1,2), FgwriteC(i,2,3), FgwriteC(i,3,1),FgwriteC(i,2,1),FgwriteC(i,3,2),FgwriteC(i,1,3));
    end
    fclose(fid5);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fid6 = fopen('FGPI_sum.txt','w');
    fprintf(fid6, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid6, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            stateSum(i,1),stateSum(i,2),stateSum(i,3),stateSum(i,4));
    end
    fclose(fid6);
    
    
    fid7 = fopen('FibreD.txt','w');
    fprintf(fid7, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid7, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
            Ffwrite(i,1,1), Ffwrite(i,1,2), Ffwrite(i,1,3), Ffwrite(i,2,1), Ffwrite(i,2,2), Ffwrite(i,2,3),...
			Ffwrite(i,3,1),Ffwrite(i,3,2),Ffwrite(i,3,3),Ffwrite(i,4,1),Ffwrite(i,4,2),Ffwrite(i,4,3));
    end
    fclose(fid7);
    
    
    
    fid9 = fopen('FibreInp.txt','w');
    fprintf(fid9, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid9, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
			fs_dir(i,1),fs_dir(i,2),fs_dir(i,3),fs_dir(i,4),fs_dir(i,5),fs_dir(i,6));
    end
    fclose(fid9);
    
    fid10 = fopen('GrowthRatio.txt','w');
    fprintf(fid10, '%i\n',size(fs_dir,1));
    for i = 1 : size(fs_dir,1)
        fprintf(fid10, '%14.10f,\t%14.10f,\t%14.10f,\t%14.10f\n', ...
			gv(i,1),gv(i,2),gv(i,3),gv(i,4));
    end
    fclose(fid10);
    
    
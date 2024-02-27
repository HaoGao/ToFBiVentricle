	 !DEC$ ATTRIBUTES ALIAS:"vumat"::VUMAT  
      subroutine vumat(
C Read only -
     *  jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *  stepTime, totalTime, dt, cmname, coordMp, charLength,
     *  props, density, strainInc, relSpinInc,
     *  tempOld, stretchOld, defgradOld, fieldOld,
     *  stressOld, stateOld, enerInternOld, enerInelasOld,
     *  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     *  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension jblock(*), props(nprops), density(*),
     *  coordMp(*),
     *  charLength(*), strainInc(*),
     *  relSpinInc(*), tempOld(*),
     *  stretchOld(*), defgradOld(*),
     *  fieldOld(*), stressOld(*),
     *  stateOld(*), enerInternOld(*),
     *  enerInelasOld(*), tempNew(*),
     *  stretchNew(*), defgradNew(*), fieldNew(*),
     *  stressNew(*), stateNew(*),
     *  enerInternNew(*), enerInelasNew(*)
C
      character*80 cmname
	  
       call vumatXtrArg(
     *  jblock(1), ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *  stepTime, totalTime, dt, cmname, coordMp, charLength,
     *  props, density, strainInc, relSpinInc,
     *  tempOld, stretchOld, defgradOld, fieldOld,
     *  stressOld, stateOld, enerInternOld, enerInelasOld,
     *  tempNew, stretchNew, defgradNew, fieldNew, 	  
     *  stressNew, stateNew, enerInternNew, enerInelasNew,
     *  jblock(5) )
	   
       return 
       end
	  
      subroutine vumatXtrArg(
C Read only -
     *  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *  stepTime, totalTime, dt, cmname, coordMp, charLength,
     *  props, density, strainInc, relSpinInc,
     *  tempOld, stretchOld, defgradOld, fieldOld,
     *  stressOld, stateOld, enerInternOld, enerInelasOld,
     *  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     *  stressNew, stateNew, enerInternNew, enerInelasNew,
     *  Nelement)
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock),
     *  F(ndir,ndir),B(ndir,ndir),C(ndir,ndir),ftf(ndir,ndir),sts(ndir,ndir),fnn(ndir),
     *  pntn(ndir,ndir),ftstf(ndir,ndir),sig(ndir,ndir),ftntf(ndir,ndir),FI(ndir,ndir),sigm(ndir,ndir),
     *  Fb(ndir, ndir), pIden(ndir,ndir),signew(ndir, ndir),U(ndir,ndir),UI(ndir,ndir),Ro(ndir,ndir),
     *  pfour(ndir, ndir, ndir, ndir),f0(ndir),fn0(ndir),sigman(ndir, ndir),RT(ndir, ndir),fNnNn(ndir,ndir),
     *  ff(ndir),fs(ndir),fn(ndir),Me(ndir, ndir),FeI(ndir, ndir),Fe(ndir, ndir),Nelement(nblock),fff(ndir),fss(ndir),
     *	Fg(ndir, ndir),Fgp(ndir, ndir),FgGI(ndir, ndir),FgMI(ndir, ndir),FgCI(ndir, ndir)
C  
      character*80 cmname
	  
      parameter ( half = 0.5d0,zero = 0.0d0, one  = 1.d0, 
     *            two  = 2.d0, three = 3.d0,index_J = 3,
     *            asmall  = 2.d-16, pi=3.1415926d0,
     *            p_lamd_cr=1.01d0,p_pre_cr=0.2d0,
     *            pf_vmax=1.5d0, pf_tao=3.2d0,pf_gama=2.0d0,
     *            ps_vmax=3.0d0, ps_tao=3.2d0,ps_gama=2.0d0)	 
      real*8   FrG(250000,9),FrM(250000,9),FrC(250000,9),Fg_G(250000,10),Fg_M(250000,10),Fg_C(250000,10),
     *         Fg_sum(250000,4),Fd(250000,12),VFrsG(3,3),VFrsM(3,3),VFrsC(3,3)
	  integer*4 SDVelement,StartAnalysis,SDVnode
      common /FibreSheetblk/FrG,FrM,FrC,Fg_G,Fg_M,Fg_C,Fg_sum,Fd,SDVelement,SDVnode,StartAnalysis
	  
C
C Read material properties
C  
      P_a  = props(1)
      P_b  = props(2)
      P_am = props(3)
      P_bm = props(4)
	  P_ac = props(5)
      P_bc = props(6)
	  P_acsn = props(7)
      P_bcsn = props(8)
	  P_lamd = props(9)

C - Compressible case
      P_D  = props(10)
      dinv=1/P_D
	  Time=totalTime
C - active 	 
	  t0= props(11)
	  p_m= props(12)
	  b_len= props(13)
	  p_l0= props(14)
	  b_sha= props(15)
	  Ca0_m= props(16)
	  Ca0= props(17)
	  Tmax= props(18)
	  a_thr= props(19)
	  p_nf= props(20)
	  p_ns= props(21)
	  p_nn= props(22)
	  
C--Define initial growth values	  
C--stateOld(k,1) for growth tensor,stateOld(k,2) for max stretch
C--write(*,*)Time,Nelement(k)-14,stateOld(k,1),StartAnalysis

		if (Time .EQ. zero .AND. StartAnalysis .EQ. 0) then
			call vextern
		end if
		
		if (Time .GT. zero .AND. Time .LE. dt .AND. StartAnalysis .EQ. 0) then
			call vextern
		end if
		
		if (Time .EQ. zero) then
C--stateOld(k,1) for PASSIVE growth tensor,stateOld(k,2) for max stretch
C--stateOld(k,3) for ACTIVE growth tensor,stateOld(k,4) for max stress
			do k=1, nblock
				stateOld(k,1)=1.0d0
				stateOld(k,2)=1.0d0
				stateOld(k,3)=1.0d0
				stateOld(k,4)=0.0d0
			end do
		end if
			  	  
c-------------------------------------
c-------------------------------------	  
C
C Please note here, initial value is not zero
      do k = 1, nblock
c-        write(*,*)Time,Nelement(k),Nelement(k)
C--Keep previous values in stepTime
C--	  write(*,*) Time, Nelement(k)
		
C--Compute deformation gradient tensor F
		
			FI(1,1)=defgradNew(k,1)
			FI(2,2)=defgradNew(k,2)
			FI(3,3)=defgradNew(k,3)
			FI(1,2)=defgradNew(k,4)
			FI(2,3)=defgradNew(k,5)
			FI(3,1)=defgradNew(k,6)
			FI(2,1)=defgradNew(k,7)
			FI(3,2)=defgradNew(k,8)
			FI(1,3)=defgradNew(k,9)
			

c- modify growth tensor at here	
			
c- growth tensor#		

c- start new method at here  ground matrix
			if (Time .LE. one) then
				ptim=Time
			else
				ptim=1.0
			end if
			
				pv=(Fg_G(Nelement(k),1)-1.0)*ptim+1.0
				FgGI(1,1)=1.0/pv
				pv=(Fg_C(Nelement(k),1)-1.0)*ptim+1.0
				FgGI(2,2)=1.0/pv
				pv=(Fg_M(Nelement(k),1)-1.0)*ptim+1.0
				FgGI(3,3)=1.0/pv

				FgGI(1,2)=zero
				FgGI(2,3)=zero
				FgGI(3,1)=zero
				FgGI(2,1)=zero
				FgGI(3,2)=zero
				FgGI(1,3)=zero
c-			write(*,*)  Time, p_gtime,FgGI(1,1),FgGI(2,2),FgGI(3,3)

c- FOR STRESS ROTATION

			U(1,1)=stretchNew(k,1)
			U(2,2)=stretchNew(k,2)
			U(3,3)=stretchNew(k,3)
			U(1,2)=stretchNew(k,4)
			U(2,3)=stretchNew(k,5)
			U(3,1)=stretchNew(k,6)
			U(2,1)=stretchNew(k,4)
			U(3,2)=stretchNew(k,5)
			U(1,3)=stretchNew(k,6)
			
		DETU = U(1,1)*U(2,2)*U(3,3) 
     *      - U(1,1)*U(2,3)*U(3,2)  
     *      - U(1,2)*U(2,1)*U(3,3)  
     *      + U(1,2)*U(2,3)*U(3,1)  
     *      + U(1,3)*U(2,1)*U(3,2)  
     *      - U(1,3)*U(2,2)*U(3,1)		
		
		UI(1,1)=1.0D0/DETU*(U(2,2)*U(3,3)-U(2,3)*U(3,2))
		UI(1,2)=1.0D0/DETU*(U(1,3)*U(3,2)-U(1,2)*U(3,3))
		UI(1,3)=1.0D0/DETU*(U(1,2)*U(2,3)-U(1,3)*U(2,2))
		UI(2,1)=1.0D0/DETU*(U(2,3)*U(3,1)-U(2,1)*U(3,3))
		UI(2,2)=1.0D0/DETU*(U(1,1)*U(3,3)-U(1,3)*U(3,1))
		UI(2,3)=1.0D0/DETU*(U(1,3)*U(2,1)-U(1,1)*U(2,3))
		UI(3,1)=1.0D0/DETU*(U(2,1)*U(3,2)-U(2,2)*U(3,1))
		UI(3,2)=1.0D0/DETU*(U(1,2)*U(3,1)-U(1,1)*U(3,2))
		UI(3,3)=1.0D0/DETU*(U(1,1)*U(2,2)-U(1,2)*U(2,1))


c  ***********************************************************
c  ***********************************************************
C- For ground matrix

		do i=1, ndir
			do j=1, ndir
				tmp=0.0d0
				do ia=1, ndir
					tmp=tmp+FI(i,ia)*FgGI(ia,j)
				end do
				Fe(i,j)=tmp
			end do
		end do 
		
			
C--Compute the determinant 	 
		DET = Fe(1,1)*Fe(2,2)*Fe(3,3) 
     *      - Fe(1,1)*Fe(2,3)*Fe(3,2)  
     *      - Fe(1,2)*Fe(2,1)*Fe(3,3)  
     *      + Fe(1,2)*Fe(2,3)*Fe(3,1)  
     *      + Fe(1,3)*Fe(2,1)*Fe(3,2)  
     *      - Fe(1,3)*Fe(2,2)*Fe(3,1)

C- F is decomposed
        do i=1, ndir
			do j=1, ndir			
				Fb(i,j)=Fe(i,j)*DET**(-1.0/3.0)
			end do
		end do

C-Compute the C, B, f tensor f, s tensor s			
		do i=1, ndir
		  do j=1, ndir
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+Fb(i,ij)*Fb(j,ij)
			end do
			B(i,j)=tmp
					
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+Fb(ij,i)*Fb(ij,j)
			end do
			C(i,j)=tmp
			
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+FI(i,ij)*UI(ij,j)
			end do
			Ro(i,j)=tmp
			
			tmp=0.0d0
			pIden(i,j)=tmp
            ftf(i,j)=Fb(i,1)*Fb(j,1)
			sts(i,j)=Fb(i,2)*Fb(j,2)
			pntn(i,j)=Fb(i,3)*Fb(j,3)
		  end do
		end do 
C--Identity tensor			
			pIden(1,1)=1.0d0
			pIden(2,2)=1.0d0
			pIden(3,3)=1.0d0
C--Invariant 
		    p_I1=C(1,1)+C(2,2)+C(3,3)
			p_I4f=C(1,1)
			p_I4s=C(2,2)
			p_I4n=C(3,3)
			p_I4fm=max(p_I4f,1.0d0)
			p_I4fc=max(p_I4f/P_lamd/P_lamd,1.0d0)
			p_I4sc=max(p_I4s/P_lamd/P_lamd,1.0d0)
			p_I4nc=max(p_I4n/P_lamd/P_lamd,1.0d0)

							
c ****************************

		do i=1, ndir
			do j=1, ndir
				sig(i,j)=2.0/DET*P_a*(B(i,j)-p_I1/3.0d0*pIden(i,j))
     *        		+2.0d0/DET*P_am*(p_I4fm-1.0d0)*exp(P_bm*(p_I4fm-1.0d0)**2)*
     *					(ftf(i,j)-p_I4fm/3.0d0*pIden(i,j))
     *        		+2.0d0/DET*P_ac*(p_I4fc-1.0d0)*exp(P_bc*(p_I4fc-1.0d0)**2)*
     *					(ftf(i,j)-p_I4fc/3.0d0*pIden(i,j))
     *        		+2.0d0/DET*P_acsn*(p_I4sc-1.0d0)*exp(P_bcsn*(p_I4sc-1.0d0)**2)*
     *					(sts(i,j)-p_I4sc/3.0d0*pIden(i,j))
     *        		+2.0d0/DET*P_acsn*(p_I4nc-1.0d0)*exp(P_bcsn*(p_I4nc-1.0d0)**2)*
     *					(pntn(i,j)-p_I4nc/3.0d0*pIden(i,j))
     *				+(DET**2-1.0d0)/DET*dinv*pIden(i,j)
			end do
		end do				
  		
			
c  ***********************************************************
      Time=totalTime
      t=0.0d0
		 p_lr=0.00185d0
         p_lff=p_lr*sqrt(p_I4f)
		 stateNew(k,3)=sqrt(p_I4f)
		 
      if (tempOld(k).GE.a_thr .AND. stateOld(k,4).LT.asmall) then
	     stateNew(k,4)=Time
c		 stateNew(k,3)=p_m*p_lff+b_len
      elseif (stateOld(k,4).GE.asmall) then
	     t=Time-stateOld(k,4)
	     stateNew(k,4)=stateOld(k,4)
c		 stateNew(k,3)=stateOld(k,3)
      else
	     t=0.0d0
		 stateNew(k,4)=stateOld(k,4)
      end if 
	  
c 	  - Ta

         t_r=p_m*p_lff+b_len
		 if (t_r .LT. zero) then
			t_r=0.001d0
		 end if
		 
c-		 if (t_r .GT. t0) then
c-			t_r=t_r/t0*0.001d0+t0
c-		 end if
		 
         t1=t_r+t0
         w_t=zero
		 
         if(t .GT. 0.0d0 .and. t .LE. t0)then 
            w_t=pi*t/t0
         end if
         if(t .GT. t0 .and. t .LE. t1)then
            w_t=pi*(t-t0+t_r)/t_r
         end if	 
         if(t .GT. t1)then
            w_t=0.0d0
			stateNew(k,4)=zero
         end if
c-         if(t .GT. 2.0*t0)then
c-            w_t=0.0d0
c-			stateNew(k,4)=zero
c-        end if		 	

         if(t .GT. zero .AND. p_lff .GT. p_l0)then
            term1=exp(b_sha*(p_lff-p_l0))
            Eca=Ca0_m/sqrt(term1-1.0d0)
			term3=Ca0**2+Eca**2
            term2=Ca0**2/term3

            Ta=half*Tmax*term2*(1.0d0-cos(w_t))

            stateNew(k,2)=Ta

c-------------------------------------

			do i=1, ndir
				do j=1, ndir
					sig(i,j)=sig(i,j)+p_nf*Ta/p_I4f*ftf(i,j)+p_ns*Ta/p_I4s*sts(i,j)+p_nn*Ta/p_I4n*pntn(i,j)
				end do
			end do
			 
         end if
c-------------------------------------
c-------------------------------------		
c  ***********************************************************
c  ***********************************************************		 
c- Do finial rotation back to initial configuration		 
		  do i=1, ndir
			do j=1, ndir
				tmp=0.0d0
				do ia=1, ndir
					do ib=1, ndir
						tmp=tmp+Ro(ia,i)*sig(ia,ib)*Ro(ib,j)
					end do
				end do
				signew(i,j)=tmp
			end do
		  end do
			


c  ***********************************************************
c-Update stress		
			stressNew(k,1) = signew(1,1)
			stressNew(k,2) = signew(2,2)
			stressNew(k,3) = signew(3,3)
			stressNew(k,4) = signew(1,2)
			stressNew(k,5) = signew(2,3)
			stressNew(k,6) = signew(3,1)


c                   write(*,*) p_I4fn,fNnNn(1,1)		
C
C
C Update the specific internal energy -
C
			
			stressPower =half * (
     * 		( stressOld(k,1) + stressNew(k,1) ) * strainInc(k,1) +
     * 		( stressOld(k,2) + stressNew(k,2) ) * strainInc(k,2) +
     * 		( stressOld(k,3) + stressNew(k,3) ) * strainInc(k,3) ) +
     * 		( stressOld(k,4) + stressNew(k,4) ) * strainInc(k,4) +
     * 		( stressOld(k,5) + stressNew(k,5) ) * strainInc(k,5) +
     * 		( stressOld(k,6) + stressNew(k,6) ) * strainInc(k,6)
	 
			enerInternNew(k) = enerInternOld(k) + stressPower / density(k)
C
C Update the dissipated inelastic specific energy -
C
C   
	 
	        plasticWorkInc = 0.0d0
			enerInelasNew(k) = enerInelasOld(k)
     * 		+ plasticWorkInc / density(k)
	 
      end do
	  
		
      return
      end
 


c read fibre direction and used for tranform stress into local material coordinate 	  
      subroutine vextern
      
      include 'vaba_param.inc'   
	   
      integer*4 i,j, elemNTotal
      real*8    FrG(250000,9),FrM(250000,9),FrC(250000,9),Fg_G(250000,10),Fg_M(250000,10),Fg_C(250000,10),
     *          Fg_sum(250000,4),Fd(250000,12)
	  integer*4 SDVelement,StartAnalysis,SDVnode
      common /FibreSheetblk/FrG,FrM,FrC,Fg_G,Fg_M,Fg_C,Fg_sum,Fd,SDVelement,SDVnode,StartAnalysis
	  
	    StartAnalysis=StartAnalysis+1
	  	
		  
C- Displacement DATA, out order

			Fg_G(:,:)=0.0D0
	        open(403,FILE='D:/DBGuan/LianTian/RatHeart/AbaqusFile/Growth/FGPI_Gff_3.txt',STATUS='OLD')
            read(403,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(403, *) Fg_G(i,1)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(403)
			
			Fg_C(:,:)=0.0D0
	        open(404,FILE='D:/DBGuan/LianTian/RatHeart/AbaqusFile/Growth/FGPI_Gss_3.txt',STATUS='OLD')
            read(404,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(404, *) Fg_C(i,1)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(404)	
			
			Fg_M(:,:)=0.0D0
	        open(405,FILE='D:/DBGuan/LianTian/RatHeart/AbaqusFile/Growth/FGPI_Gnn_3.txt',STATUS='OLD')
            read(405,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(405, *) Fg_M(i,1)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(405)	
	

      return 
      end
	  
 
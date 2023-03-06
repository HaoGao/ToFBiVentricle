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
      real*8   FrG(170000,9),FrM(170000,9),FrC(170000,9),Fg_G(170000,9),Fg_M(170000,9),Fg_C(170000,9),
     *         Fg_sum(170000,4),Fd(170000,12),VFrsG(3,3),VFrsM(3,3),VFrsC(3,3)
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
	  p_mm= props(12)
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
c-       write(*,*)Time,Nelement(k)-14,Nelement(k)
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
			
				VFrsG(1,1)=FrG(Nelement(k)-4,1)
				VFrsG(2,2)=FrG(Nelement(k)-4,2)
				VFrsG(3,3)=FrG(Nelement(k)-4,3)
				VFrsG(1,2)=FrG(Nelement(k)-4,4)
				VFrsG(2,3)=FrG(Nelement(k)-4,5)
				VFrsG(3,1)=FrG(Nelement(k)-4,6)
				VFrsG(2,1)=FrG(Nelement(k)-4,7)
				VFrsG(3,2)=FrG(Nelement(k)-4,8)
				VFrsG(1,3)=FrG(Nelement(k)-4,9)
				
				VFrsM(1,1)=FrM(Nelement(k)-4,1)
				VFrsM(2,2)=FrM(Nelement(k)-4,2)
				VFrsM(3,3)=FrM(Nelement(k)-4,3)
				VFrsM(1,2)=FrM(Nelement(k)-4,4)
				VFrsM(2,3)=FrM(Nelement(k)-4,5)
				VFrsM(3,1)=FrM(Nelement(k)-4,6)
				VFrsM(2,1)=FrM(Nelement(k)-4,7)
				VFrsM(3,2)=FrM(Nelement(k)-4,8)
				VFrsM(1,3)=FrM(Nelement(k)-4,9)
				
				VFrsC(1,1)=FrC(Nelement(k)-4,1)
				VFrsC(2,2)=FrC(Nelement(k)-4,2)
				VFrsC(3,3)=FrC(Nelement(k)-4,3)
				VFrsC(1,2)=FrC(Nelement(k)-4,4)
				VFrsC(2,3)=FrC(Nelement(k)-4,5)
				VFrsC(3,1)=FrC(Nelement(k)-4,6)
				VFrsC(2,1)=FrC(Nelement(k)-4,7)
				VFrsC(3,2)=FrC(Nelement(k)-4,8)
				VFrsC(1,3)=FrC(Nelement(k)-4,9)
c- growth tensor#		

c- start new method at here  ground matrix
				FgGI(1,1)=Fg_G(Nelement(k)-4,1)
				FgGI(2,2)=Fg_G(Nelement(k)-4,2)
				FgGI(3,3)=Fg_G(Nelement(k)-4,3)
				FgGI(1,2)=Fg_G(Nelement(k)-4,4)
				FgGI(2,3)=Fg_G(Nelement(k)-4,5)
				FgGI(3,1)=Fg_G(Nelement(k)-4,6)
				FgGI(2,1)=Fg_G(Nelement(k)-4,7)
				FgGI(3,2)=Fg_G(Nelement(k)-4,8)
				FgGI(1,3)=Fg_G(Nelement(k)-4,9)
				
c- start new method at here  Myofibre
				FgMI(1,1)=Fg_M(Nelement(k)-4,1)
				FgMI(2,2)=Fg_M(Nelement(k)-4,2)
				FgMI(3,3)=Fg_M(Nelement(k)-4,3)
				FgMI(1,2)=Fg_M(Nelement(k)-4,4)
				FgMI(2,3)=Fg_M(Nelement(k)-4,5)
				FgMI(3,1)=Fg_M(Nelement(k)-4,6)
				FgMI(2,1)=Fg_M(Nelement(k)-4,7)
				FgMI(3,2)=Fg_M(Nelement(k)-4,8)
				FgMI(1,3)=Fg_M(Nelement(k)-4,9)
		
c- start new method at here  Collagen fibre
				FgCI(1,1)=Fg_C(Nelement(k)-4,1)
				FgCI(2,2)=Fg_C(Nelement(k)-4,2)
				FgCI(3,3)=Fg_C(Nelement(k)-4,3)
				FgCI(1,2)=Fg_C(Nelement(k)-4,4)
				FgCI(2,3)=Fg_C(Nelement(k)-4,5)
				FgCI(3,1)=Fg_C(Nelement(k)-4,6)
				FgCI(2,1)=Fg_C(Nelement(k)-4,7)
				FgCI(3,2)=Fg_C(Nelement(k)-4,8)
				FgCI(1,3)=Fg_C(Nelement(k)-4,9)
			
	
			
			
			
			p_g=Fg_sum(Nelement(k)-4,1)
			p_m=Fg_sum(Nelement(k)-4,2)
			p_c=Fg_sum(Nelement(k)-4,3)
			p_fgt=Fg_sum(Nelement(k)-4,4)
			
C--		write(*,*) Time, Nelement(k), Nelement(k)-4, p_g
C--Identity tensor	
			pIden=0.0d0
			pIden(1,1)=1.0d0
			pIden(2,2)=1.0d0
			pIden(3,3)=1.0d0


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
		
		do i=1, ndir
		  do j=1, ndir
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+FI(i,ij)*UI(ij,j)
			end do
			Ro(i,j)=tmp
		  end do
		end do 

c  ***********************************************************
c  ***********************************************************
C- For ground matrix

		do i=1, ndir
			do j=1, ndir
				tmp=0.0d0
				do ia=1, ndir
					do ib=1, ndir
						tmp=tmp+FI(i,ia)*VFrsG(ia,ib)*FgGI(ib,j)
					end do
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

c-------------------------------------
c-------------------------------------		
C-Compute the C, B, f tensor f, s tensor s	ground matrix		
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
		  end do
		end do 

C--Invariant 
		        p_I1=C(1,1)+C(2,2)+C(3,3)	
C--Cauchy stress	
		do i=1, ndir
			do j=1, ndir
				sig(i,j)=p_g*P_a*exp(P_b*(p_I1-3.0d0))*(B(i,j)-p_I1/3.0d0*pIden(i,j))/DET
			end do
		end do  

c  ***********************************************************
c  ***********************************************************
c-------------------------------------	
C- For myofibre

		do i=1, ndir
			do j=1, ndir
				tmp=0.0d0
				do ia=1, ndir
					do ib=1, ndir
						tmp=tmp+FI(i,ia)*VFrsM(ia,ib)*FgMI(ib,j)
					end do
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
		
		ff(1)=Fd(Nelement(k)-4,1)
		ff(2)=Fd(Nelement(k)-4,2)
		ff(3)=Fd(Nelement(k)-4,3)
		
		fff(1)=Fb(1,1)*ff(1)+Fb(1,2)*ff(2)+Fb(1,3)*ff(3)
		fff(2)=Fb(2,1)*ff(1)+Fb(2,2)*ff(2)+Fb(2,3)*ff(3)
		fff(3)=Fb(3,1)*ff(1)+Fb(3,2)*ff(2)+Fb(3,3)*ff(3)
c-------------------------------------
c-------------------------------------		
C-Compute the C, B, f tensor f, s tensor s			
		do i=1, ndir
		  do j=1, ndir
		  
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+Fb(ij,i)*Fb(ij,j)
			end do
			C(i,j)=tmp
			
            ftf(i,j)=fff(i)*fff(j)
c
		  end do
		end do 

C--Invariant 	
				tmp1=0.0d0	
				do i=1,ndir
					do j=1, ndir
						tmp1=tmp1+ff(i)*C(i,j)*ff(j)
					end do
				end do 
				p_I4f=tmp1
				p_I4fmax=max(p_I4f,1.0d0)

                if (Time .GE. 1.0d0 .AND. Time .LE. 1.7d0) then
				    stateNew(k,2)=max(stateOld(k,2),sqrt(p_I4f))
				else
					stateNew(k,2)=sqrt(p_I4f)
				end if

C--Cauchy stress	
		do i=1, ndir
			do j=1, ndir
				sig(i,j)=sig(i,j)+p_m*2.0d0*P_am*(p_I4fmax-1.0d0)*exp(P_bm*(p_I4fmax-1.0d0)**2)*
     *          	(ftf(i,j)-p_I4fmax/3.0d0*pIden(i,j))/DET
			end do
		end do  

c  *****ACITVE  ******************************************************

      t=0.0d0
	  p_lr=0.00185d0
      p_lff=p_lr*sqrt(p_I4f)
	  p_str=0.0d0
		 
      if (tempOld(k).GE.a_thr .AND. stateOld(k,5).LT.asmall) then
	     stateNew(k,5)=Time
      elseif (stateOld(k,5).GE.asmall) then
	     t=Time-stateOld(k,5)
	     stateNew(k,5)=stateOld(k,5)
      else
	     t=0.0d0
      end if 
	  
c 	  - Ta

         t_r=p_mm*p_lff+b_len
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
			stateNew(k,5)=zero
         end if
		 	

         if(t .GT. zero .AND. p_lff .GT. p_l0)then
            term1=exp(b_sha*(p_lff-p_l0))
            Eca=Ca0_m/sqrt(term1-1.0d0)
			term3=Ca0**2+Eca**2
            term2=Ca0**2/term3

            Ta=half*Tmax*term2*(1.0d0-cos(w_t))

            stateNew(k,6)=Ta
	
c-------------------------------------
			do i=1, ndir
				do j=1, ndir
					sigm(i,j)=p_m*Ta/p_I4f*ftf(i,j)
					sig(i,j)=sig(i,j)+p_m*Ta/p_I4f*ftf(i,j)
				end do
			end do
			p_str=sigm(1,1)+sigm(2,2)+sigm(3,3)
         end if
		 
		if (Time .GE. 0.5d0 .AND. Time .LE. 1.0d0) then
		     stateNew(k,4)=max(stateOld(k,4),p_str)
		 else
			 stateNew(k,4)=p_str
		 end if
c  ***********************************************************
c  ***********************************************************
c-------------------------------------	
C- For collagen fibre

		do i=1, ndir
			do j=1, ndir
				tmp=0.0d0
				do ia=1, ndir
					do ib=1, ndir
						tmp=tmp+FI(i,ia)*VFrsC(ia,ib)*FgCI(ib,j)
					end do
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
		
		ff(1)=Fd(Nelement(k)-4,4)
		ff(2)=Fd(Nelement(k)-4,5)
		ff(3)=Fd(Nelement(k)-4,6)
		
		fs(1)=Fd(Nelement(k)-4,7)
		fs(2)=Fd(Nelement(k)-4,8)
		fs(3)=Fd(Nelement(k)-4,9)
		
		fn(1)=Fd(Nelement(k)-4,10)
		fn(2)=Fd(Nelement(k)-4,11)
		fn(3)=Fd(Nelement(k)-4,12)
		
		fff(1)=Fb(1,1)*ff(1)+Fb(1,2)*ff(2)+Fb(1,3)*ff(3)
		fff(2)=Fb(2,1)*ff(1)+Fb(2,2)*ff(2)+Fb(2,3)*ff(3)
		fff(3)=Fb(3,1)*ff(1)+Fb(3,2)*ff(2)+Fb(3,3)*ff(3)
		
		fss(1)=Fb(1,1)*fs(1)+Fb(1,2)*fs(2)+Fb(1,3)*fs(3)
		fss(2)=Fb(2,1)*fs(1)+Fb(2,2)*fs(2)+Fb(2,3)*fs(3)
		fss(3)=Fb(3,1)*fs(1)+Fb(3,2)*fs(2)+Fb(3,3)*fs(3)
		
		fnn(1)=Fb(1,1)*fn(1)+Fb(1,2)*fn(2)+Fb(1,3)*fn(3)
		fnn(2)=Fb(2,1)*fn(1)+Fb(2,2)*fn(2)+Fb(2,3)*fn(3)
		fnn(3)=Fb(3,1)*fn(1)+Fb(3,2)*fn(2)+Fb(3,3)*fn(3)
c-------------------------------------
c-------------------------------------		
C-Compute the C, B, f tensor f, s tensor s			
		do i=1, ndir
		  do j=1, ndir
					
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+Fb(ij,i)*Fb(ij,j)
			end do
			C(i,j)=tmp
			
            ftf(i,j)=fff(i)*fff(j)
			sts(i,j)=fss(i)*fss(j)
			pntn(i,j)=fnn(i)*fnn(j)
c
		  end do
		end do 

C--Invariant 		
				tmp1=0.0d0	
				tmp2=0.0d0	
				tmp3=0.0d0	
				do i=1,ndir
					do j=1, ndir
						tmp1=tmp1+ff(i)*C(i,j)*ff(j)
						tmp2=tmp2+fs(i)*C(i,j)*fs(j)
						tmp3=tmp3+fn(i)*C(i,j)*fn(j)
					end do
				end do 
				p_I4f=tmp1
				p_I4s=tmp2
				p_I4n=tmp3
   
				p_I4fc=max(p_I4f/P_lamd/P_lamd,1.0d0)
				p_I4sc=max(p_I4s/P_lamd/P_lamd,1.0d0)
				p_I4nc=max(p_I4n/P_lamd/P_lamd,1.0d0)

               	if (Time .GE. 1.0d0 .AND. Time .LE. 1.7d0) then
				    stateNew(k,3)=max(stateOld(k,3),sqrt(p_I4f))
				else
					stateNew(k,3)=sqrt(p_I4f)
				end if


C--Cauchy stress	
		do i=1, ndir
			do j=1, ndir
				sig(i,j)=sig(i,j)+p_c*2.0d0*P_ac*(p_I4fc-1.0d0)*exp(P_bc*(p_I4fc-1.0d0)**2)*
     *          	(ftf(i,j)-p_I4fc/3.0d0*pIden(i,j))/DET
     *        		+p_c*2.0d0*P_acsn*(p_I4sc-1.0d0)*exp(P_bcsn*(p_I4sc-1.0d0)**2)*
     *					(sts(i,j)-p_I4sc/3.0d0*pIden(i,j))/DET
     *        		+p_c*2.0d0*P_acsn*(p_I4nc-1.0d0)*exp(P_bcsn*(p_I4nc-1.0d0)**2)*
     *					(pntn(i,j)-p_I4nc/3.0d0*pIden(i,j))/DET
			end do
		end do  

c  ***********************************************************
c-------------------------------------	
C- For volume part
		DET = FI(1,1)*FI(2,2)*FI(3,3) 
     *      - FI(1,1)*FI(2,3)*FI(3,2)  
     *      - FI(1,2)*FI(2,1)*FI(3,3)  
     *      + FI(1,2)*FI(2,3)*FI(3,1)  
     *      + FI(1,3)*FI(2,1)*FI(3,2)  
     *      - FI(1,3)*FI(2,2)*FI(3,1)
		DET=DET/p_fgt
		
	
		
			do i=1, ndir
				do j=1, ndir
					sig(i,j)=sig(i,j)+(DET-1.0d0/DET)*dinv*pIden(i,j)
				end do
			end do

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
      real*8    FrG(170000,9),FrM(170000,9),FrC(170000,9),Fg_G(170000,9),Fg_M(170000,9),Fg_C(170000,9),
     *          Fg_sum(170000,4),Fd(170000,12)
	  integer*4 SDVelement,StartAnalysis,SDVnode
      common /FibreSheetblk/FrG,FrM,FrC,Fg_G,Fg_M,Fg_C,Fg_sum,Fd,SDVelement,SDVnode,StartAnalysis
	  
	    StartAnalysis=StartAnalysis+1
	  	
		  
C- Displacement DATA, out order

			FrG(:,:)=0.0D0
	        open(401,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/Fr_G.txt',STATUS='OLD')
            read(401,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(401, *) FrG(i,1),FrG(i,2),FrG(i,3),FrG(i,4),FrG(i,5),FrG(i,6),FrG(i,7),FrG(i,8),FrG(i,9)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(401)

			FrM(:,:)=0.0D0
	        open(409,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/Fr_M.txt',STATUS='OLD')
            read(409,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(409, *) FrM(i,1),FrM(i,2),FrM(i,3),FrM(i,4),FrM(i,5),FrM(i,6),FrM(i,7),FrM(i,8),FrM(i,9)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(409)
			
			FrC(:,:)=0.0D0
	        open(408,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/Fr_C.txt',STATUS='OLD')
            read(408,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(408, *) FrC(i,1),FrC(i,2),FrC(i,3),FrC(i,4),FrC(i,5),FrC(i,6),FrC(i,7),FrC(i,8),FrC(i,9)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(408)


			Fg_M(:,:)=0.0D0
	        open(402,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/FGPI_M.txt',STATUS='OLD')
            read(402,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(402, *) Fg_M(i,1),Fg_M(i,2),Fg_M(i,3),Fg_M(i,4),Fg_M(i,5),Fg_M(i,6),Fg_M(i,7),Fg_M(i,8),Fg_M(i,9)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(402)	


			Fg_G(:,:)=0.0D0
	        open(403,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/FGPI_G.txt',STATUS='OLD')
            read(403,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(403, *) Fg_G(i,1),Fg_G(i,2),Fg_G(i,3),Fg_G(i,4),Fg_G(i,5),Fg_G(i,6),Fg_G(i,7),Fg_G(i,8),Fg_G(i,9)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(403)	


			Fg_C(:,:)=0.0D0
	        open(404,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/FGPI_C.txt',STATUS='OLD')
            read(404,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(404, *) Fg_C(i,1),Fg_C(i,2),Fg_C(i,3),Fg_C(i,4),Fg_C(i,5),Fg_C(i,6),Fg_C(i,7),Fg_C(i,8),Fg_C(i,9)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(404)	

			Fg_sum(:,:)=0.0D0
	        open(405,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/FGPI_sum.txt',STATUS='OLD')
            read(405,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(405, *) Fg_sum(i,1),Fg_sum(i,2),Fg_sum(i,3),Fg_sum(i,4)
c-				write(*, *) Fg_sum(i,1),Fg_sum(i,2),Fg_sum(i,3),Fg_sum(i,4)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(405)	


			Fd(:,:)=0.0D0
	        open(410,FILE='/home/pgrad2/2306902g/DBGuan/DB/UMAT/Growth_heart/Con_Eccentric_Elastic/FibreD.txt',STATUS='OLD')
            read(410,*) SDVelement
            write(*,*) 'total data to read: n = ', SDVelement

            do i = 1 , SDVelement
                read(410, *) Fd(i,1),Fd(i,2),Fd(i,3),Fd(i,4),Fd(i,5),Fd(i,6),Fd(i,7),Fd(i,8),Fd(i,9),Fd(i,10),Fd(i,11),Fd(i,12)
c-		        write(*, *) Fd(i,1),Fd(i,2),Fd(i,3),Fd(i,4),Fd(i,5),Fd(i,6),Fd(i,7),Fd(i,8),Fd(i,9),Fd(i,10),Fd(i,11),Fd(i,12)
            end do 
            write(*,*) 'read ,Udata value finished!'
            close(410)			
			

      return 
      end
	  
 

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
      real*8   FrG(170000,9),FrM(170000,9),FrC(170000,9),Fg_G(170000,9),Fg_M(170000,9),Fg_C(170000,9),
     *         Fg_sum(170000,4),Fd(170000,12),VFrsG(3,3),VFrsM(3,3),VFrsC(3,3)
	  integer*4 SDVelement,StartAnalysis,SDVnode
      common /FibreSheetblk/FrG,FrM,FrC,Fg_G,Fg_M,Fg_C,Fg_sum,Fd,SDVelement,SDVnode,StartAnalysis
	  
C
C Read material properties
C  
      P_a  = props(1)
      P_b  = props(2)
      P_af = props(3)
      P_bf = props(4)
	  P_an = props(5)
      P_bn = props(6)
	  P_afs = props(7)
      P_bfs = props(8)
	  P_D  = props(9)
      dinv=1/P_D
	  
	  Time=totalTime
C - active 	 
	  t0= props(10)
	  p_mm= props(11)
	  b_len= props(12)
	  p_l0= props(13)
	  b_sha= props(14)
	  Ca0_m= props(15)
	  Ca0= props(16)
	  Tmax= props(17)
	  a_thr= props(18)
C -	  p_nf= props(19)
C -	  p_ns= props(20)
C -	  p_nn= props(21)
	  
C--Define initial growth values	  
C--stateOld(k,1) for growth tensor,stateOld(k,2) for max stretch
C--write(*,*)Time,Nelement(k)-14,stateOld(k,1),StartAnalysis
		
c-------------------------------------
c-------------------------------------	  
C
C Please note here, initial value is not zero
      do k = 1, nblock
c-       write(*,*)Time,Nelement(k)-14,Nelement(k)
C--Keep previous values in stepTime
C--	  write(*,*) Time, Nelement(k)
		
C--Compute deformation gradient tensor F
			Fe(1,1)=defgradNew(k,1)
			Fe(2,2)=defgradNew(k,2)
			Fe(3,3)=defgradNew(k,3)
			Fe(1,2)=defgradNew(k,4)
			Fe(2,3)=defgradNew(k,5)
			Fe(3,1)=defgradNew(k,6)
			Fe(2,1)=defgradNew(k,7)
			Fe(3,2)=defgradNew(k,8)
			Fe(1,3)=defgradNew(k,9)
			
			
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
				tmp=tmp+Fe(i,ij)*UI(ij,j)
			end do
			Ro(i,j)=tmp
		  end do
		end do 

c  ***********************************************************
c  ***********************************************************
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
		
		fff(1)=Fb(1,1)
		fff(2)=Fb(2,1)
		fff(3)=Fb(3,1)		
		fss(1)=Fb(1,2)
		fss(2)=Fb(2,2)
		fss(3)=Fb(3,2)
		fnn(1)=Fb(1,3)
		fnn(2)=Fb(2,3)
		fnn(3)=Fb(3,3)
		
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
			ftf(i,j)=fff(i)*fff(j)
			pntn(i,j)=fnn(i)*fnn(j)
			ftstf(i,j)=fff(i)*fss(j)+fss(i)*fff(j)
		  end do
		end do 


C--Invariant 
		p_I1=C(1,1)+C(2,2)+C(3,3)	
		p_I4f=C(1,1)
		p_I4fmax=max(p_I4f,1.0d0)
		p_I4n=C(3,3)
		p_I4nc=max(p_I4n,1.0d0)
		p_I8fs=C(1,2)
C--Cauchy stress	
		do i=1, ndir
			do j=1, ndir
				sig(i,j)=P_a*exp(P_b*(p_I1-3.0d0))*(B(i,j)-p_I1/3.0d0*pIden(i,j))/DET
     *			+2.0d0*P_af*(p_I4fmax-1.0d0)*exp(P_bf*(p_I4fmax-1.0d0)**2)*
     *				(ftf(i,j)-p_I4fmax/3.0d0*pIden(i,j))/DET
     *			+2.0d0*P_an*(p_I4nc-1.0d0)*exp(P_bn*(p_I4nc-1.0d0)**2)*
     *				(pntn(i,j)-p_I4nc/3.0d0*pIden(i,j))/DET
     *			+P_afs*p_I8fs*exp(P_bfs*p_I8fs**2)*
     *				(ftstf(i,j)-2.0d0*p_I8fs/3.0d0*pIden(i,j))/DET
			end do
		end do  


c  ***********************************************************
c  *****ACITVE  ******************************************************

      t=0.0d0
	  p_lr=0.00185d0
      p_lff=p_lr*sqrt(p_I4f)
	  p_str=0.0d0
	  stateNew(k,3)=tempOld(k)
		 
      if (tempOld(k).GE.a_thr .AND. stateOld(k,1).LT.asmall) then
	     stateNew(k,1)=Time
      elseif (stateOld(k,1).GE.asmall) then
	     t=Time-stateOld(k,1)
	     stateNew(k,1)=stateOld(k,1)
      else
	     t=0.0d0
      end if 
	  
c 	  - Ta

         t_r=p_mm*p_lff+b_len
		 if (t_r .LT. 0.0d0) then
			t_r=asmall
		 end if
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
			stateNew(k,1)=zero
         end if
		 	

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
					sigm(i,j)=Ta/p_I4f*ftf(i,j)
					sig(i,j)=sig(i,j)+Ta/p_I4f*ftf(i,j)
				end do
			end do
			p_str=sigm(1,1)+sigm(2,2)+sigm(3,3)
         end if
		 
c  ***********************************************************
c  ***********************************************************
c-------------------------------------	
	
		
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

	  
 

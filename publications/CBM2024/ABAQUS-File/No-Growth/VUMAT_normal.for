	 !DEC$ ATTRIBUTES ALIAS:"vumat"::VUMAT 
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
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
     2  enerInternNew(nblock), enerInelasNew(nblock),C(ndir,ndir),
     *  F(ndir,ndir),B(ndir,ndir),ftf(ndir,ndir),sts(ndir,ndir),ftfa(ndir,ndir),
     *  pntn(ndir,ndir),ftstf(ndir,ndir),sig(ndir,ndir),ftntf(ndir,ndir),U(ndir,ndir),UI(ndir,ndir),Ro(ndir,ndir),
     *  Fb(ndir, ndir), Fg(ndir,ndir), pIden(ndir,ndir),signew(ndir, ndir),p_value(ndir,ndir),
     *  f0(ndir),fNn(ndir),fn0(ndir),sigmaa(ndir, ndir),RT(ndir, ndir),fNnNn(ndir,ndir),
     *  fv(ndir),fvf(ndir),fcar(ndir),stnts(ndir,ndir)
C  
      character*80 cmname
	  
      parameter ( half = 0.5d0,zero = 0.0d0, one  = 1.d0, 
     *            two  = 2.d0, three = 3.d0,index_J = 3,
     *            asmall  = 2.d-16,pi=3.1415926d0)	 
      real*8 FibreSheet(10000,7), PDFdata(1500,4) 
	  integer*4 PDFelement
      common /FibreSheetblk/FibreSheet, PDFdata, PDFelement
	  
	  
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
	  
	  

C
C Please note here, initial value is not zero
      do k = 1, nblock

C--Compute deformation gradient tensor F
		F(1,1)=defgradNew(k,1)
		F(2,2)=defgradNew(k,2)
		F(3,3)=defgradNew(k,3)
		F(1,2)=defgradNew(k,4)
		F(2,3)=defgradNew(k,5)
		F(3,1)=defgradNew(k,6)
		F(2,1)=defgradNew(k,7)
		F(3,2)=defgradNew(k,8)
		F(1,3)=defgradNew(k,9)
	

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
		
C--Compute the determinant 	 
		DET = F(1,1)*F(2,2)*F(3,3) 
     *      - F(1,1)*F(2,3)*F(3,2)  
     *      - F(1,2)*F(2,1)*F(3,3)  
     *      + F(1,2)*F(2,3)*F(3,1)  
     *      + F(1,3)*F(2,1)*F(3,2)  
     *      - F(1,3)*F(2,2)*F(3,1)

C- F is decomposed
        do i=1, ndir
			do j=1, ndir			
			Fb(i,j)=F(i,j)*DET**(-1.0/3.0)
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
				tmp=tmp+F(i,ij)*UI(ij,j)
			end do
			Ro(i,j)=tmp
			
			tmp=0.0d0
			pIden(i,j)=tmp
			sigmaa(i,j)=tmp 
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
				sig(i,j)=P_a*exp(P_b*(p_I1-3.0d0))*
     *					(B(i,j)-p_I1/3.0d0*pIden(i,j))
     *        		+2.0d0*P_am*(p_I4fm-1.0d0)*exp(P_bm*(p_I4fm-1.0d0)**2)*
     *					(ftf(i,j)-p_I4fm/3.0d0*pIden(i,j))
     *        		+2.0d0*P_ac*(p_I4fc-1.0d0)*exp(P_bc*(p_I4fc-1.0d0)**2)*
     *					(ftf(i,j)-p_I4fc/3.0d0*pIden(i,j))
     *        		+2.0d0*P_acsn*(p_I4sc-1.0d0)*exp(P_bcsn*(p_I4sc-1.0d0)**2)*
     *					(sts(i,j)-p_I4sc/3.0d0*pIden(i,j))
     *        		+2.0d0*P_acsn*(p_I4nc-1.0d0)*exp(P_bcsn*(p_I4nc-1.0d0)**2)*
     *					(pntn(i,j)-p_I4nc/3.0d0*pIden(i,j))
     *				+(DET**2-1.0d0)*dinv*pIden(i,j)
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
		 
		 if (t_r .GT. t0) then
			t_r=t_r/t0*0.001d0+t0
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
			stateNew(k,4)=zero
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
					sig(i,j)=sig(i,j)/DET+p_nf*Ta/p_I4f*ftf(i,j)+p_ns*Ta/p_I4s*sts(i,j)+p_nn*Ta/p_I4n*pntn(i,j)
				end do
			end do
			 
         end if
		  
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

			stressPower = half * (
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
C
      return
      end
	  

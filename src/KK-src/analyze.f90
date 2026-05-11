!Compile with AUTODOUBLE ON!
    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535d0,grav=9.81d0
    INTEGER,PARAMETER:: NFsmooth=100
    
    REAL*8,ALLOCATABLE:: dseisn(:,:),dseise(:,:),dseisz(:,:)
    REAL*8,ALLOCATABLE:: vseisn(:,:),vseise(:,:),vseisz(:,:)
    REAL*8,ALLOCATABLE:: raseisn(:,:),raseise(:,:),raseisz(:,:)
    COMPLEX*16,allocatable:: aseisn(:,:),aseise(:,:),aseisz(:,:)
    REAL*8,allocatable:: asseisn(:,:),asseise(:,:),asseisz(:,:),freqaxis(:,:),smoothspecn(:,:),smoothspece(:,:),smoothspecz(:,:)
    REAL*8,ALLOCATABLE:: pgd(:,:),pgv(:,:),pga(:,:),Hn(:),He(:),Hz(:),AriasN(:),AriasE(:),AriasZ(:)
    REAL*8,ALLOCATABLE:: pgdrotd50(:),pgvrotd50(:),pgarotd50(:)
    REAL*8 dum,freq,df,xsta,ysta,jbdist,rotd50
    REAL*8 dt,TL,fc1,fc2,fc3,fc4
    REAL*8,ALLOCATABLE:: PSAn(:,:),PSAe(:,:),PSAz(:,:),PSVn(:,:),PSVe(:,:),PSVz(:,:),PSFreq(:),dumPS(:)
    INTEGER NS,FFT,FFTn,NR,NFMAX,PSN,NRS
    INTEGER i,j,k

    open(101,FILE='rdseisn.dat')
    open(103,FILE='rdseise.dat')
    open(105,FILE='rdseisz.dat')
    open(102,FILE='rdsseisn.dat')
    open(104,FILE='rdsseise.dat')
    open(106,FILE='rdsseisz.dat')
    open(111,FILE='vseisn.dat')
    open(113,FILE='vseise.dat')
    open(115,FILE='vseisz.dat')
    open(112,FILE='rvsseisn.dat')
    open(114,FILE='rvsseise.dat')
    open(116,FILE='rvsseisz.dat')
    open(121,FILE='raseisn.dat')
    open(123,FILE='raseise.dat')
    open(125,FILE='raseisz.dat')
    open(122,FILE='rasseisn.dat')
    open(124,FILE='rasseise.dat')
    open(126,FILE='rasseisz.dat')
    open(127,FILE='rassseisn.dat')
    open(128,FILE='rassseise.dat')
    open(129,FILE='rassseisz.dat')

    open(131,FILE='rpgd.dat')
    open(132,FILE='rpgv.dat')
    open(133,FILE='rpga.dat')
    open(134,FILE='rhousner.dat')
    open(135,FILE='rarias.dat')
    open(136,FILE='rpga2v.dat')

    open(141,FILE='rapseisn.dat')
    open(142,FILE='rapseise.dat')
    open(143,FILE='rapseisz.dat')
    open(151,FILE='rvpseisn.dat')
    open(152,FILE='rvpseise.dat')
    open(153,FILE='rvpseisz.dat')

    open(100,FILE='analyze.in')
    read(100,*);read(100,*)NRS
    close(100)
    
    open(100,FILE='input.dat')

    read(100,*)
    read(100,*)NFMAX,NFMAX
    read(100,*)
    read(100,*)TL
    df=1.d0/TL
    read(100,*)
    read(100,*)NS
    do i=1,15
      read(100,*)
    enddo
    read(100,*)FFT
    FFTn=int(log(real(FFT,8))/log(2.d0)+.99999d0)
    dt=TL/real(FFT,8)
    do i=1,5
      read(100,*)
    enddo
    read(100,*)NR
    NR=NR*NRS
    do i=1,6
      read(100,*)
    enddo
    read(100,*)fc1,fc4
    
    open(295,FILE='frequencies.txt')
    read(295,*)PSN
    allocate(PSFreq(PSN),dumPS(PSN))
    do i=1,PSN
      read(295,*)PSfreq(i)
    enddo
    PSfreq=1.d0/PSfreq
    
    do i=1,NR
      allocate(vseisn(FFT,NS),vseise(FFT,NS),vseisz(FFT,NS))
      allocate(dseisn(FFT,NS),dseise(FFT,NS),dseisz(FFT,NS))
      allocate(pgd(3,NS),pgv(3,NS),pga(3,NS),Hn(NS),He(NS),Hz(NS))
      allocate(pgdrotd50(NS),pgvrotd50(NS),pgarotd50(NS))
      allocate(AriasN(NS),AriasE(NS),AriasZ(NS))
      allocate(PSAn(PSN,NS),PSAe(PSN,NS),PSAz(PSN,NS),PSVn(PSN,NS),PSVe(PSN,NS),PSVz(PSN,NS))
      allocate(aseisn(FFT,NS),aseise(FFT,NS),aseisz(FFT,NS),asseisn(FFT,NS),asseise(FFT,NS),asseisz(FFT,NS))
      allocate(raseisn(FFT,NS),raseise(FFT,NS),raseisz(FFT,NS))
      allocate(freqaxis(NFsmooth,NS),smoothspecn(Nfsmooth,NS),smoothspece(Nfsmooth,NS),smoothspecz(Nfsmooth,NS))

! Processing velocities

      do j=1,FFT
        read(111,*)dum,vseisn(j,:)
        read(113,*)dum,vseise(j,:)
        read(115,*)dum,vseisz(j,:)
      enddo
      do j=1,NS
        pgv(1,j)=maxval(abs(vseisn(:,j)))
        pgv(2,j)=maxval(abs(vseise(:,j)))
        pgv(3,j)=maxval(abs(vseisz(:,j)))
        pgvrotd50(j)=rotd50(FFT,vseisn(:,j),vseise(:,j))
      enddo
      aseisn=vseisn;aseise=vseise;aseisz=vseisz

!$OMP parallel do private(j) DEFAULT(SHARED)
      do j=1,NS
        CALL FCOOLR(FFTn,aseisn(1:FFT,j),-1.)
        CALL FCOOLR(FFTn,aseise(1:FFT,j),-1.)
        CALL FCOOLR(FFTn,aseisz(1:FFT,j),-1.)
	      aseisn(1:FFT,j)=aseisn(1:FFT,j)*dt
	      aseise(1:FFT,j)=aseise(1:FFT,j)*dt
	      aseisz(1:FFT,j)=aseisz(1:FFT,j)*dt
      enddo
!$OMP end parallel do

      do j=1,NFMAX
        freq=df*(j-1)
        write(112,100)freq,abs(aseisn(j,:))
        write(114,100)freq,abs(aseise(j,:))
        write(116,100)freq,abs(aseisz(j,:))
      enddo

! Processing acceleration

      do j=1,FFT/2+1
        freq=df*(j-1)
        aseisn(j,:)=aseisn(j,:)*cmplx(0.d0,2.d0*PI*freq)
        aseise(j,:)=aseise(j,:)*cmplx(0.d0,2.d0*PI*freq)
        aseisz(j,:)=aseisz(j,:)*cmplx(0.d0,2.d0*PI*freq)
      enddo
      asseisn=abs(aseisn)
      asseise=abs(aseise)
      asseisz=abs(aseisz)
      do j=1,NFMAX
        freq=df*(j-1)
        write(122,100)freq,asseisn(j,:)
        write(124,100)freq,asseise(j,:)
        write(126,100)freq,asseisz(j,:)
      enddo
!$OOOOMP parallel do private(j) DEFAULT(SHARED)
      do j=1,NS
        CALL smoothspectrum(FFT,Nfsmooth,df,fc1,fc4,asseisn(:,j),freqaxis(:,j),smoothspecn(:,j))
        CALL smoothspectrum(FFT,Nfsmooth,df,fc1,fc4,asseise(:,j),freqaxis(:,j),smoothspece(:,j))
        CALL smoothspectrum(FFT,Nfsmooth,df,fc1,fc4,asseisz(:,j),freqaxis(:,j),smoothspecz(:,j))
      enddo
!$OOOOMP end parallel do
      do j=1,Nfsmooth
        write(127,100)freqaxis(j,1),smoothspecn(j,:)
        write(128,100)freqaxis(j,1),smoothspece(j,:)
        write(129,100)freqaxis(j,1),smoothspecz(j,:)
      enddo

!$OMP parallel do private(j) DEFAULT(SHARED)
      do j=1,NS
        aseisn(FFT/2+2:FFT,j)=conjg(aseisn(FFT/2:2:-1,j))
        aseise(FFT/2+2:FFT,j)=conjg(aseise(FFT/2:2:-1,j))
        aseisz(FFT/2+2:FFT,j)=conjg(aseisz(FFT/2:2:-1,j))
        CALL FCOOLR(FFTn,aseisn(1:FFT,j),1.)
        CALL FCOOLR(FFTn,aseise(1:FFT,j),1.)
        CALL FCOOLR(FFTn,aseisz(1:FFT,j),1.)
        raseisn(1:FFT,j)=real(aseisn(1:FFT,j),8)*df
        raseise(1:FFT,j)=real(aseise(1:FFT,j),8)*df
        raseisz(1:FFT,j)=real(aseisz(1:FFT,j),8)*df
      enddo
!$OMP end parallel do

      do j=1,FFT  !int(20.d0/dt),int(80.d0/dt)
        write(121,100)dt*(j-1),raseisn(j,:)
        write(123,100)dt*(j-1),raseise(j,:)
        write(125,100)dt*(j-1),raseisz(j,:)
      enddo

!$OMP parallel do private(j) DEFAULT(SHARED)
      do j=1,NS
        call PCN05(FFT,FFT,PSN,PSN,dt,.05d0,PSfreq,raseisn(1:FFT,j),dumPS(:),PSVn(:,j),PSAn(:,j),dumPS(:),dumPS(:))
        call PCN05(FFT,FFT,PSN,PSN,dt,.05d0,PSfreq,raseise(1:FFT,j),dumPS(:),PSVe(:,j),PSAe(:,j),dumPS(:),dumPS(:))
        call PCN05(FFT,FFT,PSN,PSN,dt,.05d0,PSfreq,raseisz(1:FFT,j),dumPS(:),PSVz(:,j),PSAz(:,j),dumPS(:),dumPS(:))
        call HOUSNER(raseisn(1:FFT,j),FFT,dt,.05d0,Hn(j))
        call HOUSNER(raseise(1:FFT,j),FFT,dt,.05d0,He(j))
        call HOUSNER(raseisz(1:FFT,j),FFT,dt,.05d0,Hz(j))
        pga(1,j)=maxval(abs(raseisn(:,j)))
        pga(2,j)=maxval(abs(raseise(:,j)))
        pga(3,j)=maxval(abs(raseisz(:,j)))
        pgarotd50(j)=rotd50(FFT,raseisn(:,j),raseise(:,j))
        AriasN(j)=PI/2.d0/grav*sum(raseisn(1:FFT,j)**2)
        AriasE(j)=PI/2.d0/grav*sum(raseise(1:FFT,j)**2)
        AriasZ(j)=PI/2.d0/grav*sum(raseisz(1:FFT,j)**2)
      enddo
!$OMP end parallel do

      do j=1,PSN
        write(141,100)1.d0/PSfreq(j),PSAn(j,:)
        write(142,100)1.d0/PSfreq(j),PSAe(j,:)
        write(143,100)1.d0/PSfreq(j),PSAz(j,:)
        write(151,100)1.d0/PSfreq(j),PSVn(j,:)
        write(152,100)1.d0/PSfreq(j),PSVe(j,:)
        write(153,100)1.d0/PSfreq(j),PSVz(j,:)
      enddo
      
! Processing displacement

!$OMP parallel do private(j,k) DEFAULT(SHARED)
      do j=1,NS
        do k=1,FFT
          dseisn(k,j)=sum(vseisn(1:k,j))*dt
          dseise(k,j)=sum(vseise(1:k,j))*dt
          dseisz(k,j)=sum(vseisz(1:k,j))*dt
        enddo
      enddo
!$OMP end parallel do

      do j=1,FFT
        write(101,100)dt*(j-1),dseisn(j,:)
        write(103,100)dt*(j-1),dseise(j,:)
        write(105,100)dt*(j-1),dseisz(j,:)
      enddo
      do j=1,NS
        pgd(1,j)=maxval(abs(dseisn(:,j)))
        pgd(2,j)=maxval(abs(dseise(:,j)))
        pgd(3,j)=maxval(abs(dseisz(:,j)))
        pgdrotd50(j)=rotd50(FFT,dseisn(:,j),dseise(:,j))
      enddo
      aseisn=dseisn;aseise=dseise;aseisz=dseisz

!$OMP parallel do private(j) DEFAULT(SHARED)
      do j=1,NS
        CALL FCOOLR(FFTn,aseisn(1:FFT,j),-1.)
        CALL FCOOLR(FFTn,aseise(1:FFT,j),-1.)
        CALL FCOOLR(FFTn,aseisz(1:FFT,j),-1.)
        aseisn(1:FFT,j)=aseisn(1:FFT,j)*dt
        aseise(1:FFT,j)=aseise(1:FFT,j)*dt
        aseisz(1:FFT,j)=aseisz(1:FFT,j)*dt
      enddo
!$OMP end parallel do

      do j=1,NFMAX
        freq=df*(j-1)
        write(102,100)freq,abs(aseisn(j,:))
        write(104,100)freq,abs(aseise(j,:))
        write(106,100)freq,abs(aseisz(j,:))
      enddo

! Writing peak values

!      open(99,file='stations.dat')
      open(99,file='stations.JBdist.dat')
      do j=1,NS
!        read(99,*)xsta,ysta
        read(99,*)xsta,ysta,jbdist
        write(131,'(I3,11E13.5)')j,xsta,ysta,pgd(1,j),pgd(2,j),pgd(3,j),sqrt(pgd(1,j)*pgd(2,j)),pgdrotd50(j),&
sqrt(xsta**2+ysta**2),atan2(ysta,xsta),jbdist
        write(132,'(I3,11E13.5)')j,xsta,ysta,pgv(1,j),pgv(2,j),pgv(3,j),sqrt(pgv(1,j)*pgv(2,j)),pgvrotd50(j),&
sqrt(xsta**2+ysta**2),atan2(ysta,xsta),jbdist
        write(133,'(I3,11E13.5)')j,xsta,ysta,pga(1,j),pga(2,j),pga(3,j),sqrt(pga(1,j)*pga(2,j)),pgarotd50(j),&
sqrt(xsta**2+ysta**2),atan2(ysta,xsta),jbdist
        write(134,'(I3,10E13.5)')j,xsta,ysta,Hn(j),He(j),Hz(j),max(Hn(j),He(j)),&
sqrt(xsta**2+ysta**2),atan2(ysta,xsta),jbdist
        write(135,'(I3,10E13.5)')j,xsta,ysta,AriasN(j),AriasE(j),AriasZ(j),max(AriasN(j),AriasE(j)),&
sqrt(xsta**2+ysta**2),atan2(ysta,xsta),jbdist
        write(136,'(I3,10E13.5)')j,xsta,ysta,pga(1,j)/pgv(1,j),pga(2,j)/pgv(2,j),pga(3,j)/pgv(3,j),&
maxval(pga(1:2,j)/pgv(1:2,j)),sqrt(xsta**2+ysta**2),atan2(ysta,xsta),jbdist
      enddo
close(99)

      open(400,file='SAfig10HzH.dat')
      open(444,file='SAfig10HzN.dat')
      open(445,file='SAfig10HzE.dat')
      open(446,file='SAfig10HzZ.dat')

      open(500,file='SAfig3HzH.dat')
      open(555,file='SAfig3HzN.dat')
      open(556,file='SAfig3HzE.dat')
      open(557,file='SAfig3HzZ.dat')

      open(600,file='SAfig05HzH.dat')
      open(666,file='SAfig05HzN.dat')
      open(667,file='SAfig05HzE.dat')
      open(668,file='SAfig05HzZ.dat')

      do j=1,PSN
      if(j==1)then
      write(*,*)1./PSfreq(j)
open(99,file='stations.JBdist.dat')
      do k=1,NS
        read(99,*)xsta,ysta,jbdist
        write(400,*)jbdist,atan2(ysta,xsta),sqrt(PSAn(j,k)*PSAe(j,k))
        write(444,*)jbdist,atan2(ysta,xsta),PSAn(j,k)
	write(445,*)jbdist,atan2(ysta,xsta),PSAe(j,k)
        write(446,*)jbdist,atan2(ysta,xsta),PSAz(j,k)
      enddo
close(99)
      elseif(j==18)then
      write(*,*)1./PSfreq(j)
open(99,file='stations.JBdist.dat')

      do k=1,NS
        read(99,*)xsta,ysta,jbdist
        write(500,*)jbdist,atan2(ysta,xsta),sqrt(PSAn(j,k)*PSAe(j,k))
        write(555,*)jbdist,atan2(ysta,xsta),PSAn(j,k)
	write(556,*)jbdist,atan2(ysta,xsta),PSAe(j,k)
        write(557,*)jbdist,atan2(ysta,xsta),PSAz(j,k)
      enddo
close(99)
      elseif(j==46)then
      write(*,*)1./PSfreq(j)
open(99,file='stations.JBdist.dat')

      do k=1,NS
        read(99,*)xsta,ysta,jbdist
        write(600,*)jbdist,atan2(ysta,xsta),sqrt(PSAn(j,k)*PSAe(j,k))
        write(666,*)jbdist,atan2(ysta,xsta),PSAn(j,k)
	write(667,*)jbdist,atan2(ysta,xsta),PSAe(j,k)
        write(668,*)jbdist,atan2(ysta,xsta),PSAz(j,k)
      enddo
      endif
      enddo

      close(444)
      close(445)
      close(446)
      close(555)
      close(556)
      close(557)
      close(666)
      close(667)
      close(668)
      close(600)
      close(400)
      close(500)
       
          

      close(99)
      deallocate(dseisn,dseise,dseisz)
      deallocate(vseisn,vseise,vseisz)
      deallocate(aseisn,aseise,aseisz,raseisn,raseise,raseisz,asseisn,asseise,asseisz)
      deallocate(pgd,pgv,pga,Hn,He,Hz,PSAn,PSAe,PSAz,PSVn,PSVe,PSVz)
      deallocate(AriasN,AriasE,AriasZ)
      deallocate(freqaxis,smoothspecn,smoothspece,smoothspecz)
      deallocate(pgdrotd50,pgvrotd50,pgarotd50)



    enddo

100 FORMAT(1000E13.5)

    END


    SUBROUTINE smoothspectrum(Nf,Nfsmooth,df,flo,fro,spec,freqaxis,smoothspec)
    !Smoothing spectrum by Konno & Omachi, 1998 BSSA, method
    IMPLICIT NONE
    INTEGER Nf,Nfsmooth
    REAL*8 spec(Nf),WB(Nf/2+1,Nfsmooth),freqaxis(Nfsmooth),smoothspec(Nfsmooth),flo,fro
    REAL*8 freq,df
    INTEGER i,j
    do j=1,Nfsmooth
      freqaxis(j)=10.**((log10(fro)-log10(flo))/real(Nfsmooth-1)*real(j-1)+log10(flo))
      WB(:,j)=0.
      do i=2,Nf/2+1
        freq=df*(i-1)
        if(freq.ne.freqaxis(j))then
          WB(i,j)=(sin(20.*log10(freq/freqaxis(j)))/20./log10(freq/freqaxis(j)))**4
        else
          WB(i,j)=1.
        endif
      enddo
    enddo
    do j=1,Nfsmooth
      smoothspec(j)=sum(abs(spec(1:Nf/2+1))*WB(1:Nf/2+1,j))/sum(WB(1:Nf/2+1,j))
    enddo        
    END
    
    FUNCTION rotd50(NT,a1,a2)
    IMPLICIT NONE
    REAL*8,PARAMETER:: deg2rad=3.1415926535/180.
    INTEGER,PARAMETER:: nrot=180
    REAL*8 rotd50
    INTEGER NT
    REAL*8 a1(NT),a2(NT),arot(NT),peakrot(nrot)
    INTEGER indx(nrot)
    INTEGER n50,irot
    REAL*8 angrot
    
    n50=nrot/2
    do irot=1,nrot
      angrot=(irot-1.)*deg2rad
      arot(:)=a1(:)*cos(angrot)+a2(:)*sin(angrot)
      peakrot(irot)=maxval(abs(arot(:)))
    enddo
    call indexx(nrot,peakrot,indx)
    rotd50=peakrot(indx(n50))    !median
    !rotd50=peakrot(indx(nrot))  !maximum
    
    END
    
    SUBROUTINE HOUSNER(acc,N,dt,DAMP,H)
    IMPLICIT NONE
    integer,PARAMETER:: STEPS=100
    real*8, PARAMETER:: P0=0.1d0,P1=2.5d0
    integer N,i
    real*8 acc(N),H,DAMP,dt,P(STEPS),XSV(STEPS),dumPS(STEPS)
    do i=1,STEPS
      P(i)=(P1-P0)/real(STEPS+1,8)*(real(i,8)+0.5d0)+P0
    enddo
    
    CALL PCN05(N,N,STEPS,STEPS,DT,DAMP,P,ACC,dumPS(:),XSV(:),dumPS(:),dumPS(:),dumPS(:))

    H=sum(XSV)*(P1-P0)/real(STEPS,8)

    END

!SUBROUTINE FCOOLR(K,D,SN)
!FAST FOURIER TRANSFORM OF N = 2**K COMPLEX DATA POINTS
!REPARTS HELD IN D(1,3,...2N-1), IMPARTS IN D(2,4,...2N).

        SUBROUTINE FCOOLR(K,D,SN)
!        REAL*8 INU(20),D(32768)
        REAL*8 INU(20),D(*)
        LX=2**K
        Q1=LX
        IL=LX
        SH=SN*6.28318530718/Q1
        DO 10 I=1,K
        IL=IL/2
10      INU(I)=IL
        NKK=1
        DO 40 LA=1,K
        NCK=NKK
        NKK=NKK+NKK
        LCK=LX/NCK
        L2K=LCK+LCK
        NW=0
        DO 40 ICK=1,NCK
        FNW=NW
        AA=SH*FNW
        W1=COS(AA)
        W2=SIN(AA)
        LS=L2K*(ICK-1)
        DO 20 I=2,LCK,2
        J1=I+LS
        J=J1-1
        JH=J+LCK
        JH1=JH+1
        Q1=D(JH)*W1-D(JH1)*W2
        Q2=D(JH)*W2+D(JH1)*W1
        D(JH)=D(J)-Q1
        D(JH1)=D(J1)-Q2
        D(J)=D(J)+Q1
20      D(J1)=D(J1)+Q2
        DO 29 I=2,K
        ID=INU(I)
        IL=ID+ID
        IF(NW-ID-IL*(NW/IL)) 40,30,30
30      NW=NW-ID
29      CONTINUE
40      NW=NW+ID
        NW=0
        DO 6 J=1,LX
        IF(NW-J) 8,7,7
7       JJ=NW+NW+1
        J1=JJ+1
        JH1=J+J
        JH=JH1-1
        Q1=D(JJ)
        D(JJ)=D(JH)
        D(JH)=Q1
        Q1=D(J1)
        D(J1)=D(JH1)
        D(JH1)=Q1
8       DO 9 I=1,K
        ID=INU(I)
        IL=ID+ID
        IF(NW-ID-IL*(NW/IL)) 6,5,5
5       NW=NW-ID
9       CONTINUE
6       NW=NW+ID
        RETURN
        END



      subroutine sdcomp(accg,na,omn,beta,dt,sd)
! This is a modified version of "Quake.For", written by
! Stavros A. Anagnostopoulos, Oct. 1986.  The formulation is that of
! Nigam and Jennings (BSSA, v. 59, 909-922, 1969).  This modification 
! eliminates the computation of the relative velocity and absolute 
! acceleration; it returns only the relative displacement.  
      real*8 accg(*)
      real*8 omn,beta,dt,sd
      omt=omn*dt
      d2=1-beta*beta
      d2=sqrt(d2)
      bom=beta*omn
      d3=2.*bom
      omd=omn*d2
      om2=omn*omn
      omdt=omd*dt
      c1=1./om2
      c2=2.*beta/(om2*omt)
      c3=c1+c2
      c4=1./(omn*omt)
      ss=sin(omdt)
      cc=cos(omdt)
      bomt=beta*omt
      ee=exp(-bomt)
      ss=ss*ee
      cc=cc*ee
      s1=ss/omd
      s2=s1*bom
      s3=s2+cc
      a11=s3
      a12=s1
      a21=-om2*s1
      a22=cc-s2
      s4=c4*(1.-s3)
      s5=s1*c4+c2
      b11=s3*c3-s5
      b12=-c2*s3+s5-c1
      b21=-s1+s4
      b22=-s4
      sd=0.
      n1=na-1
      y=0.
      ydot=0.
      do 1 i=1,n1
      y1=a11*y+a12*ydot+b11*accg(i)+b12*accg(i+1)
      ydot=a21*y+a22*ydot+b21*accg(i)+b22*accg(i+1)
! next two lines have been added for hardware portability
      if (y1.lt.1.e-30.and.y1.gt.0.) y1=0.
      if (ydot.lt.1.e-30.and.ydot.gt.0.) ydot=0.
      y=y1
      z=abs(y)
      z1=abs(ydot)
      if (z.gt.sd) sd=z
1     continue
      return
      end

!-----------------------------------------------------------------------
      SUBROUTINE PCN05(NMX,N,NPMX,NP,DT,DAMP,P,ACC,XSD,XSV,XSA,XPSV,XPSA)
!-----------------------------------------------------------------------
!      ORIGINAL ONLY GAVE SA. MODIFIED BY JGA TO GIVE SD,SV,SA,PSV.PSA
!      new input:
!      NMX:   dimension of array containing acceleration trace
!      P:     period
!      DAMP:  damping
!      DT:    delta t of acceleration trace
!      ACC:   acceleration trace
!      N:     number of points of acceleration trace
!      NP:    number of frequencies of SDF oscillator
!      output:
!      XSD:   spectral displacement
!      XSV:   spectral velocity
!      XSA:   spectral acceleration
!      XPSV:  pseudo velocity
!      XPSA:  pseudo acceleration
!-----------------------------------------------------------------------
      DIMENSION A(2,2),B(2,2),ACC(NMX),P(NPMX)
      DIMENSION XSD(NPMX),XSV(NPMX),XSA(NPMX),XPSV(NPMX)
      DIMENSION XPSA(NPMX),XFS(NPMX)
!-----------------------------------------------------------------------
      DO 10 IP=1,NP
         W=    6.2832/P(IP)
         DELT= P(IP)/10.
         L=    DT/DELT+1.0-1.E-05
         VERTL=1.0/FLOAT(L)
         DELT= DT*VERTL
         CALL PCN04(DAMP,W,DELT,A,B)
         XIP=    0.0
         XIPD=   0.0
         XSA(IP)=0.0
         XSV(IP)=0.0
         XSD(IP)=0.0
         NN=N-1
         DO 1 J=1,NN
            AI=ACC(J)
            SL=(ACC(J+1)-ACC(J))*VERTL
            DO 2 JJ=1,L
               AF=   AI+SL
               XIP1= XIP
               XIPD1=XIPD
               XIP=  A(1,1)*XIP1+A(1,2)*XIPD1+B(1,1)*AI+B(1,2)*AF
               XIPD= A(2,1)*XIP1+A(2,2)*XIPD1+B(2,1)*AI+B(2,2)*AF
               AI=   AF
               ABSOL=ABS(2.*DAMP*W*XIPD+W*W*XIP)
               IF(ABS(XIP).GT.XSD(IP))  XSD(IP)=ABS(XIP)
               IF(ABS(XIPD).GT.XSV(IP)) XSV(IP)=ABS(XIPD)
               IF(ABSOL.GT.XSA(IP))     XSA(IP)=ABSOL
2              CONTINUE
1           CONTINUE
         XPSV(IP)= XSD(IP)*W
         XPSA(IP)= XSV(IP)*W
         XFS(IP)=  SQRT((W*XIP)**2+XIPD**2)
10       CONTINUE
      RETURN
      END
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE PCN04(D,W,DELT,A,B)
!-----------------------------------------------------------------------
!     called by pcn05
!     D:     damping
!     W:     frequency omega
!     DELT:  something like integration step (?)
!     A:     2x2 matrix, returned
!     B:     2x2 matrix, returned
!-----------------------------------------------------------------------
      DIMENSION A(2,2),B(2,2)
!-----------------------------------------------------------------------
      DW=D*W
      D2=D*D
      A0=EXP(-DW*DELT)
      A1=W*SQRT(1.-D2)
      AD1=A1*DELT
      A2=SIN(AD1)
      A3=COS(AD1)
      A7=1.0/(W*W)
      A4=(2.0*D2-1.0)*A7
      A5=D/W
      A6= 2.0*A5*A7
      A8=1.0/A1
      A9=-(A1*A2+DW*A3)*A0
      A10=(A3-DW*A2*A8)*A0
      A11=A2*A8
      A12=A11*A0
      A13=A0*A3
      A14=A10*A4
      A15=A12*A4
      A16=A6*A13
      A17=A9*A6
      A(1,1)=A0*(DW*A11+A3)
      A(1,2)=A12
      A(2,1)=A10*DW+A9
      A(2,2)=A10
      DINV=1.0/DELT
      B(1,1)=(-A15-A16+A6)*DINV-A12*A5-A7*A13
      B(1,2)=(A15+A16-A6)*DINV+A7
      B(2,1)=(-A14-A17-A7)*DINV-A10*A5-A9*A7
      B(2,2)=(A14+A17+A7)*DINV
      RETURN
    END
!-----------------------------------------------------------------------

    
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

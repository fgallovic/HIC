! Verze pro variace mechanismu

program KingKong

implicit none

integer, parameter :: pulse=4,SUBmin=2,SUB=8 !suby stejne jako prepare2.f90
real, parameter :: twopi=6.28318530717959,PI=3.1415926535d0
real*8, parameter:: rotateto=0
real kappa
complex, parameter :: unit=(0.,1.)
integer pmax
integer ng1,ng2, nfc
integer np
integer nfmax,nfmax2
integer nr,ScenNum
real T
real*8 leng,widt,depth
real Momain
real*8 dip,strike,rake
real*8 hypo(2)
real*8 alfa,beta,rho,acko
real*8, allocatable :: staN(:),staE(:),fault(:,:),JBdist(:)
real*8 JBdistance
real dt,df
real mu,lambda
real*8 vr
integer idum,iidum
complex, allocatable :: cir(:,:,:,:),cir2(:,:,:)
complex, allocatable :: spbrune(:),srcfce(:),stf(:)
complex, allocatable :: cseis(:,:),csseis(:,:)
real, allocatable :: stftime(:)
real, allocatable :: iir(:,:,:,:),rir(:,:,:,:),fltr4(:)
real, allocatable :: rcsxy(:,:),tr(:)
real, allocatable :: x1a(:),x2a(:),y2a(:,:),ry2a(:,:),dyry2a(:),dyy2a(:)
real, allocatable :: Mo(:),fc(:),fc1(:),fc2(:)
real, allocatable :: crossking(:),crosskong(:)
integer FFTcross1,FFTcross2
real KingRiseTime
real grre,grim
real Moir,fco
integer i,j,k,l,m,n,p
integer rec
integer dum
real ddum
real*4 dumarr(6)
real, allocatable :: smogramy(:,:,:,:),pga(:,:),pga2(:,:,:)
real freq,elem,dw,dl,ww,ll,tr2,slippy
integer wi,li
real dum2,xsta,ysta
integer wmax,lmax,submax
real*8 ms,closdist
real, allocatable :: slip(:,:),ruptime(:,:)
real stressdropmain,fcmain
integer, allocatable :: fcsta(:)

logical stanice(1000)
stanice=.FALSE.
!stanice=.TRUE.;stanice(19)=.FALSE.

open(10,file='input.dat',action='read')
read(10,*)
read(10,*) nfmax2,nfmax
read(10,*)
read(10,*) T
read(10,*)
read(10,*) nr
read(10,*)
read(10,*) ng2,ng1,lmax,wmax
read(10,*)
read(10,*)
read(10,*)
read(10,*) Momain
read(10,*)
read(10,*) strike,dip,rake
read(10,*)
read(10,*)depth
read(10,*)
read(10,*) leng,widt
read(10,*)
read(10,*) hypo(2), hypo(1)
read(10,*)
read(10,*) np
read(10,*)
read(10,*) alfa,beta,rho
read(10,*)
read(10,*) vr
read(10,*)
read(10,*) pmax,iidum
read(10,*)
read(10,*) acko
read(10,*)
read(10,*) kappa
read(10,*)
read(10,*) nfc   !number of frequency bands

allocate(fc1(nfc),fc2(nfc))
do i=1,nfc
 read(10,*) fc1(i),fc2(i)
enddo

close(10)

open(10,file='stainfo.dat',action='read')
allocate(fcsta(nr))
do i=1,nr
  read(10,*)ddum,ddum,ddum,ddum,ddum,ddum,fcsta(i)
enddo
close(10)


!Calculates JB distance (needed for distance dependent cross-over frequency)
allocate(staN(nr),staE(nr),fault(2,4),JBdist(nr))
open(10,FILE='fault.dat')
do i=1,4
  read(10,*)fault(:,i)
enddo
!fault=fault/1.d3
close(10)
open(10,FILE='stations.dat')
do i=1,nr
  read(10,*)staN(i),staE(i)
  JBdist(i)=JBdistance(fault,staN(i),staE(i))
enddo
close(10)
open(10,FILE='stations.JBdist.dat')
do i=1,nr
  write(10,*)staN(i),staE(i),JBdist(i)
enddo
close(10)

allocate(x1a(ng1),x2a(ng2))
allocate(cir(ng1,ng2,3,np),iir(ng1,ng2,3,np),rir(ng1,ng2,3,np),fltr4(np))
allocate(cseis(3,np),csseis(3,np))
allocate(smogramy(pmax,nr,3,np))
allocate(crossking(nfmax),crosskong(nfmax))
allocate(stf(nfmax),stftime(np)) 

smogramy=0.

dt=T/float(np)
df=1./T
mu=rho*(beta**2)
lambda=rho*(alfa**2)-2*mu
ms=Momain/(mu*leng*widt)
Moir=Momain

if(nfmax2>0)then
  open(3,file='XYGreen.dat')
  do i=1,ng1
    read(3,*) x1a(i)
  enddo
  do i=1,ng2
    read(3,*) x2a(i)
  enddo
  close(3)
endif

if(nfmax2>0)open(2,form='binary',file='NEZsor.dat')
open(20,form='binary',file='NEZsor2.dat')
if(pmax==1)then 
  idum=iidum
  call KKgener(leng,widt,hypo(2),(widt-hypo(1)),pmax,ms,vr,beta,acko,dip,depth,lmax,wmax,SUB,SUBmin,pulse,idum)
endif

do rec=1,nr
  write(*,*) 'Stanice c.',rec
  if(nfmax2>0)then
    cir=0.
    read(2) dum
    do i=1,ng1
      do j=1,ng2
        read(2) dum
        do k=1,nfmax2
          read(2) dumarr
          cir(i,j,1,k)=cmplx(dumarr(1),dumarr(4))
          cir(i,j,2,k)=cmplx(dumarr(2),dumarr(5))
          cir(i,j,3,k)=cmplx(dumarr(3),dumarr(6))
        enddo
!        do m=1,3
!          call cosfilters(cir(i,j,m,1:np),np,T,fc1,fc2,fc3,fc4)
!        enddo
      enddo
    enddo
    rir=real(cir)
    iir=imag(cir)
    write(*,*) 'Green. fce nacteny ...'
  endif
  idum=iidum
  
  do p=1,pmax
    KingRiseTime=leng/real(Pulse)/vr/2.
    allocate(slip(wmax,lmax),ruptime(wmax,lmax))
    ScenNum=p

!    ScenNum=9;write(*,*)'Fixed source model',ScenNum

    if(pmax>1)call KKgener(leng,widt,hypo(2),(widt-hypo(1)),ScenNum,ms,vr,beta,acko,dip,depth,lmax,wmax,SUB,SUBmin,pulse,idum)

!    pause 'Zkopiruj si slipgen.dat'
    
    open(1,file='KingKong.dat')
	read(1,*)
    do i=1,lmax
      do j=1,wmax
        read(1,*) dum2,dum2,slip(j,i),ruptime(j,i),fco
      enddo
    enddo

!Crossover function
     FFTcross1=int(.3/df)/2 !int(1./fco/df)/2
     FFTcross2=int(.3/df)*2 !int(1./fco/df)*2
!   if(JBdist(rec)<10.d0)then
!     FFTcross1=int(.6/df)/2 !int(.3/df)/2
!     FFTcross2=int(.6/df)*2 !int(.3/df)*2
!   else
!     FFTcross1=int(.6/df)/2 !int(.3/df)/2
!     FFTcross2=int(.6/df)*2 !int(.3/df)*2
!   endif
    if(nfmax2==0)then
      FFTcross1=1;FFTcross2=1;write(*,*)'! WARNING! Purely composite model!'
    endif
    write(*,*)'Crossover od',df*real(FFTcross1),' do',df*real(FFTcross2)
    do k=1,nfmax
      if(k<FFTcross1)then
        crossking(k)=1.d0;crosskong(k)=0.d0
      elseif(k<FFTcross2)then
        crossking(k)=cos(real(k-FFTcross1)*PI/2./real(FFTcross2-FFTcross1))**2
        crosskong(k)=sin(real(k-FFTcross1)*PI/2./real(FFTcross2-FFTcross1))**2
      else
        crossking(k)=0.d0;crosskong(k)=1.d0
      endif
    enddo
!    crossking=1.d0;crosskong=0.d0;write(*,*)'! WARNING! Purely integral model!'
        
    slip=slip/(widt*leng)*float(wmax*lmax)
    write(*,*)'Mean slip',sum(slip)/float(wmax*lmax)

    read(1,*)
    read(1,*) submax
    allocate(Mo(submax))
    allocate(fc(submax))
    write(*,*)submax,np
    allocate(spbrune(np))
    allocate(tr(submax))
    allocate(rcsxy(submax,2))
    allocate(cir2(submax,3,np))
    write(*,*)"Alokace OK" 

    cir2=0.
    read(20) dum
    do i=1,submax
      read(20) dum
      do k=1,nfmax
        read(20) dumarr
        cir2(i,1,k)=cmplx(dumarr(1),dumarr(4))
        cir2(i,2,k)=cmplx(dumarr(2),dumarr(5))
        cir2(i,3,k)=cmplx(dumarr(3),dumarr(6))
      enddo
!      do m=1,3
!        call cosfilters(cir2(i,m,1:np),np,T,fc1,fc2,fc3,fc4)
!      enddo
    enddo

    fco=1./fco
    call Brune(spbrune(1:np),np,fco,T)

    do i=1,submax
      read(1,*) rcsxy(i,2),rcsxy(i,1),Mo(i),tr(i),fc(i)
      rcsxy(i,1)=widt-rcsxy(i,1)
      fc(i)=1./fc(i)
    enddo
    close(105)
    close(1)

    if (stanice(rec)) goto 2000

    write(*,'(A13,I3,A2,A10,I3,A2,A9,I1)') 'realizace c.',p,',','stanice c.',rec
    cseis=0.

    elem=mu*leng*widt/float(wmax*lmax)
    dw=widt/(float(wmax))
    dl=leng/(float(lmax))

!$OMP parallel  private(m,n,li,wi,i,ww,ll,tr2,slippy,freq,dyry2a,dyy2a,ry2a,y2a,grre,grim,srcfce) DEFAULT(SHARED)
    allocate(dyry2a(max(ng1,ng2)),dyy2a(max(ng1,ng2)))
    allocate(y2a(ng2,wmax),ry2a(ng2,wmax))
    allocate(srcfce(submax))
!$OMP do SCHEDULE(DYNAMIC,1)
    do m=1,3
      write(*,*) 'King+Kong'
      do n=1,nfmax2
        freq=df*float((n-1))
        do li=1,ng2
          call spline(x1a(1:ng1),rir(1:ng1,li,m,n),ng1,1.e30,1.e30,dyry2a(1:ng1))
          call spline(x1a(1:ng1),iir(1:ng1,li,m,n),ng1,1.e30,1.e30,dyy2a(1:ng1))
	  do wi=1,wmax
	    ww=float(wi-1)*dw+dw/2.
	    call splint(x1a(1:ng1),rir(1:ng1,li,m,n),dyry2a(1:ng1),ng1,ww,ry2a(li,wi))
            call splint(x1a(1:ng1),iir(1:ng1,li,m,n),dyy2a(1:ng1),ng1,ww,y2a(li,wi))
          enddo
	enddo
	do wi=1,wmax
          call spline(x2a(1:ng2),ry2a(1:ng2,wi),ng2,1.e30,1.e30,dyry2a(1:ng2))
          call spline(x2a(1:ng2),y2a(1:ng2,wi),ng2,1.e30,1.e30,dyy2a(1:ng2))
          ww=float(wi-1)*dw+dw/2.
          do li=1,lmax
            ll=float(li-1)*dl+dl/2.
	    tr2=ruptime(wmax-wi+1,li);slippy=slip(wmax-wi+1,li)
            call splint(x2a(1:ng2),ry2a(1:ng2,wi),dyry2a(1:ng2),ng2,ll,grre)
            call splint(x2a(1:ng2),y2a(1:ng2,wi),dyy2a(1:ng2),ng2,ll,grim)
            cseis(m,n)=cseis(m,n)+(cmplx(grre,grim)*exp(-unit*twopi*tr2*freq))*elem*spbrune(n)*slippy*crossking(n)
          enddo
        enddo
!        cseis(m,n)=exp((-twopi*freq*kappa/2.))*cseis(m,n)
        if(m==3)then
          cseis(m,n)=exp((-twopi*freq*kappa/2./2.))*cseis(m,n)
        else
          cseis(m,n)=exp((-twopi*freq*kappa/2.))*cseis(m,n)
        endif
      enddo

      do n=FFTcross1,nfmax
        freq=df*float((n-1))
        call Brune2(srcfce(1:submax),submax,n,fc(1:submax),T)
        do i=1,submax
          srcfce(i)=srcfce(i)*mu*Mo(i)
          cseis(m,n)=cseis(m,n)+cir2(i,m,n)*srcfce(i)*exp(-unit*twopi*tr(i)*freq)*crosskong(n)
        enddo
!        cseis(m,n)=exp((-twopi*df*float((n-1))*kappa/2.))*cseis(m,n)
        if(m==3)then
          cseis(m,n)=exp((-twopi*freq*kappa/2./2.))*cseis(m,n)
        else
          cseis(m,n)=exp((-twopi*freq*kappa/2.))*cseis(m,n)
        endif
      enddo
      do i=np/2+2,np
        cseis(m,i)=conjg(cseis(m,np+2-i))
      enddo

      cseis(m,:)=cseis(m,:)*df
      call four1(cseis(m,1:np),np,1)
    enddo
!$omp end do

    deallocate(dyry2a,dyy2a,y2a,ry2a,srcfce)
!$omp end parallel

!Calculate source time function and other parameters
    if(rec==1)then
      allocate(srcfce(submax))
      stf(:)=0.
      open(863,FILE='stf.dat')
      do n=1,nfmax
        freq=df*float(n-1)
        call Brune2(srcfce(1:submax),submax,n,fc(1:submax),T)
        do i=1,submax    
          stf(n)=stf(n)+srcfce(i)*mu*Mo(i)*exp(-unit*twopi*tr(i)*freq)
        enddo
        write(863,*)freq,abs(stf(n))
      enddo
      deallocate(srcfce)
      close(863)
      stftime=0.
      open(863,FILE='stftime.dat')
      do n=1,np
        tr2=dt*(n-1)
        do i=1,submax    
          if(tr2>tr(i))then
            stftime(n)=stftime(n)+mu*Mo(i)*(2.*pi*fc(i))**2*(tr2-tr(i))*exp(-2.*pi*(tr2-tr(i))*fc(i))
          endif
        enddo
        write(863,*)tr2,stftime(n)
      enddo
      close(863)
      open(863,FILE='stfparameters.dat')
      fcmain=sqrt(mu/Momain*sqrt(sum((Mo(:)*twopi**2*fc(:)**2)**2)))/twopi
      stressdropmain=Momain*1.e7*(fcmain/4.9e6/beta*1.e3)**3
      write(863,*)'# Corner frequency of mainshock (Hz): ',fcmain
      write(863,*)'# Stress drop of mainshock (bar): ',stressdropmain
      write(863,*)'# Brune radius of mainshock (km): ',2.34*beta/twopi/fcmain/1000.
!      allocate(srcfce(1))
!      do n=1,nfmax
!        freq=df*float(n-1)
!        call Brune2(srcfce(1),1,n,fcmain,T)
!        write(863,*)freq,abs(srcfce(1))*Momain
!      enddo
!      deallocate(srcfce)
      write(863,*)'# -----'
      write(863,*)'# Seismic moment of largest subsource (Nm): ',mu*Mo(1)
      write(863,*)'# Corner frequency of largest subsource (Hz): ',fc(1)
      write(863,*)'# Stress drop of largest subsource (bar): ',mu*Mo(1)*1.e7*(fc(1)/4.9e6/beta*1.e3)**3
      write(863,*)'# Brune radius of largest subsource (km): ',2.34*beta/twopi/fc(1)/1000.
      close(863)
    endif
    
    do i=1,3
      fltr4(1:np) = cseis(i,1:np)
      if(fc1(fcsta(rec))>0.)then
        CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'BP', fc1(fcsta(rec)), fc2(fcsta(rec)), dt, 1, np)
      else
        CALL XAPIIR(fltr4, np, 'BU', 0.0, 0.0, 4,'LP', fc1(fcsta(rec)), fc2(fcsta(rec)), dt, 1, np)
      endif
      cseis(i,1:np) = fltr4(1:np)
    enddo

    smogramy(p,rec,1,1:np)=real(cseis(1,1:np))*cos(rotateto)+real(cseis(2,1:np))*sin(rotateto)
    smogramy(p,rec,2,1:np)=real(cseis(2,1:np))*cos(rotateto)-real(cseis(1,1:np))*sin(rotateto)
    smogramy(p,rec,3,1:np)=real(cseis(3,1:np))


2000 if (stanice(rec)) write(*,*) 'Skipped'

    deallocate(tr)
    deallocate(Mo)
    deallocate(fc)
    deallocate(spbrune)
    deallocate(rcsxy)
    deallocate(slip,ruptime)
    deallocate(cir2)
  enddo
  write(*,*) idum
enddo
close(2)

open(101,file='vseisn.dat')
open(103,file='vseise.dat')
open(105,file='vseisz.dat')
do p=1,pmax
  do i=1,np
    write(101,100) dt*(i-1),smogramy(p,1:nr,1,i)
    write(103,100) dt*(i-1),smogramy(p,1:nr,2,i)
    write(105,100) dt*(i-1),smogramy(p,1:nr,3,i)
  enddo
enddo

100 format(1000E13.5)

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Brune(U,np,fc,T)
implicit none
integer np
complex U(np), unit
real pi
integer Nyq,i
real df, dt, fc, T

pi=3.141592
unit=cmplx(0.,1.)
dt=T/float(np)
df=1./T
Nyq=np/2+1

U=0.

do i=1,Nyq
U(i)=1./(1.+unit*(i-1)*df/fc)**2
enddo

do i=Nyq+1,np
  U(i)=conjg(U(np+2-i))
  enddo
  
end subroutine
  
  

subroutine Brune2(U,submax,nf,fc,T)
implicit none
integer submax,nf,i
complex U(submax), unit
real pi
real df, dt, fc(submax), T
real triangl

pi=3.141592
unit=cmplx(0.,1.)
dt=T/float(nf)
df=1./T
U=0.
do i=1,submax
U(i)=1./(1.+unit*(nf-1)*df/fc(i))**2   !Brune's pulse
!Triangular source time functions
!triangl=2.
!if(nf==1)then
!  U(i)=1.
!else
!  U(i)=sin(pi*(nf-1)*df/2./(triangl*fc(i)))**2/(pi*(nf-1)*df/2./(triangl*fc(i)))**2*exp(-unit*(nf-1)*df/(triangl*fc(i)))
!endif

enddo

end subroutine

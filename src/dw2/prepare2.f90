program PREPARE2

implicit none

integer, parameter :: pulse=4,SUBmin=2,SUB=8 !stejne jako v Kingkong2.f90
real, parameter :: t0=20.
integer, parameter :: ngmax=100000
real, parameter :: pi=3.141592654
real df,aw1,var
complex rseis,ui,freq
integer np,ng1,ng2,totsub,pmax,idum,idum2,ScenNum
real gleng,gwidt
real TM(3,3),ITM(3,3)
real NEZhypo(3),xi(3),sour(3)
real Gx10,Gx20,dx1,dx2,x1a(ngmax),x2a(ngmax)
integer i,j,k
interface
  function Transf(NEZ,smer)
    logical smer
    real Transf(3)
    real NEZ(3)
  end function Transf
end interface
integer nc,nfreq,nr,ns,ikmax
real Stat(2)
real*8 strike,dip,rake
real*8 leng,widt
real*8 hhypo,hypo(3)
real*8 ms,vr,alfa,rho,beta,acko,dum,fstrike,fdip,frake,ran2k
real tl,aw,xl,uconv,fref,mu,Momain
namelist  /input/ nc,nfreq,tl,aw,nr,ns,xl,ikmax,uconv,fref
common /transform/ TM,NEZhypo,ITM,hypo
CHARACTER*6 filename

idum2=-90;var=120.d0
!idum2=-90;var=30.d0

write(*,*)'Variations: ',var

open(1,file='input.dat')

read(1,*)
read(1,*) nfreq,nfreq
read(1,*)
read(1,*) tl
read(1,*)
read(1,*) nr
read(1,*)
read(1,*) ng2,ng1,ng2,ng1
read(1,*)
read(1,*) gleng, gwidt
read(1,*)
read(1,*) Momain
read(1,*)
read(1,*) strike,dip,rake
read(1,*)
read(1,*) hhypo
read(1,*)
read(1,*) leng,widt
read(1,*)
read(1,*) hypo(2),hypo(1)
read(1,*)
read(1,*) np
read(1,*)
read(1,*) alfa,beta,rho
read(1,*)
read(1,*) vr
read(1,*)
read(1,*) pmax,idum
read(1,*) 
read(1,*) acko

close(1)

if(ng2>ngmax.or.ng1>ngmax)stop 'Check dimensions!'

hypo(3)=0.

open(2,file='GRDAT.HED')
aw=.5;ns=1;xl=3000000.;ikmax=200000;uconv=1.E-12;fref=1. !Axitra values that does not have to be generally changed
nc=0;  !set up formal value that are actualy not used by Axitra (they are readed from elsewhere)
write(2,input)
close(2)

open(4,file='station.dat')
open(5,file='stations.dat')
write(4,*) 'Station co-ordinates'
write(4,*) 'x(N>0,km),y(E>0,km),z(km),azim.,dist.,stat.'
do i=1,nr
  read(5,*) Stat(:) !,rec
  write(4,20) Stat(:),0.
enddo
write(4,'(A1)') char(26)
20 format(2x,f12.4,2x,f12.4,2x,f5.4,2x,f9.4,2x,f9.4,2x,i3)
close(4)
close(5)
    
NEZhypo(1)=0.
NEZhypo(2)=0.
NEZhypo(3)=-hhypo

TM(1,1)=sind(strike)*cosd(dip)
TM(1,2)=-cosd(strike)*cosd(dip)
TM(1,3)=sind(dip)
TM(2,1)=cosd(strike)
TM(2,2)=sind(strike)
TM(2,3)=0.
TM(3,1)=-sind(strike)*sind(dip)
TM(3,2)=cosd(strike)*sind(dip)
TM(3,3)=cosd(dip)
ITM=transpose(TM)

open(3,file='XYGreen2.dat')
open(1,file='sources2.dat')
k=0
mu=rho*(beta**2)
ms=Momain/(mu*leng*widt)
do i=1,pmax
  ScenNum=i
  
!  ScenNum=1;write(*,*)'Fixed source model',ScenNum
  write(*,*)'IDUM: ',idum
  call KKgener(leng,widt,hypo(2),(widt-hypo(1)),ScenNum,ms,vr,beta,acko,dip,hhypo,ng2,ng1,SUB,SUBmin,pulse,idum)
      
!  pause 'Uloz si slipgen.dat'
      
  write(*,*)'IDUM: ',idum
  open(343,FILE='KingKong.dat')
  read(343,*)
  do j=1,ng1*ng2
    read(343,*)
  enddo
  read(343,*)
  read(343,*)totsub
  do j=1,totsub
    k=k+1
    if(k<10)then
      write(filename,'(A5,I1)')'00000',k
    elseif(k<100)then
      write(filename,'(A4,I2)')'0000',k
    elseif(k<1000)then
      write(filename,'(A3,I3)')'000',k
    elseif(k<10000)then
      write(filename,'(A2,I4)')'00',k
    elseif(k<100000)then
      write(filename,'(A1,I5)')'0',k
    else
      write(filename,'(I6)')k
    endif
    read(343,*) Gx20,Gx10
    Gx10=widt-Gx10
    write(3,*) Gx10,Gx20
    xi(1)=Gx10
    xi(2)=Gx20
    xi(3)=0.
    if(j>3)then
      fstrike=strike+(ran2k(idum2)-0.5d0)*var
      fdip=dip+(ran2k(idum2)-0.5d0)*var
      frake=rake+(ran2k(idum2)-0.5d0)*var
    else
      fstrike=strike
      fdip=dip
      frake=rake
    endif
    sour=Transf(xi,.FALSE.)
    sour(3)=-sour(3)
    sour=sour/1000.    
    write(1,'(A6,6E13.5)') filename,sour,fstrike,fdip,frake
  enddo
  close(343)
enddo
                                                      
close(1)
open(1,file='fault.dat')
xi=0.
write(1,*) Transf(xi,.FALSE.)/1000.
xi(1)=widt
xi(2)=0.
xi(3)=0.
write(1,*) Transf(xi,.FALSE.)/1000.
xi(1)=widt
xi(2)=leng
xi(3)=0.
write(1,*) Transf(xi,.FALSE.)/1000.
xi(1)=0.
xi(2)=leng
xi(3)=0.
write(1,*) Transf(xi,.FALSE.)/1000.
xi=0.
write(1,*) Transf(xi,.FALSE.)/1000.
close(1)    

open(1,file='dirac.dat')
ui=cmplx(0.,1.)
df=1./tl
aw1=-aw/(2.*tl)
do i=1,nfreq
  freq=cmplx(df*(i-1),aw1)
  rseis=exp(-ui*2.*pi*t0*freq)
  write(1,*) real(rseis),imag(rseis)
enddo
do i=nfreq+1,np
  rseis=0.d0
  write(1,*) real(rseis),imag(rseis)
enddo
close(1)        

end

function Transf(NEZ,smer)
implicit none
logical smer
real Transf(3)
real NEZhypo(3)
real TM(3,3),ITM(3,3)
real*8 hypo(3)
real NEZ(3)
common /transform/ TM,NEZhypo,ITM,hypo

if (smer) then
  Transf=matmul(TM,(NEZ-NEZhypo))+hypo
else
  Transf=matmul(ITM,(NEZ-hypo))+NEZhypo
endif
end function


      FUNCTION ran2k(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2k,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2k=min(AM*iy,RNMX)
      return
      END

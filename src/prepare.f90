program PREPARE

implicit none

real, parameter :: t0=20.
integer, parameter :: ngmax=10000
real, parameter :: pi=3.141592654
real df,aw1
complex rseis,ui,freq
integer ng1,ng2,np
real strike,dip,rake
real leng,widt
real hhypo
real gleng,gwidt
real TM(3,3),ITM(3,3)
real NEZhypo(3),hypo(3),xi(3),sour(3)
real Gx10,Gx20,dx1,dx2,x1a(ngmax),x2a(ngmax)
integer i,j,k
interface
  function Transf(NEZ,smer)
    logical smer
    real Transf(3)
    real NEZ(3)
  end function Transf
end interface
common /transform/ TM,NEZhypo,ITM,hypo
integer nc,nfreq,nr,ns,ikmax
real Stat(2)
real tl,aw,xl,uconv,fref
namelist  /input/ nc,nfreq,tl,aw,nr,ns,xl,ikmax,uconv,fref
CHARACTER*6 filename

open(1,file='input.dat')

read(1,*)
read(1,*) nfreq
read(1,*)
read(1,*) tl
read(1,*)
read(1,*) nr
read(1,*)
read(1,*) ng2,ng1
read(1,*)
read(1,*) gleng,gwidt
read(1,*)
read(1,*)
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
close(1)

if(ng2>ngmax.or.ng1>ngmax)stop 'Check dimensions!'

hypo(3)=0.

open(2,file='GRDAT.HED')
aw=1.;ns=1;xl=5.747506e+06;ikmax=100000;uconv=0.1E-03;fref=1. !Axitra values that does not have to be generally changed
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

! Change 9.12.2002

Gx10=(widt-gwidt)/2.
Gx20=(leng-gleng)/2.
dx1=gwidt/float((ng1-1))
dx2=gleng/float((ng2-1))

open(3,file='XYGreen.dat')
do i=1,ng1
  x1a(i)=Gx10+float(i-1)*dx1
  write(3,*) x1a(i)
enddo
do i=1,ng2
  x2a(i)=Gx20+float(i-1)*dx2
  write(3,*) x2a(i)
enddo

open(1,file='sources.dat')
k=0
do i=1,ng1
  do j=1,ng2
    xi(1)=x1a(i)
    xi(2)=x2a(j)
    xi(3)=0.
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
    sour=Transf(xi,.FALSE.)
    sour(3)=-sour(3)
    sour=sour/1000.
    write(1,'(A6,6E13.5)') filename,sour,strike,dip,rake
  enddo
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

ui=cmplx(0.,1.)
df=1./tl
aw1=-aw/(2.*tl)
open(1,file='dirac.dat')
do i=1,nfreq
  freq=cmplx(df*(i-1),aw1)
  rseis=exp(-ui*2.*pi*t0*freq)
  write(1,*) real(rseis),imag(rseis)
enddo
do i=nfreq+1,np
  rseis=0.
  write(1,*) real(rseis),imag(rseis)
enddo
close(1)        

end

function Transf(NEZ,smer)
implicit none
logical smer
real Transf(3)
real NEZhypo(3)
real TM(3,3),ITM(3,3),hypo(3)
real NEZ(3)
common /transform/ TM,NEZhypo,ITM,hypo

if (smer) then
  Transf=matmul(TM,(NEZ-NEZhypo))+hypo
else
  Transf=matmul(ITM,(NEZ-hypo))+NEZhypo
endif
end function

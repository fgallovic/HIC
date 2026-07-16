program vytvor

implicit none

integer, parameter :: NF=74, NR=18 

integer :: i,j,statweightz(NR),statweighte(NR),statweightn(NR),Z,E,N
real*8 :: rzobs(NF,NR+1),reobs(NF,NR+1),rnobs(NF,NR+1),rzsyn(NF,NR+1),resyn(NF,NR+1),rnsyn(NF,NR+1)
real*8 :: rzkomp(NF,NR+3),rekomp(NF,NR+3),rnkomp(NF,NR+3),avgofavgz,avgofavge,avgofavgn,avgall
real*8 :: avgz,avge,avgn
character, dimension(4) :: statnamez(NR),statnamee(NR),statnamen(NR)

open(101,file='observed/rapseisz.dat')
open(102,file='observed/rapseise.dat')
open(103,file='observed/rapseisn.dat')
open(104,file='rapseisz.dat')
open(105,file='rapseise.dat')
open(106,file='rapseisn.dat')
open(107,file='testZ.dat')
open(108,file='testE.dat')
open(109,file='testN.dat')

rzkomp=0.
rekomp=0.
rnkomp=0.
Z=0
E=0
N=0

do i=1,NR
  read(107,*)statweightz(i),statnamez(i)
  if (statweightz(i)==1) then
    Z=Z+1
  endif
  read(108,*)statweighte(i),statnamee(i)
  if (statweighte(i)==1) then
    E=E+1
  endif
  read(109,*)statweightn(i),statnamen(i)
  if (statweightn(i)==1) then
    N=N+1
  endif
enddo

close(107);close(108);close(109)

do i=1,NF
  read(101,*)rzobs(i,:)
  read(102,*)reobs(i,:)
  read(103,*)rnobs(i,:)
  read(104,*)rzsyn(i,:)
  read(105,*)resyn(i,:)
  read(106,*)rnsyn(i,:)
  do j=2,NR+1 
    rzkomp(i,j)=log(rzsyn(i,j)/rzobs(i,j))*statweightz(j-1)
    rekomp(i,j)=log(resyn(i,j)/reobs(i,j))*statweighte(j-1)
    rnkomp(i,j)=log(rnsyn(i,j)/rnobs(i,j))*statweightn(j-1)
  enddo
enddo

close(101);close(102);close(103);close(104);close(105);close(106)

rzkomp(:,1)=rzobs(:,1)
rekomp(:,1)=reobs(:,1)
rnkomp(:,1)=rnobs(:,1)

avgz=0
avge=0
avgn=0

do i=1,NF
  avgz=avgz+sum((rzkomp(i,2:NR+1)))/real(Z)
  avge=avge+sum((rekomp(i,2:NR+1)))/real(E)
  avgn=avgn+sum((rnkomp(i,2:NR+1)))/real(N)
  !avgz=avge+sum(rzkomp(i,2:NR+1)**2)/real(Z)
  !avge=avge+sum(rekomp(i,2:NR+1)**2)/real(E)
  !avgn=avgn+sum(rnkomp(i,2:NR+1)**2)/real(N)
  rzkomp(i,NR+2)=sum(rzkomp(i,2:NR+1))/real(Z)
  rekomp(i,NR+2)=sum(rekomp(i,2:NR+1))/real(E)
  rnkomp(i,NR+2)=sum(rnkomp(i,2:NR+1))/real(N)
enddo

avgofavgz=avgz/real(NF)
avgofavge=avge/real(NF)
avgofavgn=avgn/real(NF)


avgall=(avgofavgz+avgofavgn+avgofavge)/3.
open(400,file='avgofallfreq.dat')
write(400,'(23F9.5)')avgofavgz, avgofavge,avgofavgn, avgall

!do i=1,NF
!  rzkomp(i,NR+3)=sum((rzkomp(i,NR+2)-avgofavgz)**2)

do i=1,NF
  rzkomp(i,NR+3)=sqrt(sum((rzkomp(i,2:NR+1)-rzkomp(i,NR+2)*statweightz(1:NR))**2)/real(Z))
  rekomp(i,NR+3)=sqrt(sum((rekomp(i,2:NR+1)-rekomp(i,NR+2)*statweighte(1:NR))**2)/real(E))
  rnkomp(i,NR+3)=sqrt(sum((rnkomp(i,2:NR+1)-rnkomp(i,NR+2)*statweightn(1:NR))**2)/real(N))
enddo

open(201,file='rzkomp.dat')
open(202,file='rekomp.dat')
open(203,file='rnkomp.dat')

do i=1,NF
    write(201,'(1000F9.5)')rzkomp(i,:)
    write(202,'(1000F9.5)')rekomp(i,:)
    write(203,'(1000F9.5)')rnkomp(i,:)
enddo

close(201);close(202);close(203)

end program vytvor

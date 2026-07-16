    SUBROUTINE KKgener(LL,WW,epicL,epicW,ScenNum,ms,vr,vs,a,dip,depth,NL,NW,SUB,SUBmin,Puls,idum)
    IMPLICIT NONE
    REAL,PARAMETER:: gamma=1,beta=2
    REAL*8,PARAMETER::vs2=3100.d0,top=2000.d0
    REAL*8,PARAMETER:: PI=3.1415926535
    REAL*8 LL,WW,epicL,epicW,ms,vr,Vs,ran2,dum,dip,depth,h1
    INTEGER NL,NW,SUB,SUBmin,Puls,idum,pocet,asperL,asperW
    REAL*8 dumL,dumW,sizeL,sizeW,RTime,ms1,slipnorm,KingRiseTime,a,compRT
    INTEGER i,j,odL,doL,odW,doW,stL,stW,k,l,m,ii,jj,n,NLdum,NWdum
    REAL*8,ALLOCATABLE:: m0(:,:),kompM0(:),kompRT(:),posL(:),posW(:)
    REAL, ALLOCATABLE:: AA(:,:),AB(:,:)
    REAL slipvalue,krokL,krokW
    INTEGER idum2,ScenNum
!PDF for subsource position
    INTEGER,PARAMETER:: pdfNL=20000,pdfNW=1000
    INTEGER fileNL,fileNW,pdfOption,ml(2)
    REAL*8,ALLOCATABLE:: pdf2D(:,:),cpdf2D(:,:)
    REAL*8 pdfDL,pdfDW,pdfGaussL,pdfGaussW,pdfGaussSL,pdfGaussSW
    CHARACTER filename*256

!Preparing PDF for distribution of subsources
    pdfGaussL=LL/2.;pdfGaussW=WW/2.;pdfGaussSL=LL/4.;pdfGaussSW=WW/4.
    fileNL=134;fileNW=101;filename='GreeceDeep-slip2D-r.dat'
    pdfOption=1 !PDF for subsource distribution (1: uniform, 2: Gaussian - provide location and width, 3: read from file - specify discretization and filename)
    write(*,*)'Preparing PDF for subsource distribution...'
    ALLOCATE(pdf2D(pdfNL,pdfNW),cpdf2D(pdfNL,pdfNW))
    CALL fillpdf(pdfNL,pdfNW,LL,WW,pdf2D,cpdf2D,pdfOption,filename,fileNL,fileNW,pdfGaussL,pdfGaussW,pdfGaussSL,pdfGaussSW)
    pdfDL=LL/real(pdfNL)
    pdfDW=WW/real(pdfNW)
    write(*,*)'... done.'
    
    idum2=-10
    h1=(depth-top)/sin(dip/180.d0*3.1415926d0)+WW-epicW
    write(*,*)'h1=',h1
    OPEN(485,FILE='KingKong.dat')
    OPEN(486,FILE='slipgen.dat')
    ms1=ms

    pocet=(SUB**2-(SUBmin-1)**2)*int(LL/WW)

    ALLOCATE(m0(NL,NW),kompM0(pocet),kompRT(pocet),posL(pocet),posW(pocet))
    m0=0.

    write(485,*)'KING', NL, NW

    KingRiseTime=LL/real(Puls)/vr
    KingRiseTime=.1d0;write(*,*)'Rise time fixed at 0.1sec'

    m=0;i=(SUBmin-1)**2*int(LL/WW)
    do n=SUBmin,SUB !loop over sizes of the subsources
      write(*,*) n,i,(n**2)*int(LL/WW)-i
      sizeL=WW/real(n)
      sizeW=WW/real(n)
      do j=1,(n**2)*int(LL/WW)-i       !loop over subsources of the same size
        m=m+1
        kompM0(m)=ms1/real(n)**gamma*sizeL*sizeW
        kompRT(m)=sqrt(sizeL*sizeW/PI)/0.37/Vs    !Brune
        kompRT(m)=sqrt(sizeL*sizeW)/vr/a

!Definition of subsource centers:

!        if(n==2)then
        if(n==1)then
            
!Asperity in the center
          if(ScenNum==1.or.ScenNum==2)then
            if(m<=2)then
              dumL=LL/2.;dumW=WW/2.
            else

71            dumL=ran2(idum);dumL=dumL*(LL-sizeL)+sizeL/2.
              dumW=ran2(idum);dumW=dumW*(WW-sizeW)+sizeW/2.
              !write(*,*)'dumL: ', dumL, 'dumW: ', dumW, 'IDUM: ',idum
              if((dumL>3./8.*LL.and.dumL<5./8.*LL).and.(dumW>3./8.*WW.and.dumW<5./8.*WW))goto 71
            endif

!Asperity at one side
          elseif(ScenNum==3.or.ScenNum==4)then
            if(m<=3)then
              dumL=3.*LL/4.;dumW=WW/4.*m
            else
72            dumL=ran2(idum);dumL=dumL*(LL-sizeL)+sizeL/2.
              dumW=ran2(idum);dumW=dumW*(WW-sizeW)+sizeW/2.
              if(dumL>1./2.*LL)goto 72
            endif

!Asperity at the other side
          elseif(ScenNum==5.or.ScenNum==6)then
            if(m<=3)then
              dumL=1.*LL/4.;dumW=WW/4.*m
            else
73            dumL=ran2(idum);dumL=dumL*(LL-sizeL)+sizeL/2.
              dumW=ran2(idum);dumW=dumW*(WW-sizeW)+sizeW/2.
              if(dumL<1./2.*LL)goto 73
            endif

!Asperity at both sides
          elseif(ScenNum==7)then
            if(m<=2)then
              dumL=1.*LL/4.;dumW=WW/3.*m
            elseif(m==3)then
              dumL=3.*LL/4.;dumW=WW/2.
            else
74            dumL=ran2(idum);dumL=dumL*(LL-sizeL)+sizeL/2.
              dumW=ran2(idum);dumW=dumW*(WW-sizeW)+sizeW/2.
              if(dumL<1./2.*LL)goto 74
            endif
    
!Asperity at both sides Shallower
          elseif(ScenNum==8)then
            if(m==1)then
              dumL=1.*LL/4.;dumW=WW/3.
            elseif(m==2)then
              dumL=1.*LL/2.;dumW=WW/2.
            elseif(m==3)then
              dumL=3.*LL/4.;dumW=WW/2.
            else
75            dumL=ran2(idum);dumL=dumL*(LL-sizeL)+sizeL/2.
              dumW=ran2(idum);dumW=dumW*(WW-sizeW)+sizeW/2.
              if(dumL<1./2.*LL)goto 75
            endif

!Model of L'Aquila from inversions
          elseif(ScenNum==9)then
            if(m==1)then
              dumL=8000.;dumW=15000.-9000.
            elseif(m==2)then
              dumL=12000.;dumW=15000.-4000.
            elseif(m==3)then
              dumL=15000.;dumW=15000.-4000.
            endif
          endif
  
        else
!          dumL=ran2(idum);dumL=dumL*(LL-sizeL)+sizeL/2.
!          dumW=ran2(idum);dumW=dumW*(WW-sizeW)+sizeW/2.
          do
            ml=minloc(abs(cpdf2D(:,:)-ran2(idum)))
            dumL=(real(ml(1))-.5)*pdfDL
            dumW=(real(ml(2))-.5)*pdfDW
            if(dumL-sizeL/2.>=0..and.dumW-sizeW/2.>=0..and.dumL+sizeL/2.<=LL.and.dumW+sizeW/2.<=WW)exit
          enddo
        endif
  
        posL(m)=dumL;posW(m)=dumW

!oprava na skok v rupture velocity
        if(dumW+sizeW/2.<WW-h1)then
          kompRT(m)=kompRT(m)*vs/vs2;kompM0(m)=kompM0(m)*vs2**2/vs**2
        elseif(dumW-sizeW/2.<WW-h1)then
          kompRT(m)=kompRT(m)*(1+(WW-h1-dumW+sizeW/2.)/sizeW*(vs/vs2-1))
          kompM0(m)=kompM0(m)*(1+(WW-h1-dumW+sizeW/2.)/sizeW*(vs2**2/vs**2-1))
        endif

        odL=int((dumL-sizeL/2.)/LL*real(NL))+1;doL=int((dumL+sizeL/2.)/LL*real(NL)+.999999);stL=real(doL+odL)/2.
        odW=int((dumW-sizeW/2.)/WW*real(NW))+1;doW=int((dumW+sizeW/2.)/WW*real(NW)+.999999);stW=real(doW+odW)/2.

        if((doL-odL)<3.or.(doW-odW)<3)then
          m0((odL+doL)/2,(odW+doW)/2)=m0((odL+doL)/2,(odW+doW)/2)+ms1*LL*WW/real(NL*NW)/real(n)**gamma           
        else
          NLdum=2**int(log(real(doL-odL+1))/log(2.d0)+1)
          NWdum=2**int(log(real(doW-odW+1))/log(2.d0)+1)
          ALLOCATE(AA(NLdum,NWdum),AB(doL-odL+1,doW-odW+1))
          CALL slip(sizeL,sizeW,NLdum,NWdum,1.d0,idum2,AA)
          krokL=real(NLdum-1)/real(doL-odL);krokW=real(NWdum-1)/real(doW-odW)
          do k=1,doL-odL+1
            ii=real(k-1)*krokL+1
            do l=1,doW-odW+1
              jj=real(l-1)*krokW+1
              AB(k,l)=AA(ii,jj)
            enddo
          enddo
          dum=sum(AB)
          AB=AB/dum*real((doL-odL+1)*(doW-odW+1))
          ii=0
          do k=odL,doL
            ii=ii+1;jj=0
            do l=odW,doW
              jj=jj+1
              m0(k,l)=m0(k,l)+ms1*LL*WW/real(NL*NW)/real(n)**gamma*AB(ii,jj)
            enddo
          enddo
          DEALLOCATE(AA,AB)
        endif
      enddo
!      write(*,*)'Level n=',n,'//',SUB,' done.'
      i=n**2*int(LL/WW)
! for illustration figures (individual subsources):
!      open(375,file='slipdum.txt')
!      do i=1,NL
!        do j=1,NW
!          write(375,*)real(i-1)*LL/real(NL-1),real(j-1)*WW/real(NW-1),m0(i,j)/LL/WW*real(NL*NW)
!        enddo
!      enddo
!      close(375)
!      pause 'copy slipdum.txt somewhere...';m0=0.d0
      
    enddo

!Correction to the jump in the rupture velocity
    m0(1:NL,1:int((WW-h1)/WW*NW))=m0(1:NL,1:int((WW-h1)/WW*NW))*vs2**2/vs**2
!write(*,*)'Nekoriguji moment na zmenu mu u povrchu!'
    ms1=ms/(sum(m0)/LL/WW);m0=m0*ms1
    ms1=ms/(sum(KompM0)/LL/WW);KompM0=KompM0*ms1       !correcting the final mean slip
    write(*,*)'Mean slip:        ',sum(m0)/LL/WW
    write(*,*)'Slip on asperity: ',sum(m0(1:NL/2,NW/2:NW))*4/LL/WW
    ms1=sqrt(sum(kompM0(1:pocet)**2/kompRT(1:pocet)**4))/(vr**2*ms);write(*,*)'Adjusting a=',sqrt(ms1),'to a=',a;kompRT=kompRT*sqrt(ms1)/a

        
    do k=1,NL
      do l=1,NW
        RTime=compRT(abs(real(k,8)*LL/real(NL,8)-EpicL),real(l,8)*WW/real(NW,8),EpicW,vr,vs,vs2,WW-h1)

!Model of L'Aquila from inversions
        if(ScenNum==9)then
          if(real(k,8)*LL/real(NL,8)>9000..and.real(l,8)*WW/real(NW,8)>15000.-7500.)RTime=RTime+3.d0
        endif

        write(485,'(100E14.6)')real(k-1)*LL/real(NL-1),real(l-1)*WW/real(NW-1),m0(k,l),RTime,KingRiseTime
        write(486,'(100E14.6)')real(k-1)*LL/real(NL-1),real(l-1)*WW/real(NW-1),m0(k,l)*NL*NW/LL/WW,RTime
!        write(486,'(100E14.6)')real(k-1)*LL/real(NL-1),real(l-1)*WW/real(NW-1),m0(k,l),RTime
      enddo
      write(486,*)
    enddo
        

    write(485,*)'KONG'
    write(485,*)pocet

    do m=1,pocet
      RTime=compRT(abs(posL(m)-EpicL),posW(m),EpicW,vr,vs,vs2,WW-h1)

!Model of L'Aquila from inversions
      if(ScenNum==9)then
        if(posL(m)>9000..and.posW(m)>15000.-7500.)RTime=RTime+3.d0
      endif

      write(485,'(100E14.6)')posL(m),posW(m),KompM0(m),RTime,kompRT(m)
    enddo

    CLOSE(485)
    CLOSE(486)

    END

    FUNCTION slipvalue(x,y)
    REAL x,y,slipvalue
    REAL*8,PARAMETER:: PI=3.1415926535
    
    !slipvalue=cos(x*PI)*cos(y*PI)
    !if(.5<sqrt(x**2+y**2))then
    !  slipvalue=0.
    !else 
    !  slipvalue=.5-sqrt(x**2+y**2)
    !endif
    slipvalue=(sqrt(x**2+y**2)+1)*exp(-sqrt(x**2+y**2)*4.)
    slipvalue=1.
    if(.24<x**2+y**2)then
      slipvalue=0.
    else 
      slipvalue=sqrt(.24-x**2-y**2)
    endif

!    if(.24<x**2+y**2)then   !Savage(1972)
!      slipvalue=0.
!    else 
!      slipvalue=1-exp(-(.5-sqrt(x**2+y**2))/.5)
!    endif


    END

    
    SUBROUTINE slip(L,W,M,N,ms,idum,A)

    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535
    COMPLEX*16 wh
    INTEGER i,j,k,M,N,FNM,FNN,idum
    COMPLEX speq1(N),AC(M/2,N)
    REAL A(M,N)
    REAL*8 dkx,dky,L,W,kx,ky,dx,dy,fkx,fky,KCx,KCy,corner
    REAL*8 cx,cy,ms,mini,dum,gasdev,ran2s
    INTEGER ci,cj,times,z
    REAL*8 arx,ary,phx,phy

    FNN=N/2+1;FNM=M/2+1

!Generating spectrum
    corner=.8d0
    dkx=1./L;dky=1./W
    dx=1./dkx/real(M)
    dy=1./dky/real(N)
    fkx=PI*L;fky=PI*W
    KCx=corner/L
    KCy=corner/W
!    arx=dkx;ary=dky
    arx=KCx;ary=KCy
    arx=dkx
    ary=dky

      do i=1,M
        do j=1,N
!          a(i,j)=gasdev()
          a(i,j)=ran2s(idum)
        enddo
      enddo

      CALL rlft3(A,speq1,M,N,1,1)
!      mini=A(1,1);speq1=speq1/A(1,1);A=A/mini
      speq1=exp(cmplx(0.,atan2(imag(speq1),real(speq1))))
      do i=1,M/2
        AC(i,1:N)=exp(cmplx(0.,atan2(A(2*i,1:N),A(2*i-1,1:N))))
      enddo

      do j=1,N
        if(j<=N/2+1)then
          ky=dky*real(j-1)
        else
          ky=-dky*real(N-j+1)
        endif
        do i=1,M/2+1
          kx=dkx*real(i-1)
          if(kx**2+ky**2<=arx**2+ary**2)then
            AC(i,j)=abs(AC(i,j))*exp(cmplx(0.,(kx*fkx+ky*fky)))
          endif
          wh=ms*L*W/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**2)
          if(i<M/2+1)then
            AC(i,j)=AC(i,j)*wh
          else
            speq1(j)=speq1(j)*wh
          endif
        enddo
      enddo

      CALL rlft3(AC,speq1,M,N,1,-1)
      do i=1,M/2
        A(2*i-1,1:N)=real(AC(i,1:N))/M/N*2.
        A(2*i,1:N)=imag(AC(i,1:N))/M/N*2.
      enddo


!Modification in spatial domain:

      mini=minval(A)
!      write(*,*)'Total of samples:  ',M*N
!      write(*,*)'Minimum: ',mini
!      A=A-mini
      do i=1,M
        do j=1,N
          if(A(i,j)<0.d0) A(i,j)=0.d0
        enddo
      enddo

      
!Cutting edges in spatial domain 
      cy=W/8.;cx=L/8.
      ci=int(cx/dx);cj=int(cy/dy)
      do j=1,cj+1
        A(1:M,j)=A(1:M,j)*(.5+.5*cos(PI*real(cj-j+1)/real(cj+1)))
        k=N-cj+j-1
        A(1:M,k)=A(1:M,k)*(.5+.5*cos(PI*real(k-N+cj)/real(cj+1)))
      enddo
      do i=1,ci+1
        A(i,1:N)=A(i,1:N)*(.5+.5*cos(PI*real(ci-i+1)/real(ci+1)))
        k=M-ci+i-1
        A(k,1:N)=A(k,1:N)*(.5+.5*cos(PI*real(k-M+ci)/real(ci+1)))
      enddo

!imposing the mean slip ms
      dum=sum(A)
      A=A*ms/dum*real(M*N)
      

!Writing spatial domain
!    do i=1,M
!      do j=1,N
!        write(101,'(3E13.6)') real(i-1)*dx,real(j-1)*dy,a(i,j)
!        write(101,'(3E13.6)') real(i-1)*dx,real(j-1)*dy,1.d0
!      enddo
!    enddo

!    j=N/2
!    do i=1,M
!      write(120,*)real(i-1)*dx,a(i,j)
!    enddo
!    i=M/2
!    do j=1,N
!      write(121,*)real(j-1)*dy,a(i,j)
!    enddo


!    do i=1,M/2
!      AC(i,:)=cmplx(A(2*i-1,:),A(2*i,:))
!    enddo
!    CALL rlft3(AC,speq1,M,N,1,1)
!    AC=AC/real(M*N/2);speq1=speq1/real(M*N/2)

!Writing spec. along y:
!    do i=1,N/2+1
!      write(106,*)(i-1)*dky,abs(AC(1,i))
!    enddo
!Writes spec. along x:
!    do i=1,M/2
!      write(105,*)(i-1)*dkx,abs(AC(i,1))
!    enddo
!      write(105,*)(M/2)*dkx,abs(speq1(1))

    END


      FUNCTION gasdev()
      DOUBLE PRECISION gasdev
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       CALL RANDOM_NUMBER(ran1)
        v1=2.d0*ran1-1.d0
        CALL RANDOM_NUMBER(ran1)
        v2=2.d0*ran1-1.d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END


    FUNCTION compRT(x,y,hypoy,vr1,vs1,vs2,h2)
    IMPLICIT NONE
    REAL*8 h2,vs1,vs2,vr1,vr2
    REAL*8 x,y,compRT,rtbis,xx,func,hypoy
    EXTERNAL func
    if(hypoy<h2)then
      write(*,*)'Hypocentre lies in the low-velocity zone!';stop
    endif
    vr2=vr1/vs1*vs2
    if(y>h2)then
      compRT=sqrt(x**2+(y-hypoy)**2)/vr1
      return
    elseif(x<.0001d0)then
      xx=0.
    else
      xx=rtbis(func,0.d0,x,.0001d0,vr1,vr2,h2,y)
!      write(*,*)x,xx
    endif
    compRT=sqrt(xx**2+(h2-hypoy)**2)/vr1+sqrt((x-xx)**2+(y-h2)**2)/vr2
    END
    
    FUNCTION func(x,xx,vr1,vr2,h2,y)
    IMPLICIT NONE
    REAL*8 x,xx,func,h2,vr1,vr2,y
    func=vr2/vr1*x/sqrt((h2)**2+x**2)-(xx-x)/sqrt((y-h2)**2+(xx-x)**2)
    END
    
      
      FUNCTION rtbis(func,x1,x2,xacc,vr1,vr2,h2,y)
      IMPLICIT NONE
      INTEGER JMAX
      DOUBLE PRECISION rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      DOUBLE PRECISION dx,f,fmid,xmid,vr1,vr2,h2,y
      fmid=func(x2,x2,vr1,vr2,h2,y)
      f=func(x1,x2,vr1,vr2,h2,y)
      if(f*fmid.ge.0.d0) pause 'root must be bracketed in rtbis'
      if(f.lt.0.d0)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5d0
        xmid=rtbis+dx
        fmid=func(xmid,x2,vr1,vr2,h2,y)
        if(fmid.le.0.d0)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.d0) return
11    continue
      pause 'too many bisections in rtbis'
    END

    
    SUBROUTINE fillpdf(pdfNL,pdfNW,LF,WF,pdf2D,cpdf2D,pdfOption,filename,fileNL,fileNW,pdfGaussL,pdfGaussW,pdfGaussSL,pdfGaussSW)  ! creates pdf and cumulative pdf for subsource distribution
    IMPLICIT NONE
    INTEGER pdfNL,pdfNW,pdfOption
    REAL*8 pdf2D(pdfNL,pdfNW),cpdf2D(pdfNL,pdfNW)
    REAL*8 LF,WF,L,W,smL,smW,pdfGaussL,pdfGaussW,pdfGaussSL,pdfGaussSW
    CHARACTER*256 filename
    INTEGER i,j,k,fileNL,fileNW
    REAL*8 cumul,pdfDL,pdfDW,slipDL,slipDW
    REAL*8,ALLOCATABLE:: slip(:,:)
    INTEGER pifrom,pito,pjfrom,pjto
    smL=0.;smW=0.;L=LF;W=WF
    pdfDL=LF/real(pdfNL)
    pdfDW=WF/real(pdfNW)
    pifrom=int(smL/pdfDL)+1
    pito=int((smL+L)/pdfDL+0.999)
    pjfrom=int(smW/pdfDW)+1
    pjto=int((smW+W)/pdfDW+0.999)
    pdf2D(:,:)=0.
    SELECT CASE(pdfOption)
    CASE(1)
      write(*,*)'Uniform spatial PDF for subsources'
      pdf2D(pifrom:pito,pjfrom:pjto)=1.
    CASE(2)
      write(*,*)'Gaussian PDF for subsources'
      do j=pjfrom,pjto
        do i=pifrom,pito
          pdf2D(i,j)=exp(-.5*(((real(i)-.5)*pdfDL-pdfGaussL)**2/pdfGaussSL**2+((real(j)-.5)*pdfDW-pdfGaussW)**2/pdfGaussSW**2))
        enddo
      enddo
    CASE(3)
      write(*,*)'Reading spatial PDF from',trim(filename)
      allocate(slip(fileNL,fileNW))
      slipDL=LF/real(fileNL)
      slipDW=WF/real(fileNW)
      open(329,FILE=trim(filename))
      do j=1,fileNW
        read(329,*)(slip(i,j),i=1,fileNL)
      enddo
      close(329)
      do j=pjfrom,pjto
        do i=pifrom,pito
          pdf2D(i,j)=slip(int(pdfDL/slipDL*(float(i)-0.5))+1,int(pdfDW/slipDW*(float(j)-0.5))+1)
        enddo
      enddo
      pdf2D=pdf2D+.1*maxval(pdf2D)
      deallocate(slip)
    CASE DEFAULT
      write(*,*)'Wrong pdfOption!'
      stop
    END SELECT
!normalize and calculate cumulative distribution
    k=0
    cumul=0
    do j=1,pdfNW
      do i=1,pdfNL
        k=k+1
        cumul=cumul+pdf2D(i,j)
        cpdf2D(i,j)=cumul
      enddo
    enddo
    pdf2D=pdf2D/cumul
    cpdf2D=cpdf2D/cumul
    END SUBROUTINE
    

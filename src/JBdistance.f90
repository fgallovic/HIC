! Calculates Joyner-Boore distance from a fault
    FUNCTION jbdistance(f,staN,staE)
    IMPLICIT NONE
    REAL*8 jbdistance
    REAL*8 f(2,4),staN,staE,p1(2),p2(2),distout
    LOGICAL inside
    EXTERNAL inside
    REAL*8 t,bot,pn(2)
    INTEGER i

    if(inside(f,staN,staE))then
      jbdistance = 0.0
    else
      do i=1,4
        if(i==4)then
          p1(:) = f(:,i)
          p2(:) = f(:,1)
        else
          p1(:) = f(:,i)
          p2(:) = f(:,i+1)
        endif
        !calcola la distanza di un punto da un segmento
        if ( p1(1) == p2(1) .and. p1(2) == p2(2) )then
          t = 0.d0;
        else
          bot = sum ( ( p2(1:2) - p1(1:2) )**2 )
          t = ((staN-p1(1))*(p2(1)-p1(1))+(staE-p1(2))*(p2(2)-p1(2))) / bot;
          t = max ( t, 0.0 );
          t = min ( t, 1.0 );
        endif
        pn(:) = p1(:) + t * ( p2(:) - p1(:) );
        distout=sqrt((staN-pn(1))**2+(staE-pn(2))**2)
        if(i==1)then
          jbdistance=distout
        else
          jbdistance=min(distout,jbdistance)
        endif
      enddo
    endif
    END
    
    FUNCTION inside(f,staN,staE)
    IMPLICIT NONE
    LOGICAL inside
    REAL*8 f(2,4),staN,staE
    REAL*8 x1,y1,x2,y2,p(2)
    INTEGER i
    inside = .false.;
    p(1)=staN
    p(2)=staE
    do i=1,4
      x1=f(1,i)
      y1=f(2,i)
      if(i<4)then
        x2 = f(1,i+1)
        y2 = f(2,i+1)
      else
        x2 = f(1,1);
        y2 = f(2,1);
      endif
      if ( ( y1 < p(2) .and. p(2) <= y2 ) .or. ( p(2) <= y1 .and. y2 < p(2) ) ) then
        if ( ( p(1) - x1 ) - ( p(2) - y1 ) * ( x2 - x1 ) / ( y2 - y1 ) < 0.0d0 ) then
          inside = not(inside)
        endif
      endif
    enddo
    END
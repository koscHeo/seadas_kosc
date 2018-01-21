      subroutine fd_tausq(idayjl,xlat,xlon,tht, itau, tausq) 
        implicit none

 
        integer(4), intent(in)  :: idayjl
        real(4)   , intent(in)  :: xlat,xlon,tht
        integer(2), intent(in)  :: itau(360,180,12,91)
        real(4)   , intent(out) :: tausq

        integer(4) istart
        integer(4) i1,i2,j1,j2,k1,k2,l1,l2
        real(4) a1,a2,b1,b2,c1,c2,d1,d2,brief,xmon

        data istart/1/

        if(istart.eq.1) then
        istart=0
        if(maxval(itau).ne.10000) stop 'itau table was not read in, pgm stopped'
        endif

!      Change stop to write(*,*)  JMG 09/21/11
!      if(idayjl.lt.1 .or. idayjl.gt.366.)  stop 'idayjl oob in fd_vege_tran, pgm stopped'
!      if(abs(xlat).gt.90.)                 stop 'xlat   oob in fd_vege_tran, pgm stopped'
!      if(xlon.lt.0..or.xlon.gt.360.)       stop 'xlon   oob in fd_vege_tran, pgm stopped'
!      if(tht.lt.0. .or.tht.gt.  90.)       stop 'tht    oob in fd_vege_tran, pgm stopped'

      if(idayjl.lt.1 .or. idayjl.gt.366.)  write(*,*) 'idayjl oob in fd_vege_tran, pgm stopped'
      if(abs(xlat).gt.90.)                 write(*,*) 'xlat   oob in fd_vege_tran, pgm stopped'
      if(xlon.lt.0..or.xlon.gt.360.)       write(*,*) 'xlon   oob in fd_vege_tran, pgm stopped'
      if(tht.lt.0. .or.tht.gt.  90.)       write(*,*) 'tht    oob in fd_vege_tran, pgm stopped'

        xmon=12.*(idayjl - 0.5)/365.25
        if(xmon.gt.11.9999) xmon=11.9999

      brief=xmon-0.5
        i1=1+brief
        i2=i1+1
        a1=i1-brief
        a2=1-a1
        if(i1.eq. 0) i1=12
        if(i2.eq.13) i2= 1 
  
        brief=xlat+89.5
      j1=1+brief
      j2=j1+1
      b1=j1-brief
      b2=1.-b1
      if(j1.eq.  0) j1=  1
        if(j2.eq.181) j2=180

      brief=xlon-0.5
      k1=1+brief
      k2=k1+1
      c1=k1-brief
      c2=1-c1
        if(k1.eq.  0) k1=360
      if(k2.eq.361) k2=  1
 
        brief=tht
        if(brief.gt.89.999) brief=89.999
        l1=1+brief
      l2=l1+1
      d1=l1-brief
      d2=1-d1
 
      tausq=
     & d1*(a1*b1*(c1*itau(k1,j1,i1,l1)+c2*itau(k2,j1,i1,l1))+
     &     a1*b2*(c1*itau(k1,j2,i1,l1)+c2*itau(k2,j2,i1,l1))+
     &     a2*b1*(c1*itau(k1,j1,i2,l1)+c2*itau(k2,j1,i2,l1))+
     &     a2*b2*(c1*itau(k1,j2,i2,l1)+c2*itau(k2,j2,i2,l1)))+
     & d2*(a1*b1*(c1*itau(k1,j1,i1,l2)+c2*itau(k2,j1,i1,l2))+
     &     a1*b2*(c1*itau(k1,j2,i1,l2)+c2*itau(k2,j2,i1,l2))+
     &     a2*b1*(c1*itau(k1,j1,i2,l2)+c2*itau(k2,j1,i2,l2))+
     &     a2*b2*(c1*itau(k1,j2,i2,l2)+c2*itau(k2,j2,i2,l2)))

        tausq=1.e-4*tausq
  
        return
        end

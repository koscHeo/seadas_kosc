      subroutine fd_sss_sigma0(tht,sst,surtb_mea,em0,irad,sss)
        use salinity_module

        implicit none

        integer(2), intent(in )  :: irad
        real(4),    intent(in )  :: tht,sst,surtb_mea(2)
        real(4),    intent(out)  :: em0(2)
        real(4),    intent(out)  :: sss

        integer(4), parameter    :: itermax=51
        real(4), parameter       :: xtol=0.0005
        real(4),    parameter :: missing_val=-9999. ! indicator for missing value

        integer(4) iter,itersv
        real(4) freq,sss0,sss1,surtep,surtb0(2),surtb1(2),a,y,y0,y1,wt(3)
        
        data freq/1.413/
        
        data wt/0.,0.0,0.0/
        
        y=surtb_mea(1) - wt(irad)*surtb_mea(2)
        
        surtep=273.15+sst
        
        sss0=35
        sss1=sss0+0.5
        
        do iter=1,itermax

           if (iter==itermax) then ! not converged
              sss=missing_val
              return
           endif

           itersv=iter

           call fdem0_meissner_wentz_salinity(freq,tht,sst,sss0, em0);     surtb0=surtep*em0 
           call fdem0_meissner_wentz_salinity(freq,tht,sst,sss1, em0);     surtb1=surtep*em0 
        
           y0=surtb0(1) - wt(irad)*surtb0(2)
           y1=surtb1(1) - wt(irad)*surtb1(2)
        
           a=2.*(y1-y0)
        
           if (abs(a).lt.1.0E-8) then
              sss=missing_val
              return
           endif

           sss=sss0 + (y-y0)/a
           if(sss.lt.0) sss=0
        
           if(abs(sss-sss0).lt.xtol) exit

           sss0=sss
           sss1=sss0+0.5

        enddo                   !iter
        
        return
        end

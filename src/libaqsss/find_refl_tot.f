       subroutine find_refl_tot(idayjl,xlat,xlon,tht,phir,frac_land,frac_ice,
     1                          sss,surtep,wind,sm,staticdata,
     2                          refl_tot,tbsur) 
        use static_data_module

        implicit none

        integer(4), intent(in)  :: idayjl
        real(4),    intent(in)  :: xlat,xlon,tht,phir,frac_land,frac_ice,sss,surtep,wind,sm
        type(static_data_struct), intent(in) :: staticdata
        real(4),    intent(out) :: refl_tot(2),tbsur(2)

        real(4) sst,frac_ice_adj,frac_water,surtep_k
        real(4) refl_water(2),refl_ice(2),refl_land(2),tausq,sm_adj

        ! surtep (degrees C)
        sst=surtep
        if(sst.lt.-2) sst=-2
        if(sst.gt.35) sst=35

        ! surtep_k (degrees K)
        surtep_k = surtep + 273.15
  
        frac_ice_adj=frac_ice
        if(frac_ice.gt.1-frac_land) frac_ice_adj=1-frac_land

        frac_water=1 - frac_land - frac_ice_adj

        refl_water=0
        refl_ice  =0
        refl_land =0

        if(frac_water  .ne.0) then
              call fd_water_refl_exact(tht,sss,sst,wind,phir,staticdata,refl_water) 
        endif

        if(frac_ice_adj.ne.0) call fd_seaice_refl(tht, refl_ice) 
        
        if(frac_land   .ne.0) then
                                    
c     using 1.01 rather than 1 accounts for possible roundoff error
           sm_adj=sm
           if(sm_adj.gt.1.01) sm_adj=0.35 !fill in missing values, looking at the sm images 0.35 is a typical value for islands
           if(sm_adj.gt.1.)   sm_adj=1. !fixes any roundoff error

           ! Temp check in Kelvin
           ! Was incorrectly using C temp causing land refl to use snow
           ! rather than soil for "warmer" temps  JMG  07/15/11
           if(surtep_k.ge.273.15) then
              call fd_land_refl(1,tht,sm_adj, refl_land) !1 denotes soil
           else
              call fd_land_refl(2,tht,sm_adj, refl_land) !2 denotes snow
           endif

           call fd_tausq(idayjl,xlat,xlon,tht, staticdata%itau, tausq)

           refl_land = refl_land*0.74*tausq !0.74=exp(-.3) 0.3 is rough hght,b param as in pellarin et al. 2003
           ! Change stop to write JMG  09/20/11
           if(minval(refl_land).lt.0 .or. maxval(refl_land).gt.1) write(*,*) 'land refl oob, pgm stopped'

        endif

        refl_tot=frac_water*refl_water + frac_ice_adj*refl_ice + frac_land*refl_land

        tbsur=frac_water*(1-refl_water)*(sst+273.15) + 
     1        frac_ice_adj*(1-refl_ice)*surtep_k + 
     2        frac_land*(1-refl_land)*surtep_k
        return
        end
 

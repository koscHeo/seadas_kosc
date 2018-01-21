      subroutine adjust_tagal_ref(irad,refl0,refl,transq0,transq,ta_earth,tagal_ref, apc_matrix, apc_matrix_inv, tagal_ref_adj)
        implicit none

        integer(4), intent(in)  :: irad
        real(4),    intent(in)  :: refl0(2),refl(2),transq0,transq,ta_earth(3),tagal_ref(3)  !ta_earth has tasun_dir, tasun_ref, tagal_dir removed
        real(4),    intent(in)  :: apc_matrix(3,3,3), apc_matrix_inv(3,3,3)
        real(4),    intent(out) :: tagal_ref_adj(3)

        integer(4) istokes
        real(4) tbtoi(3),tb_gal_toa(3),tb_gal_toi(3),faraday_deg 
      real(4) tbv,tbh

      real(4) atand, cosd, sind

!     ===============================================================================================
!     ==================== find estimates of earth tb toi and faraday rotation angle ================
!     ===============================================================================================
        do istokes=1,3
        tbtoi(istokes)=dot_product(apc_matrix(:,istokes,irad),ta_earth(1:3))
        enddo

      faraday_deg=0.5*atand(tbtoi(3)/tbtoi(2))

!     =================================================================================================================
!     ==================== convert normalized tagal to tbgal at toa (tagal_ref is found assuming no faraday rotation) ====
!     =================================================================================================================
        do istokes=1,2
        tb_gal_toa(istokes)=dot_product(apc_matrix(:,istokes,irad),tagal_ref(1:3))  
        enddo
        tb_gal_toa(3)=0


!     =================================================================================================================
!     ==================== convert normalized tbgal to tbgal corresponding to actual env conditions ===================
!     =================================================================================================================

        tbv=0.5*(tb_gal_toa(1)+tb_gal_toa(2))
        tbh=0.5*(tb_gal_toa(1)-tb_gal_toa(2))

        tbv=tbv*refl(1)/refl0(1)
        tbh=tbh*refl(2)/refl0(2)

        tb_gal_toa(1)=tbv+tbh
        tb_gal_toa(2)=tbv-tbh

      tb_gal_toa(1:2)=tb_gal_toa(1:2)*transq/transq0

!     =================================================================================================================
!     =============================== apply faraday rotation to tbgal =================================================
!     =================================================================================================================

        tb_gal_toi(1)=tb_gal_toa(1)
        tb_gal_toi(2)=tb_gal_toa(2)*cosd(2*faraday_deg)
        tb_gal_toi(3)=tb_gal_toa(2)*sind(2*faraday_deg)

!     =================================================================================================================
!     =============================== convert tb toi galatic radiation to antenna temperature =========================
!     =================================================================================================================

        do istokes=1,3
        tagal_ref_adj(istokes)=dot_product(apc_matrix_inv(:,istokes,irad),tb_gal_toi(1:3))
        enddo

        return
        end

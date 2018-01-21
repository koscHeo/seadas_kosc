! RFI Detection
! T. Meissner, June 2010
! adapted from C. Ruf and S. Misra, April 2010

        subroutine get_rfi_auxiliary_data(rfi)
        use l2gen_aquarius_module
        implicit none

        integer(4), parameter :: inunit=3
        integer(4)                        :: i, j, kpol, irad

        type(rfi_struct) :: rfi

        character datadir*255
        integer len, lenstr
        call getenv('OCDATAROOT', datadir)
        len = lenstr(datadir)

!               Reading in WM
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_WM,status='old',form='formatted')
        read(inunit,223) ((rfi%WM(i,j),j=1,360),i=1,181)
                close(inunit)

!               Reading in WD
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_WD,status='old',form='formatted')
        read(inunit,223) ((rfi%WD(i,j),j=1,360),i=1,181)
                close(inunit)


!               Reading in TD
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_TD,status='old',form='formatted')
        read(inunit,224) ((rfi%TD(i,j),j=1,360),i=1,181)
                close(inunit)

!               Reading in TD_A
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_TD_A,status='old',form='formatted')
        read(inunit,224) ((rfi%TD_A(i,j),j=1,360),i=1,181)
                close(inunit)

!               Reading in TD_D
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_TD_D,status='old',form='formatted')
        read(inunit,224) ((rfi%TD_D(i,j),j=1,360),i=1,181)
                close(inunit)


!               Reading in TM
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_TM,status='old',form='formatted')
        read(inunit,224) ((rfi%TM(i,j),j=1,360),i=1,181)
                close(inunit)

!               Reading in WM_C
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_WM_C,status='old',form='formatted')
        read(inunit,223) ((rfi%WM_C(i,j),j=1,360),i=1,181)
                close(inunit)


!               Reading in WD_C
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_WD_C,status='old',form='formatted')
        read(inunit,223) ((rfi%WD_C(i,j),j=1,360),i=1,181)
                close(inunit)


!               Reading in TD_C
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_TD_C,status='old',form='formatted')
        read(inunit,224) ((rfi%TD_C(i,j),j=1,360),i=1,181)
                close(inunit)


!               Reading in TM_C
                open(unit=inunit,file=datadir(1:len)//'/'//rfi_aux_TM_C,status='old',form='formatted')
        read(inunit,224) ((rfi%TM_C(i,j),j=1,360),i=1,181)
                close(inunit)

  223   format(360I8)
  224   format(360F8.0)

!               single accumulation NEDT for typical ocean or CND observation after front end losses 
        do irad=1,n_rad
           do kpol=1,npol
              rfi%stdta_CND(kpol,irad) = (ta_CND_nom(kpol)   + receiver_noise)/bw_tau
           enddo
        enddo

        rfi%idata=1

        return
        end subroutine get_rfi_auxiliary_data 


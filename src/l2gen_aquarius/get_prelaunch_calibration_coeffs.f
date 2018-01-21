        subroutine get_prelaunch_calibration_coeffs(coeff_loss_file, coeff_nl_file, plc)

        use l2gen_aquarius_module
        implicit none

        character*256 coeff_loss_file
        character*256 coeff_nl_file

        type(prelaunch_struct) :: plc

        open(unit=3,form='formatted',file=coeff_loss_file,status='old',action='read')
        read(3,*) ! header strings
        read(3,*) ! header strings
        read(3,*) ! header strings
        read(3,*) ! header strings
        read(3,*) ! header strings
        read(3,*) ! header strings

        read(3,3001) plc%CL1    
        read(3,3001) plc%CL2A
        read(3,3001) plc%CL2B
        read(3,3001) plc%CL3,     plc%Tref_3
        read(3,3001) plc%CL4,     plc%Tref_4
        read(3,3001) plc%dL4_dT4              !parts/million
        read(3,3001) plc%CL5,     plc%Tref_5
        read(3,3001) plc%dL5_dT5              !parts/million
        read(3,3001) plc%CLMM
        read(3,3001) plc%TCND_0,    plc%Tref_CND
        read(3,3001) plc%dTCND_dT4
        read(3,3001) plc%dTCND_dT
        read(3,3001) plc%TND_0,     plc%Tref_ND
        read(3,3001) plc%dTND_dT
        read(3,3001) plc%TND_offset_0, plc%Tref_ND_offset
        read(3,3001) plc%dTND_offset_dT
        read(3,3001) plc%TDL_offset_0, plc%Tref_DL_offset
        read(3,3001) plc%dTDL_offset_dT

 3001 format(24x,1x,7(f8.4,1x))

        plc%dL4_dT4=1.0E-6*plc%dL4_dT4
        plc%dL5_dT5=1.0E-6*plc%dL5_dT5

        plc%Tref_3   =plc%Tref_3   +273.15
        plc%Tref_4   =plc%Tref_4   +273.15
        plc%Tref_5   =plc%Tref_5   +273.15
        plc%Tref_CND =plc%Tref_CND +273.15
        plc%Tref_ND  =plc%Tref_ND  +273.15
        plc%Tref_ND_offset  =plc%Tref_ND_offset +273.15
        plc%Tref_DL_offset  =plc%Tref_DL_offset +273.15

        close(3)

        open(unit=3,form='formatted',file=coeff_nl_file,status='old',action='read')
        read(3,*)       ! header strings

        read(3,3002) plc%Tref_nl
        read(3,3002) plc%Dnl20
        read(3,3002) plc%Dnl21
        read(3,3002) plc%Dnl22
        read(3,3002) plc%Dnl30
        read(3,3002) plc%Dnl31
        read(3,3002) plc%Dnl32
 3002 format(12x,12f16.1)

        close(3)

        return
        end subroutine get_prelaunch_calibration_coeffs

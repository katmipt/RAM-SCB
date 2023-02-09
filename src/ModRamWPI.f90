!============================================================================
!    Copyright (c) 2016, Los Alamos National Security, LLC
!    All rights reserved.
!============================================================================

MODULE ModRamWPI
! Contains subroutines related to wave particle interactions
! and electron lifetimes

  implicit none

  contains

!**************************************************************************
!                              WAVEPARA1
!               Life time due to WPI inside plasmasphere
!**************************************************************************
  SUBROUTINE WAVEPARA1(S)

    use ModRamMain,      ONLY: DP
    use ModRamGrids,     ONLY: NE, NR
    use ModRamVariables, ONLY: EKEV, LZ, WALOS1

    implicit none
    integer, intent(in) :: S

    integer :: i, ii, j, k
    real(DP):: TAU_WAVE,xE,xL,xlife, rEa(5),rL(8),rlife(5,8),clife(5)
    DATA rEa/0.2,0.5,1.0,1.5,2.0/
    DATA rL/5.0,4.5,4.0,3.5,3.0,2.5,2.0,1.65/

    ! life time (as a function of energy and L)
    rlife(1,1)=6.80
    rlife(1,2)=16.44
    rlife(1,3)=13.75
    rlife(1,4)=17.38
    rlife(1,5)=53.08
    rlife(1,6)=187.06
    rlife(1,7)=93.72
    rlife(1,8)=101571.57
    rlife(2,1)=23.38
    rlife(2,2)=55.98
    rlife(2,3)=43.43
    rlife(2,4)=31.75
    rlife(2,5)=38.20
    rlife(2,6)=104.90
    rlife(2,7)=164.86
    rlife(2,8)=185.67
    rlife(3,1)=343.16
    rlife(3,2)=475.15
    rlife(3,3)=99.87
    rlife(3,4)=62.46
    rlife(3,5)=98.82
    rlife(3,6)=134.95
    rlife(3,7)=171.96
    rlife(3,8)=73.63
    rlife(4,1)=619.62
    rlife(4,2)=356.89
    rlife(4,3)=139.64
    rlife(4,4)=130.32
    rlife(4,5)=210.25
    rlife(4,6)=283.46
    rlife(4,7)=359.03
    rlife(4,8)=159.19
    rlife(5,1)=1062.13
    rlife(5,2)=381.88
    rlife(5,3)=210.37
    rlife(5,4)=231.97
    rlife(5,5)=370.61
    rlife(5,6)=498.14
    rlife(5,7)=638.07
    rlife(5,8)=473.75

    ! Calculates the losses due to the w-p interaction
    DO K=2,NE
      DO II=2,NR
        xE=EKEV(K)/1000.
        xL=LZ(II)
        if (xL.ge.1.65.and.xL.le.5.0) then
          do i=8,2,-1
            if (xL.ge.rL(i).and.xL.lt.rL(i-1)) then
              do j=1,5
                clife(j)=(log10(rlife(j,i-1))-log10(rlife(j,i)))/(rL(i-1)-rL(i))*(xL-rL(i))+log10(rlife(j,i))
                clife(j)=10.**clife(j)
              enddo
              EXIT
            endif
          enddo
        elseif (xL.gt.5.0) then
          do j=1,5
            clife(j)=(log10(rlife(j,1))-log10(rlife(j,2)))/(rL(1)-rL(2))*(xL-rL(1))+log10(rlife(j,1))
            clife(j)=10.**clife(j)
          enddo
        elseif (xL.lt.1.65) then
          do j=1,5
            clife(j)=(log10(rlife(j,7))-log10(rlife(j,8)))/(rL(7)-rL(8))*(xL-rL(8))+log10(rlife(j,8))
            clife(j)=10.**clife(j)
          enddo
        endif

        if (xE.ge.0.2.and.xE.lt.2.0) then
          do i=1,4
            if (xE.ge.rEa(i).and.xE.lt.rEa(i+1)) then
              xlife=(log10(clife(i+1))-log10(clife(i)))/(log10(rEa(i+1)) &
                    -log10(rEa(i)))*(log10(xE)-log10(rEa(i)))+log10(clife(i))
              xlife=10.**xlife
              EXIT
            endif
          enddo
        elseif (xE.lt.0.2) then
          xlife=(log10(clife(2))-log10(clife(1)))/(log10(rEa(2)) &
                -log10(rEa(1)))*(log10(xE)-log10(rEa(1)))+log10(clife(1))
          xlife=10.**xlife
        elseif (xE.ge.2.0) then
          xlife=(log10(clife(5))-log10(clife(4)))/(log10(rEa(5)) &
                -log10(rEa(4)))*(log10(xE)-log10(rEa(5)))+log10(clife(5))
          xlife=10.**xlife
        endif

        tau_wave=xlife*60.*60.*24.  ! day -> sec
        WALOS1(II,K)=tau_wave
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE WAVEPARA1

!**************************************************************************
!                              WAVEPARA2
!       Another life time due to diffusion not everywhere strong
!**************************************************************************
  SUBROUTINE WAVEPARA2(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: HMIN, RE, PI
    use ModRamGrids,     ONLY: NE, NR
    use ModRamVariables, ONLY: EKEV, LZ, RLZ, V, WALOS2, WALOS3
    !!!! Module Subroutines/Functions
    use ModRamFunctions, ONLY: asind

    implicit none
    integer, intent(in) :: S

    integer :: i, k
    real(DP):: TAU_WAVE,EMEV,R1,R2,CLC
    real(DP), ALLOCATABLE :: CONE(:)

    ALLOCATe(CONE(nR+4))
    CONE = 0.0

    DO K=2,NE
      DO I=2,NR
        EMEV=EKEV(K)*0.001
        R1=0.08*EMEV**(-1.32)
        R2=0.4*10.**(2.*LZ(I)-6.+0.4*log10(29.*EMEV))
        tau_wave=min(R1,R2)
        tau_wave=1.0/tau_wave
        tau_wave=tau_wave*60.*60.*24.  ! day -> sec
        WALOS2(I,K)=tau_wave
      ENDDO
    ENDDO

    ! CONE - in degree
    DO I=1,NR
      CLC=(RE+HMIN)/RLZ(I)
      CONE(I)=ASIND(SQRT(CLC**3/SQRT(4.-3.*CLC)))
    END DO
    CONE(NR+1)=2.5
    CONE(NR+2)=1.5
    CONE(NR+3)=1.
    CONE(NR+4)=0.

    DO I=2,NR
      DO K=2,NE
        WALOS3(I,K)=64.0*LZ(I)*RE/35./(1-0.25)/ &
        SIN(CONE(I)*PI/180.)/SIN(CONE(I)*PI/180.)/V(S,K)
      ENDDO
    ENDDO

    DEALLOCATE(CONE)
    RETURN
  END SUBROUTINE WAVEPARA2

!*************************************************************************
!                              WAPARA_EMIC 
!       Routine reading normalized PA diffusion coeff
!          based on B. Ni's diffusion coefficients.
!**************************************************************************
  SUBROUTINE WAPARA_EMIC

    use ModRamVariables, ONLY: Daa_emic_h, Daa_emic_he, EKEV_emic, &
                               fp2c_emic, PAbn, DP, LZ, Ihs_emic, Ihes_emic
    use MOdRamGrids,     ONLY: ENG_emic, NCF_emic, NPA_emic, NLZ_emic, NR, NPA
    use ModIoUnit,       ONLY: UnitTMP_
    use ModRamMain,      ONLY: PathRamIn
    use ModRamGSL,       ONLY: GSL_Interpolation_1D

    implicit none
    integer :: jfpc, iR, L, iE, IER,I, ipa, iu,im,ih
    character(len=2) :: fpetofce
    character :: LL
    character*100 :: filename
    character :: ind1
    real(DP) :: Y, YY1, YY2, Lz_emic(NLZ_emic), PA_emic(NPA_emic)
    real(DP) :: Daahe_emic_tmp(NPA_emic), Daah_emic_tmp(NPA_emic)
    real(DP), allocatable :: Daahe_emic(:,:,:), Daah_emic(:,:,:)


    ! read and store the EMIC related pitch angle diffusion coefficients
    allocate(Daahe_emic(NLZ_emic, ENG_emic, NPA_emic), &
             Daah_emic( NLZ_emic, ENG_emic, NPA_emic))
    
    ! pitch angle in the Daa file
    do L=1, NPA_emic
       PA_emic(L) = L
    end do

    ! Lz in the Daa file
    do iR=1, NLZ_emic
       Lz_emic(iR) = 3+(iR-1)
    end do
    
    ! energies in the Daa file
    do iE=1, ENG_emic
       EkeV_emic(iE) = 10.**(-1.+0.1*(iE-1))
    end do
    
26  format(89(1PE16.7))
    
    ! fp2c in the Daa file
    do jfpc = 1, NCF_emic
       fp2c_emic(jfpc) = jfpc*2
       
       write(fpetofce,'(I2.2)')jfpc*2
       
       do iR=3, 7
          write(LL,'(I1.1)')iR
          filename = trim(PathRamIn)//'/HBand/EMIC_Proton_HBand_L'//LL//&
               '_fpetofce'//fpetofce//'_Daa.txt'
          
          open(unit=unittmp_,file=filename, status='old')
          do iE=1,ENG_emic
             read(unittmp_,26)(DaaH_emic_tmp(ipa), ipa=1,NPA_emic)
             ! daa_emic: interpolate it into ram pitch angle grid
             where(Daah_emic_tmp < 1.0e-31)Daah_emic_tmp = 1.0e-31
             Daah_emic(iR,iE,1) = DaaH_emic_tmp(NPA_emic) ! around 89 degree
             do L=2,NPA-1
                if(PAbn(L) .GE. 1.0 .and. PAbn(L) .LE. 89)then
                   call GSL_Interpolation_1D(PA_emic, DaaH_emic_tmp, &
                        PAbn(L), Daah_emic(iR,iE,L), IER)                  
                else
                   Daah_emic(iR,iE,L) = 1.0e-31
                end if
             end do
          end do
          
          close(unittmp_)
          
          open(unit=unittmp_,file=trim(PathRamIn)//'/HeBand/EMIC_Proton_HeBand_L'//LL//&
               '_fpetofce'//fpetofce//'_Daa.txt', status='old')
          do iE=1,ENG_emic
             read(unittmp_,26)(DaaHe_emic_tmp(ipa), ipa=1,NPA_emic)
             ! daa_emic: interpolate it into ram pitch angle grid
             where(Daahe_emic_tmp < 1.0e-31)Daahe_emic_tmp = 1.0e-31
             Daahe_emic(iR,iE,L) = DaaHe_emic_tmp(NPA_emic)                   
             do L=2,NPA-1
                if(PAbn(L) .GE. PA_emic(1) .and. PAbn(L) .LE. PA_emic(NPA_emic))then
                   call GSL_Interpolation_1D(PA_emic, DaaHe_emic_tmp, &
                        PAbn(L), Daahe_emic(iR,iE,L), IER)
                else
                   Daahe_emic(iR,iE,L) = 1.0e-31
                end if
             end do
          end do
          close(unittmp_)
       end do
       
       do iE=1,ENG_emic
          do L=2,NPA-1
             do I=2,NR
                if(maxval(Daah_emic(:,IE,L)) > 1.0e-31 .and. Lz(I) .GE.Lz_emic(1) .and. LZ(I) .LE. Lz_emic(NLz_emic))then
                   call GSL_Interpolation_1D(Lz_emic, daah_emic(:,iE, L),  Lz(I), YY1, IER)
                   Daa_emic_h(I,iE,L,jfpc) = YY1
                   
                else
                   Daa_emic_h(I,iE,L,jfpc) = 1.0e-31
                end if
                
                if(maxval(Daahe_emic(:,IE,L)) > 1.0e-31 .and. Lz(I) .GE.Lz_emic(1) .and. LZ(I) .LE. Lz_emic(NLz_emic))then
                   call GSL_Interpolation_1D(Lz_emic, daahe_emic(:,iE, L),  Lz(I), YY1, IER)
                   Daa_emic_he(I,iE,L,jfpc) = YY1
                else
                   Daa_emic_he(I,iE,L,jfpc) = 1.0e-31
                end if
                
             end do
          end do
       end do
       
    end do
    
    where(Daa_emic_h  ==0.0)Daa_emic_h = 1.0e-31
    where(Daa_emic_he ==0.0)Daa_emic_he= 1.0e-31

    deallocate(Daahe_emic, Daah_emic)


    !Read and store EMIC WAVE Intensity from Saikin model
    do iu=1,4
       write(ind1,'(I1.1)')iu
       filename =  trim(PathRamIn)//'/EMIC_model/EMIC_H_intensity_AE_'//ind1//'.txt'
       open(unit=unittmp_,file=filename, status='old')
       do im=1,20
          read(unittmp_,*)(Ihs_emic(iu,im,ih), ih=1,25)
       enddo
       close(unittmp_)

       filename = trim(PathRamIn)//'/EMIC_model/EMIC_He_intensity_AE_'//ind1//'.txt'
       open(unit=unittmp_,file=filename, status='old')
       do im=1,20
          read(unittmp_,*)(Ihes_emic(iu,im,ih), ih=1,25)
       enddo
       close(unittmp_)       
    enddo
    
  END SUBROUTINE WAPARA_EMIC

!*************************************************************************
!                              WAPARA_HISS
!       Routine reading normalized Energy & PA hiss diffusion coeff
!**************************************************************************
  SUBROUTINE WAPARA_HISS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: PathRamIn
    use ModRamGrids,     ONLY: NR, NCF, ENG, NPA
    use ModRamParams,    ONLY: HissFilePath
    use ModRamVariables, ONLY: LZ, fpofc, ndaaj, ENOR, ndvvj
    !!!! Share Modules
    use ModIoUnit,   ONLY: UNITTMP_

    implicit none
    integer, intent(in) :: S

    integer :: i, ix, kn, l
    character(len=80) HEADER
    character(len=3) ST4
    character(len=2) ST3, ST2

    ST2 = '_e'
    DO I=1,NR
      fpofc(1)=2.
      write(ST4,'(I3.3)') INT(LZ(I)*100)
      DO IX=1,NCF
        write(ST3,'(I2.2)') INT(fpofc(ix))
        OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/whis_L'//ST4//'_'//ST3//ST2//'.aan',STATUS='old')
        READ(UNITTMP_,20) HEADER
        DO KN=1,ENG
          read(UNITTMP_,17) ENOR(KN)
          read(UNITTMP_,27)(ndaaj(i,kn,l,ix),l=1,npa)
        ENDDO
        ndvvj = 0
        IF (IX.LT.NCF) fpofc(ix+1)=fpofc(ix)+4.

        DO KN=1,ENG
          DO L=1,NPA
            if (ndaaj(i,kn,l,ix).lt.1e-20) ndaaj(i,kn,l,ix)=1e-20
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(UNITTMP_)

17  FORMAT(E13.4)
20  FORMAT(A80)
27  FORMAT(80(1PE12.3))

    RETURN
  END SUBROUTINE WAPARA_HISS

!*************************************************************************
!                              WAPARA_CHORUS
!       Routine reading bounce-aver PA wave diffusion coeff
!**************************************************************************
  SUBROUTINE WAPARA_CHORUS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, PathRamIn
    use ModRamGrids,     ONLY: NR, NT, ENG, NPA
    use ModRamVariables, ONLY: KP, ECHOR, BDAAR
    !!!! Share Modules
    use ModIoUnit, ONLY: UNITTMP_

    implicit none
    integer, intent(in) :: S

    integer :: i, j, kn, l, ikp
    real(DP), ALLOCATABLE :: RLDAA(:,:),RUDAA(:,:)
    character(len=1) ST3
    character(len=2) ST2
    character(len=80) HEADER

    ST2 = '_e'

    ALLOCATE(RLDAA(ENG,NPA),RUDAA(ENG,NPA))
    RLDAA = 0.0; RUDAA = 0.0

    ikp=INT(KP)
    IF (ikp.gt.4) ikp=4
    write(ST3,'(I1.1)') ikp

    OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/wlowcho_Kp'//ST3//ST2//'.aan',STATUS='old')
    READ(UNITTMP_,20) HEADER
    print*,'in WAPARA_CHORUS: ',trim(HEADER)

    DO I=1,NR
      DO J=1,NT
        READ(UNITTMP_,20) HEADER
        DO KN=1,ENG
          read(UNITTMP_,17) ECHOR(KN)
          read(UNITTMP_,27)(RLDAA(kn,l),l=1,npa)
        ENDDO
        DO KN=1,ENG
          DO L=1,NPA
            if (RLDAA(kn,l).lt.1e-20) RLDAA(kn,l)=1e-20
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(UNITTMP_)

    OPEN(UNIT=UNITTMP_,FILE=trim(PathRamIn)//'/wuppcho_Kp'//ST3//ST2//'.aan',STATUS='old')
    READ(UNITTMP_,20) HEADER

    DO I=1,NR
      DO J=1,NT
        READ(UNITTMP_,20) HEADER
        DO KN=1,ENG
          read(UNITTMP_,17) ECHOR(KN)
          read(UNITTMP_,27)(RUDAA(kn,l),l=1,npa)
        ENDDO
        DO KN=1,ENG
          DO L=1,NPA
            if (RUDAA(kn,l).lt.1e-20) RUDAA(kn,l)=1e-20
            BDAAR(i,j,kn,l)=(RLDAA(kn,l)+RUDAA(kn,l))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CLOSE(UNITTMP_)

17  FORMAT(E13.4)
20  FORMAT(A80)
27  FORMAT(80(1PE12.3))

    DEALLOCATE(RLDAA,RUDAA)

    RETURN
  END SUBROUTINE WAPARA_CHORUS

! *************************************************************************
!                              WAPARA_Kp
!       Interpolate Kp-dependent diffusion cofficients
!**************************************************************************
  SUBROUTINE WAPARA_Kp(S)

    use ModRamVariables, ONLY: KP, CDAAR, CDAAR_chorus, NKpDiff, Kp_chorus

    implicit none
    integer, intent(in) :: S

    integer :: i1,i2

    if (Kp.gt.maxval(Kp_chorus)) then
      CDAAR(:,:,:,:) = CDAAR_chorus(:,:,1:35,:,NKpDiff)
    else
      i1 = minloc(abs(Kp-Kp_chorus),dim=1)
      if (Kp.lt.Kp_chorus(i1)) then
        i1 = i1-1
      end if
      i2 = i1+1
      ! Linear interpolation of Kp
      CDAAR(:,:,:,:) = (Kp-Kp_chorus(i1))*CDAAR_chorus(:,:,1:35,:,i2) &
                     + (Kp_chorus(i2)-Kp)*CDAAR_chorus(:,:,1:35,:,i1)
      CDAAR = CDAAR/(Kp_chorus(i2) - Kp_chorus(i1))
    end if
    CDAAR(:,25,:,:) = CDAAR(:,1,:,:)

    return
  end SUBROUTINE WAPARA_Kp

! *************************************************************************
!                              WAPARA_BAS
!       Routine reading normalized Energy & PA wave diffusion coeff
!**************************************************************************
  SUBROUTINE WAPARA_BAS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: DoUseKpDiff, BASFilePath
    use ModRamGrids,     ONLY: NPA, NT, NE, NR
    use ModRamVariables, ONLY: MU, nR_Dxx, nE_Dxx, nPa_Dxx, RCHOR_Dxx, &
                               TCHOR_Dxx, ECHOR_Dxx, CDAAR_chorus, &
                               PACHOR_Dxx, CDAAR, nKpDiff, nT_Dxx
    !!!! Module Subroutines/Functions
    use ModRamFunctions, ONLY: ACOSD
    !!!! Share Modules
    use ModIoUnit, ONLY: UNITTMP_

    implicit none
    integer, intent(in) :: S

    integer :: i,j,k,l,nkp,nloop
    character(len=32) :: H1,H2,H3,nchar
    character(len=200) :: fname
    real(DP), ALLOCATABLE :: Dxx_hold(:,:,:), PA(:)

    ALLOCATE(Dxx_hold(NR_Dxx,NE_Dxx,NPA_Dxx), PA(NPA))
    Dxx_hold = 0.0; Pa = 0.0

    write(*,*) "Starting WAPARA_BAS"

    DO L=1,NPA
      PA(L)=ACOSD(MU(L))
    END DO

    if (DoUseKpDiff) then
      nloop = NKpDiff
    else
      nloop = 1
    endif

    do nkp=1,nloop
      write(nchar,'(i1)') nkp-1
      fname = trim(BASFilePath)//'bav_diffcoef_chorus_rpa_Kp'//trim(nchar)//'.PAonly.dat'
      OPEN(UNIT=UNITTMP_,FILE=trim(fname),STATUS='old')
      ! First skip over header
      do i=1,12,1
        read(UNITTMP_,*)
      end do
      do i=1,NR_Dxx
        do j=1,NT_Dxx
          do k=1,NE_Dxx
            read(UNITTMP_,'(A19,F6.4,A9,F6.4,A12,F8.2)') H1,RCHOR_Dxx(i), H2, &
                           TCHOR_Dxx(j), H3, ECHOR_Dxx(k)
            do l=1,NPA_Dxx
              read(UNITTMP_,'(F15.6,3E18.6)') PACHOR_Dxx(l),CDAAR_chorus(i,j,k,l,nkp)
            end do
            ! skip over blank lines
            do l=1,4
              read(UNITTMP_,*)
            end do
          end do
        end do
      end do
      CLOSE(UNIT=UNITTMP_)
      ! Interpolate onto L-shell, energy and pitch angle, assuming time is
      ! normal
      if ((NT.ne.25).and.(NE.ne.35)) then
        write(*,*) "We have a problem... assuming NT=25 & NE=35"
        stop
      end if

      if ((NR_Dxx.eq.NR).and.(NPA.eq.NPA_Dxx)) then
        write(*,*) "No interpolation of diffusion coeff for", fname
      end if

    end do   ! end NKpDiff loop

    ! Initialization
    CDAAR(:,:,:,:) = CDAAR_chorus(:,:,1:35,:,1)

    write(*,*) "Finished WAPARA_BAS"

    DEALLOCATE(Dxx_hold, PA)
    RETURN
  END SUBROUTINE WAPARA_BAS

!************************************************************************
!                       WAVELO_LEGACY
!       Calculate loss due to waves everywhere using electron lifetimes
!************************************************************************
  SUBROUTINE WAVELO_LEGACY(S)

    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: DoUsePlasmasphere
    use ModRamGrids,     ONLY: NE, NR, NT, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: F2, KP, Kpmax12, LZ, IP1, IR1, EKEV, NECR, &
                                WALOS1, WALOS2, WALOS3

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l, j1, i1
    real(DP) :: TAU_LIF,Bw
    real(DP), ALLOCATABLE :: RLpp(:)

    ALLOCATE(RLpp(nT))
    RLpp = 0.0

    Bw=30.
    IF (KP.GE.4.0) Bw=100.
    DO J=1,NT
      RLpp(J)=5.39-0.382*Kpmax12  ! PP from Moldwin et al. [2002]
      DO I=2,NR
        IF (DoUsePlasmasphere.and.NECR(I,J).gt.50.) RLpp(J)=LZ(I)
      ENDDO
    ENDDO

    DO K=2,NE
      DO I=2,NR
        DO J=1,NT
          DO L=2,NPA
            IF (LZ(I).LE.RLpp(J)) THEN
              TAU_LIF=WALOS1(I,K)*((10./Bw)**2)
            ELSEIF (LZ(I).GT.RLpp(J)) THEN
              IF (EKEV(K).LE.1000.) THEN
                TAU_LIF=WALOS2(I,K)*(1+WALOS3(I,K)/WALOS2(I,K))
                if (ekev(k).le.1.1) then
                  tau_lif=tau_lif*37.5813*exp(-1.81255*ekev(k))
                elseif (ekev(k).gt.1.1.and.ekev(k).le.5.) then
                  tau_lif=tau_lif*(7.5-1.15*ekev(k))
                else
                  tau_lif=tau_lif
                endif
              ELSEIF(EKEV(K).GT.1000.) THEN
                TAU_LIF=5.*3600*24/KP
              ENDIF
            ENDIF
            F2(S,I,J,K,L)=F2(S,I,J,K,L)*EXP(-DTs/TAU_LIF)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(RLpp)
    RETURN
  END SUBROUTINE WAVELO_LEGACY

 
!************************************************************************
!                       WAVELO
!        To fix: for chorus waves we use here Orlov & Shprits (2014)
!        Note that the model is valid for R0 = 3-8,Kp = 0-6, 
!        E = 1 keV-2 MeV. For now for hight Kp we use Kp=6 losses,
!        for R0 < 3, R0 = 3, etc. 
!
!************************************************************************
  subroutine wavelo(s)

    use ModRamMain,      ONLY: DP
    use ModRamParams,    ONLY: DoUsePlasmasphere
    use ModRamGrids,     ONLY: NE, NR, NT, NPA
    use ModRamTiming,    ONLY: Dts
    use ModRamVariables, ONLY: F2, Kp, Kpmax12, LZ, MLT, EKEV, NECR, &
                                WALOS1, WALOS2, WALOS3

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l
    real(DP) :: tau_lif, logtau, Bw
    real(DP), ALLOCATABLE :: RLpp(:)
    real(DP), dimension(33) :: coeff, rke_array

    ALLOCATE(RLpp(nT))
    RLpp = 0.0
    Bw=30.
    IF (KP.GE.4.0) Bw=100.
 
    do j=1,NT
      RLpp(j) = 5.39 - 0.382*Kpmax12  ! PP from Moldwin et al. [2002]
      do i = 2,NR
        if (DoUsePlasmasphere .and. NECR(i,j) < 50.) then
            RLpp(j) = LZ(i-1)
            exit
        end if
      end do
    end do

    do k=2,NE
      do i=2,NR
        do j=1,NT
          do l=2,NPA
            if (Lz(i) <= RLpp(j)) then
              tau_lif=WALOS1(I,K)*((10./Bw)**2)
 !             E = log10(EkeV(k)*0.001)
 !             fL = 0.1328*Lz(i)**2 - 2.1463*Lz(i) + 3.7857
 !             if (E >= fL) then
 !               call hiss_tauel(Lz(i), MLT(j), E, Kp, tau_hiss)
 !               tau_lif = tau_hiss*3600*24
 !             else
 !               tau_lif = 3000.0*3600*24
 !             end if
               F2(s,i,j,k,l)=F2(s,i,j,k,l)*EXP(-DTs/tau_lif)
            else if (Lz(i) > RLpp(j)) then
              call chorus_tauel_coeff(MLT(j), EkeV(k), Kp, coeff)
              call rke_terms(Lz(i), Kp + 1, EkeV(k)*0.001, rke_array)
              logtau = sum(coeff*rke_array)
              tau_lif = 10**logtau ! lifetime in days
              tau_lif = tau_lif*3600*24
            end if
            if (MLT(j) > 15 .and. MLT(j) < 21) then 
              F2(s,i,j,k,l) = F2(s,i,j,k,l)
            else
              F2(s,i,j,k,l) = F2(s,i,j,k,l)*exp(-DTs/tau_lif)
            end if
          end do
        end do
      end do
    end do

    deallocate(RLpp)
    return

  end subroutine wavelo

  subroutine hiss_tauel(R, mlt, E, Kp, tau_hiss)
    use ModRamMain, ONLY: DP

    implicit none

    real(DP), intent(in) :: R, mlt, E, Kp
    real(DP), intent(out) :: tau_hiss

    real(DP) :: tau_hiss_av, g, h

    call hiss_tauel_av(R, E, tau_hiss_av)
    call gmlt(mlt, g)
    call hKp(Kp, h)
    tau_hiss = tau_hiss_av/(g*h)

  end subroutine hiss_tauel

  subroutine hiss_tauel_av(R, E, tau_hiss_av)
    use ModRamMain, ONLY: DP

    implicit none

    real(DP), intent(in) :: R, E
    real(DP), intent(out) :: tau_hiss_av

    real(DP), dimension(20) :: a_coeff, EL_coeff

    a_coeff(:) = (/ 77.323, -92.641, -55.754, 44.497, 48.981, 8.9067, -10.704,&
                  -15.711, -3.3326, 1.5189, 1.294, 2.2546, 0.31889, -0.85916,&
                  -0.22182, 0.034318, 0.097248, -0.12192, -0.062765, 0.0063218 /)

    EL_coeff(:) = (/ 1.0, R, E, R**2, R*E, E**2, R**3, R**2*E, R*E**2, E**3,&
                   R**4, R**3*E, R**2*E**2, R*E**3, E**4, R*E**4, R**2*E**3,&
                   R**4*E, R**5, E**5 /)

    tau_hiss_av = 10**(sum(a_coeff*EL_coeff))
  end subroutine hiss_tauel_av

  subroutine gmlt(mlt, g)
    use ModRamMain, ONLY: DP

    implicit none

    real(DP), intent(in) :: mlt
    real(DP), intent(out) :: g

    real(DP) :: g0mlt, G0, b2, b1, b0

    G0 = 782.3
    b2 = -0.007338
    b1 = 0.1773
    b0 = 2.080

    g0mlt = b2*mlt**2 + b1*mlt + b0
    g = (1/G0) * (10**g0mlt)

  end subroutine gmlt

  subroutine hKp(Kp, h)
    use ModRamMain, ONLY: DP

    implicit none

    real(DP), intent(in) :: Kp
    real(DP), intent(out) :: h

    real(DP) :: h0Kp, H0, c2, c1, c0, K

    H0 = 1315
    c2 = -0.01414
    c1 = 0.2321
    c0 = 2.598

    if (Kp <= 5) then
      K = Kp
    else
      K = 5
    end if
    
    h0Kp = c2*K**2 + c1*K + c0
    h = (1/H0) * (10**h0Kp)

  end subroutine hKp
  

  !***********************************************************************
  !   Calculate coefficients for Orlova & Shprits (2014) model of
  !   elctron lifetimes due to scattering by chorus waves 
  !***********************************************************************
  subroutine chorus_tauel_coeff(mlt, E, Kp, coeff)
    use ModRamMain, ONLY: DP

    implicit none

    real(DP), intent(in) :: mlt, E, Kp
    real(DP), dimension(33), intent(out) :: coeff

    real(DP) :: coeffNight(5,33), coeffDawn(3,33), coeffPrenoon(3,33),&
                coeffPostnoon(3,33)

    coeffNight(1,:) = (/ -6.298, 2.475, -0.72639, 86.973, -1571.8, -0.015397,&
                       -1.253, 0.0011662, -319.94, -5.9866, 0.0, -0.32751,&
                       0.090541, -4.7742, 0.0015693, 11.896, 1073.2, 0.01501,&
                       -0.0039606, 0.0, 2.1137, -424.7, 26441.0, -988010.0,&
                       -0.088225, 8.6697, 419.1, 0.0047009, -0.79409, 2025.4,&
                       -147880.0, 5629200.0, 0.0 /)
    
    coeffNight(2,:) = (/ 9.2456, -3.7685, 1.2559, 0.040296, 1.1833, -0.10843,&
                       0.0, 0.013698, 16.075, -19.866, 8.134, 0.54406, -0.16616,&
                       -0.091841, 0.0, -1.6971, 0.99102, -0.028145, 0.0092606,&
                       0.058892, -3.5674, -9.0573, -17.686, 10.029, 0.6322,&
                       7.368, 0.0, -0.085533, -1.0005, -31.398, 88.895, -59.737, 0.0 /)
                      
    coeffNight(3,:) = (/ -3.8064, -1.0616, -0.50282, -0.93478, 0.0, 0.15969,&
                       -0.081421, -0.0081014, 27.262, -10.554, 2.3565, 0.4421,&
                       -0.059054, 0.15542, -0.0027346, -4.05, 0.9056, -0.028424,&
                       0.0049602, 0.14818, 4.1008, -55.812, 10.392, 0.0, -0.85125,&
                       9.9087, -0.83409, 0.044016, -0.53404, 82.36, -76.472,&
                       158.86, -148.75 /)

    coeffNight(4,:) = (/ 9.4941, -4.608, 1.9814, -0.35458, 0.10723, 0.15431,&
                       -0.015408, -0.011301, 3.1246, -0.067317, -0.088908,&
                       0.91137, -0.37474, 0.0063232, -0.0063599, -0.36853, 0.0,&
                       -0.061111, 0.021043, 0.018011, -7.9589, 2.0, 0.93541,&
                       -0.42571, 2.037, 1.0634, 0.0, -0.30541, -0.13225, 9.441,&
                       -12.231, 5.0008, -0.5057 /)
                       
    coeffNight(5,:) = (/ -22.631, 11.693, -4.6909, -0.24701, 0.033458, 0.98105,&
                       0.0, -0.055868, 3.8379, -0.85239, 0.0, -0.090658, -0.15869,&
                       0.0090052, -0.0079675, -0.35965, 0.068028, -0.034443,&
                       0.018795, 0.0051132, 10.647, -11.332, 0.13953, -0.06843,&
                       -2.1683, 2.0908, 0.0, 0.12855, -0.12191, 19.956, 0.39078,&
                       -0.27118, 0.19286 /)
                            


    coeffDawn(1,:) = (/ 13.023, -5.1397, 0.36629, -114.33, 5231.2, -0.004326,&
                      1.7987, 0.0, 2364.7, -237100.0, 8168000.0, 0.66154,&
                      -0.035386, 4.475, 0.0, -223.63, 9035.4, -0.029805,&
                      0.0013554, 8.0424, -1.788, 568.06, -63536.0, 2307700.0,&
                      0.076799, -16.838, 1334.8, 0.0, 0.0, -6518.8, 923510.0,&
                      -42788000.0, 0.0 /)

    coeffDawn(2,:) = (/ 3.0872, 0.44092, -0.34824, 2.835, -10.838, 0.013725,&
                    -0.06885, 0.0, 23.283, -426.74, 1219.8, -0.14119, 0.039748,&
                    -0.100111, -0.00074614, -1.0401, 24.774, 0.010743, -0.0016182,&
                    -0.091657, 0.015684, -12.221, 22.634, 396.6, 0.0037368, 0.79364,&
                    -3.8069, 0.0011614, 0.0, -170.6, 4947.6, -49439.0, 181330.0 /)

    coeffDawn(3,:) = (/ 4.0112, -0.39094, -0.03256, -0.037959, -0.0062679, 0.0052271,&
                      0.001263, -0.00017601, 0.10032, 0.07421, 0.011347, 0.040335,&
                      0.0055705, 0.003349, -0.00038591, -0.02609, -0.0065557,&
                      -0.0024257, -0.0002716, 0.0015915, -0.93675, 0.057597,&
                      0.056129, -0.0087605, 0.059057, -0.0068603, 0.0014772,&
                      0.0010756, -0.00031907, 3.4193, -2.9549, 1.257, -0.24123 /)   

    
    coeffPrenoon(1,:) = (/ 9.3791, -3.3801, 0.2295, -21.044, 3237.3, -0.014983,&
                         -1.3745, 0.0, 1956.3, -260190.0, 12551000.0, 0.40549,&
                         -0.016404, 0.0, 0.0014003, -180.88, 6908.2, -0.018219,&
                         0.0, 7.9534, -1.4482, 280.49, -32982.0, 0.0, 0.10963,&
                         -19.436, 2541.5, 0.0, 0.0, -4987.2, 843540.0, -44455000.0, 0.0 /)
    
    coeffPrenoon(2,:) = (/ 6.8865, -1.6632, -0.018769, -0.4048, -4.0194, -0.0047973,&
                         0.0, 0.00049759, 37.493, -184.03, 0.0, 0.14847, 0.0099013,&
                         0.088663, 0.0, -3.6665, 16.119, -0.0025242, -0.00098309,&
                         0.045456, -0.30072, 3.9614, -49.548, 376.0, 0.018935,&
                         -0.073137, 0.0, -0.0012688, 0.0, -217.92, 3729.6,&
                         -34125.0, 127070.0 /)
    
    coeffPrenoon(3,:) = (/ 0.37772, 0.31909, -0.097241, 0.052829, -0.0054941, 0.01252,&
                         -0.0045574, -0.00014461, -0.11659, 0.090307, 0.004568,&
                         -0.056582, 0.0079453, -0.001957, -0.00090294, -0.016334,&
                         -0.0049014, 0.0014537, 0.0, 0.0030733, 0.051883, -0.23224,&
                         0.083002, -0.018977, -0.033829, 0.018878, -0.001079,&
                         0.00042082, 0.0, 4.7071, -4.0683, 1.8424, -0.33505 /)

    
    coeffPostnoon(1,:) = (/5.4927, -1.2239, 0.1676, -32.412, 4524.0, -0.013993,&
                          -1.1361, 0.0 ,1181.8, -211750.0, 13923000.0, 0.079532 ,&
                          -0.0075973, 0.0, 0.0009959, -81.463, 1780.8, -0.0038066,&
                          0.0, 4.754, -1.4015, 478.46, -84497.0, 3569700.0, 0.14088,&
                          -26.347, 3597.9, -0.0022453 ,0.0, -3681.7, 858490.0,&
                          -60968000.0, 0.0 /)

    coeffPostnoon(2,:) = (/ 5.2583, -0.99398, -0.10343, 0.2155, -5.8931, -0.0060298,&
                          0.01832, 0.00051372, 34.572, -167.45, -208.86, 0.06291 ,&
                          0.022672, 0.040741, 0.0, -3.9241, 19.809, 0.001487,&
                          -0.0015407, 0.05622, -0.050128, -0.13125, -12.771, 293.38,&
                          0.02575, -0.016043, -1.2529, -0.0017619, 0.0, -174.2,&
                          3267.4, -32412.0, 130380.0 /)

    coeffPostnoon(3,:) = (/ 0.85672, 0.41383, -0.11487, 0.052481, -0.0057176,&
                          0.013318, -0.0052097, -0.00012685, -0.066638, 0.15784,&
                          -0.0042214, -0.075598, 0.010877, -0.0015561, -0.00097701,&
                          -0.035035, -0.008202, 0.0025751, -0.00015869, 0.0044586,&
                          0.14084, -0.25711, 0.092002, -0.017677, -0.036853,&
                          0.025471, -0.0033818, 0.00032439, 0.0, 4.8355, -4.46,&
                          1.9887, -0.35605 /)

    if (mlt >= 21 .or. mlt <= 3) then
       
       if (E <= 10) then
          coeff(:) = coeffNight(1,:)
       else if (E > 10 .and. E < 500 .and. Kp <= 3) then
          coeff(:) = coeffNight(2,:)
       else if (E > 10 .and. E < 500 .and. Kp > 3) then
          coeff(:) = coeffNight(3,:)
       else if (E >= 500 .and. Kp <= 3) then
          coeff(:) = coeffNight(4,:)
       else if (E >= 500 .and. Kp > 3) then
          coeff(:) = coeffNight(5,:)
       end if
       
    else if (mlt > 3 .and. mlt <=9) then
       
       if (E < 7) then
          coeff(:) = coeffDawn(1,:)
       else if(E >= 7 .and. E < 90) then
          coeff(:) = coeffDawn(2,:)
       else if (E >=90) then
          coeff(:) = coeffDawn(3,:)
       end if
       
    else if (mlt > 9 .and. mlt <= 12) then
       
       if (E < 7) then
          coeff(:) = coeffPrenoon(1,:)
       else if (E >= 7 .and. E < 100) then
          coeff(:) = coeffDawn(2,:)
       else if (E >= 100) then
          coeff(:) = coeffDawn(3,:)
       end if
       
    else if (mlt > 12 .and. mlt <=15) then
       
       if (E < 6) then
          coeff(:) = coeffPostnoon(1,:)
       else if (E>=6 .and. E < 100) then
          coeff(:) = coeffPostnoon(2,:)
       else if (E >= 100) then
          coeff(:) = coeffPostnoon(3,:)
       end if
       
    else
       coeff(:) = 0
    end if
    
  end subroutine chorus_tauel_coeff

  subroutine rke_terms(R, Kp1, E, rke_array)
    use ModRamMain,      ONLY: DP

    implicit none

    real(DP), intent(in) :: R, Kp1, E
    real(DP), dimension(33), intent(out) :: rke_array
    real(DP) :: K, R1, E1

    if (Kp1 > 7) then
       K = 7
    else
       K = Kp1
    end if
    
    if (R < 3) then
        R1 = 3
    else
        R1 = R
    end if

    if (E < 1e-3) then
        E1 = 1e-3
    else
        E1 = E
    end if
    
    rke_array(:) = (/ 1.0, R1, R1*K, R1*K*E1, R1*K*E1**2, R1*K**2, R1*K**2*E1, R1*K**3, R1*E1,&
                    R1*E1**2, R1*E1**3, R1**2, R1**2*K, R1**2*K*E1, R1**2*K**2, R1**2*E1,&
                    R1**2*E1**2, R1**3, R1**3*K, R1**3*E1, K, K*E1, K*E1**2,&
                    K*E1**3, K**2, K**2*E1, K**2*E1**2, K**3, K**3*E1, E1,&
                    E1**2, E1**3, E1**4 /)
    
  end subroutine rke_terms
  
      
   
!*************************************************************************
!                               WPADIF
!     Routine calculates the decay of the distribution function
!        due to WPI pitch angle diffusion using implicit scheme
!*************************************************************************
  SUBROUTINE WPADIF(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, PathRamOut
    use ModRamGrids,     ONLY: NT, NR, NE, NPA
    use ModRamTiming,    ONLY: Dts, T
    use ModRamParams,    ONLY: DoUseEMIC
    use ModRamVariables, ONLY: F2, FNHS, MU, DMU, WMU, ATAC, ATAW, species,&
                               ATAW_emic_h, ATAW_emic_he        
    !!!! Share Modules
    use ModIoUnit, ONLY: UNITTMP_

    implicit none

    integer, intent(in) :: S
    integer :: i, j, k, l
    real(DP) :: AN,BN,GN,RP,DENOM
    real(DP), ALLOCATABLE :: F(:), RK(:), RL(:), FACMU(:)

    ALLOCATE(F(NPA),RK(NPA),RL(NPA),FACMU(NPA))
    F = 0.0; RK = 0.0; RL = 0.0; FACMU = 0.0

    DO J=1,NT
      DO I=2,NR
        DO K=2,NE
          DO L=2,NPA
            FACMU(L)=FNHS(I,J,L)*MU(L)
            F(L)=F2(S,I,J,K,L)/FACMU(L)
          END DO
          FACMU(1)=FNHS(I,J,1)*MU(1)
          F(1)=F(2) ! lower b.c.
          RK(1)=0.
          RL(1)=-1.
          DO L=2,NPA-1
            IF(species(S)%s_name.eq.'Electron')then
               AN=(ATAW(I,J,K,L)+ATAC(I,J,K,L))/DMU(L)       ! Hiss & chorus
               GN=(ATAW(I,J,K,L-1)+ATAC(I,J,K,L-1))/DMU(L-1) !  "
            ELSE
               IF(DoUseEMIC)then
                  AN = (ATAW_emic_h(I,J,K,L)+ATAW_emic_he(I,J,K,L))/DMU(L)       ! EMIC waves for ions
                  GN = (ATAW_emic_h(I,J,K,L-1)+ATAW_emic_he(I,J,K,L-1))/DMU(L-1) !  "   
               END IF
            END IF
            AN=AN*DTs/FACMU(L)/WMU(L)
            GN=GN*DTs/FACMU(L)/WMU(L)
            BN=AN+GN
            if (abs(-1-bn).lt.(abs(an)+abs(gn))) then
              open(UNITTMP_,file=trim(PathRamOut)//'diffcf_e.dat',status='unknown',position='append')
              write(UNITTMP_,*) ' T=',T/3600,' hr'
              write(UNITTMP_,*) 'i=',i,' j=',j,' k=',k,' l=',l
              write(UNITTMP_,*) 'an=',AN,' -1-bn=',(-1-BN),' gn=',GN
              close(UNITTMP_)
            endif
            RP=F(L)
            DENOM=BN+GN*RL(L-1)+1
            RK(L)=(RP+GN*RK(L-1))/DENOM
            RL(L)=-AN/DENOM
          END DO
          F2(S,I,J,K,NPA-1)=RK(NPA-1)/(1+RL(NPA-1)) ! upper b.c.
          DO L=NPA-2,1,-1
            F2(S,I,J,K,L)=RK(L)-RL(L)*F2(S,I,J,K,L+1)
          END DO
          F2(S,I,J,K,NPA)=F2(S,I,J,K,NPA-1)
          DO L=1,NPA
            F2(S,I,J,K,L)=F2(S,I,J,K,L)*FACMU(L)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(F,RK,RL,FACMU)
    RETURN
  END SUBROUTINE WPADIF

!*************************************************************************
!                               EMIC wave amplitude
!                           Based on Empirical Model by Saikin et al.
!*************************************************************************
  subroutine I_emic(AEind, xl, nmlt, I_H, I_He)

    use ModKind
    use ModRamVariables, ONLY: Ihs_emic, Ihes_emic
    implicit none

    integer, intent(in) :: xl, nmlt, AEind
    real(kind=real8_), intent(out) :: I_H,I_He

    I_H=0.0
    I_He=0.0

    if(AEind.ge.0.and.AEind.lt.100) then
       I_H = Ihs_emic(1,xl,nmlt)
       I_He= Ihes_emic(1,xl,nmlt)

    elseif(AEind.ge.100.and.AEind.lt.300) then
       I_H = Ihs_emic(2,xl,nmlt)
       I_He= Ihes_emic(2,xl,nmlt)
       
    elseif(AEind.ge.300.and.AEind.lt.400) then
       I_H = Ihs_emic(3,xl,nmlt)
       I_He= Ihes_emic(3,xl,nmlt)
       
    elseif(AEind.ge.400) then
       I_H = Ihs_emic(4,xl,nmlt)
       I_He= Ihes_emic(4,xl,nmlt)
       
    endif

  end subroutine I_emic
  
!==============================================================================  
END MODULE ModRamWPI

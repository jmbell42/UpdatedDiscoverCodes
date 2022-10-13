!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\
! ------------------------------------------------------------
! advance
! ------------------------------------------------------------
!/

subroutine advance_vertical_1d

  use ModVertical
  use ModPlanet, only: Mass
  use ModGITM, ONLY : Dt, iCommGITM, iProc, iUp_
  use ModInputs, only: UseBarriers, iDebugLevel, Is1D
  implicit none
  !-----------------------------------------------------------
  integer :: iError, iAlt, iSpecies, iDir
  !!!!! Variables for the Runga-Kutta 4th Order Time-stepping
  real :: OrigLogNS(-1:nAlts+2,1:nSpecies)
  real :: OrigLogINS(-1:nAlts+2,1:nIons)
  real :: OrigLogRho(-1:nAlts+2)
  real :: OrigVel_GD(-1:nAlts+2,1:3)
  real :: OrigTemp(-1:nAlts+2)
  real :: OrigVS(-1:nAlts+2,1:nSpecies)

  real :: UpdatedLogNS(-1:nAlts+2,1:nSpecies)
  real :: UpdatedLogINS(-1:nAlts+2,1:nIons)
  real :: UpdatedLogRho(-1:nAlts+2)
  real :: UpdatedVel_GD(-1:nAlts+2,1:3)
  real :: UpdatedTemp(-1:nAlts+2)
  real :: UpdatedVS(-1:nAlts+2,1:nSpecies)

  real :: FinalLogNS(-1:nAlts+2,1:nSpecies)
  real :: FinalLogINS(-1:nAlts+2,1:nIons)
  real :: FinalLogRho(-1:nAlts+2)
  real :: FinalVel_GD(-1:nAlts+2,1:3)
  real :: FinalTemp(-1:nAlts+2)
  real :: FinalVS(-1:nAlts+2,1:nSpecies)

!!! RK-4 Coefficients
  real :: K1LogNS(-1:nAlts+2,1:nSpecies)
  real :: K1LogINS(-1:nAlts+2,1:nIons)
  real :: K1LogRho(-1:nAlts+2)
  real :: K1Vel_GD(-1:nAlts+2,1:3)
  real :: K1Temp(-1:nAlts+2)
  real :: K1VS(-1:nAlts+2,1:nSpecies)

  real :: K2LogNS(-1:nAlts+2,1:nSpecies)
  real :: K2LogINS(-1:nAlts+2,1:nIons)
  real :: K2LogRho(-1:nAlts+2)
  real :: K2Vel_GD(-1:nAlts+2,1:3)
  real :: K2Temp(-1:nAlts+2)
  real :: K2VS(-1:nAlts+2,1:nSpecies)

  real :: K3LogNS(-1:nAlts+2,1:nSpecies)
  real :: K3LogINS(-1:nAlts+2,1:nIons)
  real :: K3LogRho(-1:nAlts+2)
  real :: K3Vel_GD(-1:nAlts+2,1:3)
  real :: K3Temp(-1:nAlts+2)
  real :: K3VS(-1:nAlts+2,1:nSpecies)

  real :: K4LogNS(-1:nAlts+2,1:nSpecies)
  real :: K4LogINS(-1:nAlts+2,1:nIons)
  real :: K4LogRho(-1:nAlts+2)
  real :: K4Vel_GD(-1:nAlts+2,1:3)
  real :: K4Temp(-1:nAlts+2)
  real :: K4VS(-1:nAlts+2,1:nSpecies)

!!! RStage-4 Coefficients
  real :: Stage1LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage1LogINS(-1:nAlts+2,1:nIons)
  real :: Stage1LogRho(-1:nAlts+2)
  real :: Stage1Vel_GD(-1:nAlts+2,1:3)
  real :: Stage1Temp(-1:nAlts+2)
  real :: Stage1VS(-1:nAlts+2,1:nSpecies)

  real :: Stage2LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage2LogINS(-1:nAlts+2,1:nIons)
  real :: Stage2LogRho(-1:nAlts+2)
  real :: Stage2Vel_GD(-1:nAlts+2,1:3)
  real :: Stage2Temp(-1:nAlts+2)
  real :: Stage2VS(-1:nAlts+2,1:nSpecies)

  real :: Stage3LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage3LogINS(-1:nAlts+2,1:nIons)
  real :: Stage3LogRho(-1:nAlts+2)
  real :: Stage3Vel_GD(-1:nAlts+2,1:3)
  real :: Stage3Temp(-1:nAlts+2)
  real :: Stage3VS(-1:nAlts+2,1:nSpecies)

  real :: Stage4LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage4LogINS(-1:nAlts+2,1:nIons)
  real :: Stage4LogRho(-1:nAlts+2)
  real :: Stage4Vel_GD(-1:nAlts+2,1:3)
  real :: Stage4Temp(-1:nAlts+2)
  real :: Stage4VS(-1:nAlts+2,1:nSpecies)

  real :: DtOriginal
  real :: DtIn

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 6) write(*,*) "=======> vertical bcs 1", iproc

  DtOriginal = Dt  !!! Store this so that it doesn't change

!!! =================
!!! General RK4 Update:
!!! Y(n+1) = Y(n) + Dt/6*(k1 + 2k2 + 2k3 + k4)
!!! Time(n+1) = Time(n) + Dt
!!! 
!!! k1 = f(tn,yn)
!!! k2 = f(tn + Dt/2, Yn + Dt/2*k1)
!!! k3 = f(tn + Dt/2, Yn + Dt/2*k2)
!!! k4 = f(tn + Dt, Yn + Dt*k3)
!!! =================

  ! Step 1, Fill in Ghost Cells
  call set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp,LogINS,IVel,VertVel)

!!! Set the Original State -> Orig = Y(n)
   OrigLogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  OrigLogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  OrigLogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  OrigVel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    OrigTemp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      OrigVS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage1LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage1LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage1LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage1Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage1Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage1VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage2LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage2LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage2LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage2Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage2Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage2VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage3LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage3LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage3LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage3Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage3Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage3VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage4LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage4LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage4LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage4Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage4Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage4VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel

  DtIn = 0.5*Dt  !!! Store this so that it doesn't change
  call advance_vertical_1stage(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)


!!! note that Stage 1 -> updated by a 1/2 step
!!! (NewLogNS - LogNS) = f(tn + Dt/2, yn + dt/2)
       Stage2LogNS(1:nAlts,1:nSpecies)      = &
       Stage1LogNS(1:nAlts,1:nSpecies)      + &
     (NewLogNS(1:nAlts,1:nSpecies) - LogNS(1:nAlts,1:nSpecies))

      Stage2LogINS(1:nAlts,1:nIons)  = &
      Stage1LogINS(1:nAlts,1:nIons)  + &
    (NewLogINS(1:nAlts,1:nIons) - LogINS(1:nAlts,1:nIons))

      Stage2LogRho(1:nAlts)                = &
      Stage1LogRho(1:nAlts)                + &
    (NewLogRho(1:nAlts) - LogRho(1:nAlts))

      Stage2Vel_GD(1:nAlts,1:3)            = & 
      Stage1Vel_GD(1:nAlts,1:3)            + & 
    (NewVel_GD(1:nAlts,1:3) - Vel_GD(1:nAlts,1:3))

        Stage2Temp(1:nAlts)                  = &
        Stage1Temp(1:nAlts)                  + &
     (NewTemp(1:nAlts) -  Temp(1:nAlts))

          Stage2VS(1:nAlts,1:nSpecies)      = &
          Stage1VS(1:nAlts,1:nSpecies)      + &
   (NewVertVel(1:nAlts,1:nSpecies) - VertVel(1:nAlts,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage2
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage2Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage2LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage2LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage2LogRho(-1:nAlts+2) 

   UpdatedTemp(-1:nAlts+2) = &
         Stage2Temp(-1:nAlts+2) 

   UpdatedVS(-1:nAlts+2,1:nSpecies) = &
         Stage2VS(-1:nAlts+2,1:nSpecies) 

  DtIn = 0.5*Dt
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                        UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

   LogNS  = UpdatedLogNS
   LogINS = UpdatedLogINS
   LogRho = UpdatedLogRho
   Vel_GD = UpdatedVel_GD
     Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS  = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
     NewTemp = Temp
  NewVertVel = VertVel

!!!!! Calculate K2
  call advance_vertical_1stage(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)


       Stage3LogNS(1:nAlts,1:nSpecies)      = &
       Stage1LogNS(1:nAlts,1:nSpecies)      + &
     (NewLogNS(1:nAlts,1:nSpecies) - LogNS(1:nAlts,1:nSpecies))

      Stage3LogINS(1:nAlts,1:nIons)  = &
      Stage1LogINS(1:nAlts,1:nIons)  + &
    (NewLogINS(1:nAlts,1:nIons) - LogINS(1:nAlts,1:nIons))

      Stage3LogRho(1:nAlts)                = &
      Stage1LogRho(1:nAlts)                + &
    (NewLogRho(1:nAlts) - LogRho(1:nAlts))

      Stage3Vel_GD(1:nAlts,1:3)            = & 
      Stage1Vel_GD(1:nAlts,1:3)            + & 
    (NewVel_GD(1:nAlts,1:3) - Vel_GD(1:nAlts,1:3))

        Stage3Temp(1:nAlts)                  = &
        Stage1Temp(1:nAlts)                  + &
      (NewTemp(1:nAlts) -  Temp(1:nAlts))

          Stage3VS(1:nAlts,1:nSpecies)      = &
          Stage1VS(1:nAlts,1:nSpecies)      + &
   (NewVertVel(1:nAlts,1:nSpecies) - VertVel(1:nAlts,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage3
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage3Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage3LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage3LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage3LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage3Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
            Stage3VS(-1:nAlts+2,1:nSpecies) 

  DtIn = Dt
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)
!
  LogNS  = UpdatedLogNS
  LogINS = UpdatedLogINS
  LogRho = UpdatedLogRho
  Vel_GD = UpdatedVel_GD
  Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel
!
!
!!!!!! Calculate K3
  call advance_vertical_1stage(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
!!!! K3 Coefficients for RK-4

       Stage4LogNS(1:nAlts,1:nSpecies)      = &
       Stage1LogNS(1:nAlts,1:nSpecies)      + &
     (NewLogNS(1:nAlts,1:nSpecies) - LogNS(1:nAlts,1:nSpecies))

      Stage4LogINS(1:nAlts,1:nIons)  = &
      Stage1LogINS(1:nAlts,1:nIons)  + &
    (NewLogINS(1:nAlts,1:nIons) - LogINS(1:nAlts,1:nIons))

      Stage4LogRho(1:nAlts)                = &
      Stage1LogRho(1:nAlts)                + &
    (NewLogRho(1:nAlts) - LogRho(1:nAlts))

      Stage4Vel_GD(1:nAlts,1:3)            = & 
      Stage1Vel_GD(1:nAlts,1:3)            + & 
    (NewVel_GD(1:nAlts,1:3) - Vel_GD(1:nAlts,1:3))

        Stage4Temp(1:nAlts)                  = &
        Stage1Temp(1:nAlts)                  + &
      (NewTemp(1:nAlts) -  Temp(1:nAlts))

          Stage4VS(1:nAlts,1:nSpecies)      = &
          Stage1VS(1:nAlts,1:nSpecies)      + &
   (NewVertVel(1:nAlts,1:nSpecies) - VertVel(1:nAlts,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage4
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage4Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage4LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage4LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage4LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage4Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
          Stage4VS(-1:nAlts+2,1:nSpecies) 


!!!! Update Boundary Conditions
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

  LogNS  = UpdatedLogNS
  LogINS = UpdatedLogINS
  LogRho = UpdatedLogRho
  Vel_GD = UpdatedVel_GD
  Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel
  
!! Calculate K4 (Final Coefficient)

  DtIn = 0.5*Dt
  call advance_vertical_1stage(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
  
  ! This section ensures that our lower boundary conditions are maintained
  ! and not overwritten.
  FinalLogNS  = Stage1LogNS
  FinalLogINS = Stage1LogINS
  FinalLogRho = Stage1LogRho
  FinalVel_GD = Stage1Vel_GD
  FinalTemp   = Stage1Temp
  FinalVS     = Stage1VS

!!! Set the Updated State:  Stage 2
  FinalLogNS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1LogNS(1:nAlts,:) + Stage2LogNS(1:nAlts,:)  +  &
                 2.0*Stage3LogNS(1:nAlts,:) + Stage4LogNS(1:nAlts,:)  +  &
                       (NewLogNS(1:nAlts,:) - LogNS(1:nAlts,:)) )

  FinalLogINS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1LogINS(1:nAlts,:) + Stage2LogINS(1:nAlts,:)  +  &
                 2.0*Stage3LogINS(1:nAlts,:) + Stage4LogINS(1:nAlts,:)  +  &
                       (NewLogINS(1:nAlts,:) - LogINS(1:nAlts,:)) )

  FinalLogRho(1:nAlts) = &
     (1.0/3.0)*(-1.0*Stage1LogRho(1:nAlts) + Stage2LogRho(1:nAlts)  +  &
                 2.0*Stage3LogRho(1:nAlts) + Stage4LogRho(1:nAlts)  +  &
                       (NewLogRho(1:nAlts) - LogRho(1:nAlts)) )

  FinalVel_GD(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1Vel_GD(1:nAlts,:) + Stage2Vel_GD(1:nAlts,:)  +  &
                 2.0*Stage3Vel_GD(1:nAlts,:) + Stage4Vel_GD(1:nAlts,:)  +  &
                       (NewVel_GD(1:nAlts,:) - Vel_GD(1:nAlts,:)) )

  FinalTemp(1:nAlts) = &
     (1.0/3.0)*(-1.0*Stage1Temp(1:nAlts) + Stage2Temp(1:nAlts)  +  &
                 2.0*Stage3Temp(1:nAlts) + Stage4Temp(1:nAlts)  +  &
                       (NewTemp(1:nAlts) - Temp(1:nAlts)) )

  FinalVS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1VS(1:nAlts,:) + Stage2VS(1:nAlts,:)  +  &
                 2.0*Stage3VS(1:nAlts,:) + Stage4VS(1:nAlts,:)  +  &
                  (NewVertVel(1:nAlts,:) - VertVel(1:nAlts,:)) )

  if(.not. Is1D) then
     call set_vertical_bcs_relax(FinalLogRho, FinalLogNS, FinalVel_GD, &
                                  FinalTemp, FinalLogINS, IVel, FinalVS)
  else
     call set_vertical_bcs(FinalLogRho, FinalLogNS, FinalVel_GD, &
                             FinalTemp, FinalLogINS, IVel, FinalVS)
  endif

   LogNS = FinalLogNS
  LogINS = FinalLogINS
  LogRho = FinalLogRho
  Vel_GD = FinalVel_GD
    Temp = FinalTemp
 VertVel = FinalVS

   VerticalHeating1D(1:nAlts) = &
       (FinalTemp(1:nAlts) - OrigTemp(1:nAlts))/Dt
  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d

!=============================================================================
subroutine advance_vertical_1stage( DtIn, &
     LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
     LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! With fluxes and sources based on LogRho..Temp, update NewLogRho..NewTemp

  use ModGITM, only: &
       Dt, iEast_, iNorth_, iUp_, ThermalDiffCoefS
  use ModPlanet
  use ModSizeGitm
  use ModVertical, only : &
       Heating, EddyCoef_1D, ViscCoef_1d, Coriolis, &
       MeanMajorMass_1d, Gamma_1d, InvRadialDistance_C, &
       ChemSources_1d, &
       Centrifugal, &
       CD_MeshCoefs, &
       HydroHeatCool_1d, &
       AdiabaticHeatCool_1d, &
       Gravity_G, Altitude_G,Cv_1D, dAlt_F
  use ModTime
  use ModInputs
  use ModConstants
  use ModSources, only : EddyCondAdia
  implicit none

  real, intent(in) :: DtIn
  real, intent(in) :: LogRho(-1:nAlts+2)
  real, intent(in) :: LogNS(-1:nAlts+2,nSpecies)
  real, intent(in) :: LogINS(-1:nAlts+2,nIons)
  real, intent(in) :: Vel_GD(-1:nAlts+2,3)
  real, intent(in) :: IVel(-1:nAlts+2,3)
  real, intent(in) :: Temp(-1:nAlts+2)
  real, intent(in) :: VertVel(-1:nAlts+2,nSpecies)

  real, intent(inout) :: NewLogRho(-1:nAlts+2)
  real, intent(inout) :: NewLogNS(-1:nAlts+2,nSpecies)
  real, intent(inout) :: NewLogINS(-1:nAlts+2,nIons)
  real, intent(inout) :: NewVel_GD(-1:nAlts+2,3)
  real, intent(inout) :: NewTemp(-1:nAlts+2)
  real, intent(out) :: NewVertVel(-1:nAlts+2,nSpecies)

  integer :: iAlt, iSpecies, jSpecies, iDim

  real :: NS(-1:nAlts+2,nSpecies)
  real :: Rho(-1:nAlts+2)
  real, dimension(1:nAlts)    :: DivVel, GradTemp, DiffTemp, &
                                GradTmp, DiffTmp, DiviVel
  real, dimension(1:nAlts,3) :: GradVel_CD, DiffVel_CD
  real, dimension(1:nAlts,3) :: GradiVel_CD, DiffiVel_CD
  real, dimension(1:nAlts,nIons) :: GradLogINS, DiffLogINS
  real :: NewSumRho, NewLogSumRho
  real, dimension(-1:nAlts+2)    :: NT, Press

  real :: nVel(1:nAlts,1:nSpecies)
  integer :: nFilter, iFilter
  real :: LowFilter

!! WAVEDRAG Heating  Hickey et al [2000]
  real, dimension(1:nAlts)    :: StressHeating
!\
! Parameters Used for the Sponge
! This Sponge is useful to dampen out spurious modes
! oscillating between the bottom and top of the model.

  integer ::   nAltsSponge = 12
  integer :: nAltsRayleigh = 12
  real :: kSP, NuSP, AmpSP
  real :: kRay, NuRay, AmpRay

  ! JMB:  Adding Eddy Diffusion Variables here
  ! Note:  These are used in the calc_neutral_friction
  !--------------------------------------------------------------------------
  !! Eddy Diffusion Variables
  real, dimension(1:nAlts,nSpecies)    :: GradLogConS
  real, dimension(-1:nAlts+2,nSpecies)    :: ConS, LogConS
  real, dimension(1:nAlts,nSpecies)    :: EddyCoefRatio_1d

  ! ----------------------------------------------------
  ! JMB:  AUSM Variables
  ! ----------------------------------------------------
  real ::    RhoS(-1:nAlts+2,1:nSpecies),&
          NewRhoS(-1:nAlts+2,1:nSpecies),&
      AUSMRhoSFluxes(1:nAlts,1:nSpecies)

  ! Momentum Fluxes and Variables
  real ::    MomentumS(-1:nAlts+2,1:nSpecies),&
          NewMomentumS(-1:nAlts+2,1:nSpecies),&
              Momentum(-1:nAlts+2,1:3),&        ! Bulk Momentum
           NewMomentum(-1:nAlts+2,1:3),&        ! Bulk Momentum
       AUSMMomentumSFluxes(1:nAlts,1:nSpecies), &
       AUSMMomentumFluxes(1:nAlts,3)

  real :: PressureS(-1:nAlts+2,1:nSpecies), &
          NewNS(-1:nAlts+2,1:nSpecies), &
          NewNT(-1:nAlts+2)

  real :: TotalEnergy(-1:nAlts+2),&
       NewTotalEnergy(-1:nAlts+2),&
       AUSMTotalEnergyFluxes(1:nAlts), &
       NewPress(-1:nAlts+2), NewRho(-1:nAlts+2)

  real :: RadialDistance_C(-1:nAlts+2)
  real :: EffectiveGravity(-1:nAlts+2)
  real :: EffectiveGravity_NoCoriolis(-1:nAlts+2)

  ! JMB:  Use these as Limiters on Winds for an initial startup
  real :: TimeFactor, Vel0, DeltaV, VelocityCap

  !--------------------------------------------------------------------------

  Vel0 = 1.0 ! initial velocity in m/s
  TimeFactor = exp(-tSimulation/10.0/86400.0)
  !TimeFactor = 0.0
  DeltaV = MaximumVerticalVelocity - Vel0
  VelocityCap = Vel0 + DeltaV*(1.0 - TimeFactor)

  !TimeFactor  = exp(-tSimulation/10.0/86400.0)
  !TimeFactor  = 0.0
  nAltsSponge = 0 + floor(nAlts*TimeFactor) 

  do iAlt = -1, nAlts + 2
     EffectiveGravity(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) + & 
        (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
        * InvRadialDistance_C(iAlt) + & 
        Coriolis * Vel_GD(iAlt,iEast_)

     EffectiveGravity_NoCoriolis(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) 
  enddo 

  NS = exp(LogNS)
  Rho = exp(LogRho)
  nFilter = 10
  do iAlt = -1, nAlts + 2
     RadialDistance_C(iAlt) = 1.0/InvRadialDistance_C(iAlt)
  enddo 
  
  NT(-1:nAlts+2) = 0.0
  do iSpecies = 1, nSpecies
     NT(-1:nAlts+2) = NT(-1:nAlts+2) + &
                      NS(-1:nAlts+2,iSpecies)
  enddo 

  do iAlt = -1, nAlts + 2
    Press(iAlt) = NT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
  enddo

  call calc_rusanov_alts(Temp   ,GradTemp,    DiffTemp)
  do iDim = 1, 3
     call calc_rusanov_alts(Vel_GD(:,iDim), &
          GradVel_CD(:,iDim),DiffVel_CD(:,iDim))
     call calc_rusanov_alts(iVel(:,iDim), &
          GradiVel_CD(:,iDim),DiffiVel_CD(:,iDim))
  enddo

  ! Add geometrical correction to gradient and obtain divergence
  DivVel = GradVel_CD(:,iUp_) + 2*Vel_GD(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)
  DiviVel = GradiVel_CD(:,iUp_) + 2*iVel(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)

  do iSpecies=1,nIons-1
     call calc_rusanov_alts(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo

  !! -------- Check the Method for Eddy Diffusion 
  ! Use Colegrove Method
  ! Step 1:  Calculate Ln(rho_{s}/Rho)
  do iSpecies = 1, nSpecies
    LogConS(-1:nAlts+2,iSpecies) = &
         alog(Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)/Rho(-1:nAlts+2))
  enddo 

  do iAlt = 1, nAlts
     do iSpecies = 1, nSpecies
        GradLogConS(iAlt,iSpecies) =  &
           CD_MeshCoefs(iAlt,1)*LogConS(iAlt-2,iSpecies)&
        +  CD_MeshCoefs(iAlt,2)*LogConS(iAlt-1,iSpecies)&
        +  CD_MeshCoefs(iAlt,3)*LogConS(iAlt  ,iSpecies)&
        +  CD_MeshCoefs(iAlt,4)*LogConS(iAlt+1,iSpecies)&
        +  CD_MeshCoefs(iAlt,5)*LogConS(iAlt+2,iSpecies)
     enddo  ! iSpecies Loop
  enddo ! iAlt Loop

  do iSpecies = 1, nSpecies
     RhoS(-1:nAlts+2,iSpecies) =  &
        Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)
     NewRhoS(-1:nAlts+2,iSpecies) = RhoS(-1:nAlts+2,iSpecies)
     MomentumS(-1:nAlts+2,iSpecies) =  &
        Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)*&
         VertVel(-1:nAlts+2,iSpecies)
     PressureS(-1:nAlts+2,iSpecies) =  &
        NS(-1:nAlts+2,iSpecies)*Temp(-1:nAlts+2)*&
        Boltzmanns_Constant
  enddo 

  do iAlt = -1, nAlts + 2
   TotalEnergy(iAlt) = &
        Press(iAlt)/(Gamma_1d(iAlt) - 1.0) + &
        0.5*Rho(iAlt)*(Vel_GD(iAlt,iUp_)**2.0 + &
                       Vel_GD(iAlt,iNorth_)**2.0 + &
                       Vel_GD(iAlt,iEast_)**2.0) 
  enddo 
  do iDim = 1, 3
     Momentum(-1:nAlts+2,iDim) = Rho(-1:nAlts+2)*Vel_GD(-1:nAlts+2,iDim)
  enddo 
!--------------------------
  NewRho = Rho
  NewPress = Press
  NewTotalEnergy = TotalEnergy
  ! Call the AUSM Solvers

  call calc_all_fluxes_hydro(DtIn, RhoS, PressureS,  &
        AUSMRhoSFluxes,AUSMMomentumSFluxes, &
        AUSMTotalEnergyFluxes, AUSMMomentumFluxes, RadialDistance_C )


  AmpSP  = (1.0/(10.0*DtIn))
  kSP    = nAltsSponge + 1

  do iAlt = 1,nAlts
     do iSpecies=1,nSpecies
        NewRhoS(iAlt,iSpecies) = RhoS(iAlt,iSpecies) &
               - DtIn*(AUSMRhoSFluxes(iAlt,iSpecies))
        NewLogNS(iAlt,iSpecies) = alog( NewRhoS(iAlt,iSpecies)/Mass(iSpecies) )
     enddo

     do iSpecies=1,nIonsAdvect
        if (UseImprovedIonAdvection) then
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
                (DiviVel(iAlt) * LogINS(iAlt,iSpecies) + &
                IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + DtIn * DiffLogINS(iAlt,iSpecies)
        else
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
                (IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + DtIn * DiffLogINS(iAlt,iSpecies)
!
!           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
!               ( DiviVel(iAlt)  + &
!                IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
!                + DtIn * DiffLogINS(iAlt,iSpecies)
        endif
     enddo
  enddo !iAlt = 1,nAlts

  NewNS  = 0.0
  NewNT  = 0.0
  NewRho = 0.0

  do iAlt = -1, nAlts+2
     do iSpecies = 1, nSpecies
         NewNS(iAlt,iSpecies) = NewRhoS(iAlt,iSpecies)/Mass(iSpecies)

         NewRho(iAlt) = NewRho(iAlt) + &
                  NewRhoS(iAlt,iSpecies)

         NewNT(iAlt) = NewNT(iAlt) + &
               NewNS(iAlt,iSpecies)
     enddo 
  enddo 


  if (UseDamping) then
     VertTau(iAlt) = &
          15 - (1 - exp(-1.0*altitude_G(ialt)/1000.0/40.0))*5.0
  endif

  do iAlt = 1,nAlts
     do iSpecies=1,nSpecies
           NewMomentumS(iAlt,iSpecies) = MomentumS(iAlt,iSpecies) - &
                 DtIn*(AUSMMomentumSFluxes(iAlt,iSpecies)) + &
                 DtIn*RhoS(iAlt,iSpecies)*EffectiveGravity(iAlt) 

        NewMomentumS(iAlt,iSpecies) = NewMomentumS(iAlt,iSpecies) + &
              DtIn*ChemSources_1d(iAlt,iSpecies)*&
              VertVel(iAlt,iSpecies)*Mass(iSpecies)
!------------------
        NewVertVel(iAlt,iSpecies) = &
           NewMomentumS(iAlt,iSpecies)/NewRhoS(iAlt,iSpecies) 

        ! Thermal Diffusion Effects (For Light Species H2, H, and He) 
        ! ThermalDiffCoefS is set in calc_rates
        ! Note:  ThermalDiffCoefS is < 0.0 for light species
        ! This forces light species into hot zones and heavy species into cold zones
        NewVertVel(iAlt,iSpecies) = NewVertVel(iAlt,iSpecies) - &
          DtIn*(ThermalDiffCoefS(iSpecies)*Boltzmanns_Constant*&
             GradTemp(iAlt))/Mass(iSpecies)
     enddo ! iSpecies

  enddo ! iAlt

  if (UseNeutralFriction) then
     nVel(1:nAlts,1:nSpecies) = NewVertVel(1:nAlts,1:nSpecies)
     call calc_neutral_friction(DtIn,nVel(1:nAlts,1:nSpecies), &
                         EddyCoef_1d(1:nAlts), &
                               NewNT(1:nAlts), &
                               NewNS(1:nAlts,1:nSpecies), &
                         GradLogConS(1:nAlts,1:nSpecies), &
                                Temp(1:nAlts))
     NewVertVel(1:nAlts,1:nSpecies) = nVel(1:nAlts,1:nSpecies)
  endif 


!  do iAlt = -1, nAlts+2
!     do iSpecies=1,nSpecies
!        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!     enddo
!     NewVertVel(iAlt, iN4S_) = max(-10.0, NewVertVel(iAlt, iN4S_))
!     NewVertVel(iAlt, iN4S_) = min( 10.0, NewVertVel(iAlt, iN4S_))
!  enddo

 
  ! Lower Boundary Rayleigh Friction (removes spurious oscillations from the lower boundary)
  !do iAlt = -1, nAlts+2
  !   if (iAlt <= nAltsRayleigh) then
  !      NuRay = AmpRay*(1.0 - cos( pi*(kRay - iAlt)/kRay) )
  !   else
  !      NuRay = 0.0
  !   endif
  !   NewVertVel(iAlt,1:nSpecies) = NewVertVel(iAlt,1:nSpecies)/&
  !      (1.0 + DtIn*NuRay)
  !enddo 


  NewVel_GD(-1:nAlts+2,iUp_) = 0.0
  do iAlt = -1, nAlts+2
     do iSpecies=1,nSpecies
!        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!
        NewVertVel(iAlt, iSpecies) = max(-VelocityCap, &
             NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( VelocityCap, &
             NewVertVel(iAlt, iSpecies))
!
        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies)*NewRhoS(iAlt,iSpecies)/NewRho(iAlt)
     enddo

     NewVertVel(iAlt, iN4S_) = max(-10.0, NewVertVel(iAlt, iN4S_))
     NewVertVel(iAlt, iN4S_) = min( 10.0, NewVertVel(iAlt, iN4S_))
  enddo


  StressHeating = 0.0
  if (UseStressHeating) then
     do iAlt = 1, nAlts 
       StressHeating(iAlt) = ViscCoef_1d(iAlt)* &
        (  (  (Gamma_1d(iAlt) - 1.0)/ ( NT(iAlt)*Boltzmanns_Constant) ) * &
            (  &
               (4.0/3.0)*GradVel_CD(iAlt,iUp_)**2 +    &
                         GradVel_CD(iAlt,iNorth_)**2 + &
                         GradVel_CD(iAlt,iEast_)**2    &
            )  )
     enddo
  endif

  do iAlt = 1, nAlts

      NewMomentum(iAlt,iEast_) = Momentum(iAlt,iEast_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iEast_)) 
 
      NewVel_GD(iAlt,iEast_) = NewMomentum(iAlt,iEast_)/NewRho(iAlt)
      NewVel_GD(iAlt,iEast_) = NewVel_GD(iAlt,iEast_) &
          + (0.10 + 0.40*TimeFactor)*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
            (20.0e+03)) * DtIn*DiffVel_CD(iAlt,iEast_)
! 
      NewMomentum(iAlt,iNorth_) = Momentum(iAlt,iNorth_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iNorth_)) 
 
      NewVel_GD(iAlt,iNorth_) = NewMomentum(iAlt,iNorth_)/NewRho(iAlt)

      NewVel_GD(iAlt,iNorth_) = NewVel_GD(iAlt,iNorth_) &
          + (0.10 + 0.40*TimeFactor)*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
            (20.0e+03)) * DtIn*DiffVel_CD(iAlt,iNorth_)

     ! dT/dt = -(V.grad T + (gamma - 1) T div V +  &
     !        (gamma - 1) * g  * grad (KeH^2  * rho) /rho 
     ! AUSM Method
     NewTotalEnergy(iAlt)   = TotalEnergy(iAlt) - &
         DtIn*AUSMTotalEnergyFluxes(iAlt) + &
         DtIn*Rho(iAlt)*Vel_GD(iAlt,iUp_)*&
         EffectiveGravity_NoCoriolis(iAlt) 


     NewPress(iAlt) = &
        (NewTotalEnergy(iAlt) - &
         0.5*NewRho(iAlt)*(NewVel_GD(iAlt,iUp_)**2.0 + &
                           NewVel_GD(iAlt,iNorth_)**2.0 + &
                           NewVel_GD(iAlt,iEast_ )**2.0 ) )*&
                           (Gamma_1d(iAlt) - 1.0)

     NewTemp(iAlt) = NewPress(iAlt)/(Boltzmanns_Constant*NewNT(iAlt))

     AdiabaticHeatCool_1d(iAlt) = &
     AdiabaticHeatCool_1d(iAlt) + &
         -1.0*DtIn*Press(iAlt)*DivVel(iAlt)

     HydroHeatCool_1d(iAlt) = &
     HydroHeatCool_1d(iAlt) + &
         -1.0*DtIn*Rho(iAlt)*Cv_1d(iAlt)*Vel_GD(iAlt,iUp_)*&
               GradTemp(iAlt)

     NewTemp(iAlt) = NewTemp(iAlt) + &
       (0.10 + 0.40*TimeFactor)*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
       (20.0e+03)) * DtIn*DiffTemp(iAlt)

  enddo ! iAlt


  do iAlt = -1, nAlts+2
     if (iAlt >= (nAlts - nAltsSponge)) then
        NuSP = AmpSP*(1.0 - cos( pi*(kSP - (nAlts - iAlt))/kSP) )
     else
        NuSP = 0.0
     endif
     NewVertVel(iAlt,1:nSpecies) = NewVertVel(iAlt,1:nSpecies) - &
         DtIn*NuSP*VertVel(iAlt,1:nSpecies)

     NewVel_GD(iAlt,1:3) = NewVel_GD(iAlt,1:3) - &
         DtIn*NuSP*Vel_GD(iAlt,1:3)
  enddo 

  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

end subroutine advance_vertical_1stage

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov_alts(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

  call calc_facevalues_alts(Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt_C(1:nAlts)

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt_C(1:nAlts)

end subroutine calc_rusanov_alts

!\
! ------------------------------------------------------------
! calc_facevalues_alts
! ------------------------------------------------------------
!/

subroutine calc_facevalues_alts(Var, VarLeft, VarRight)

  use ModVertical, only: dAlt_F, InvDAlt_F
  use ModSizeGITM, only: nAlts
  use ModLimiterGitm

  implicit none
  
  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: VarLeft(1:nAlts+1), VarRight(1:nAlts+1)

  real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)

  real, parameter :: Factor1=0.6250000 ! 15/24
  real, parameter :: Factor2=0.0416667 !  1/24
  real :: h

  integer :: i

  do i=1,nAlts

     ! 4th order scheme for calculating face values

     h  = InvDAlt_F(i+1)*2.0
     dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
     h  = InvDAlt_F(i)*2.0
     dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))

!     ! This is Gabor's scheme
!     dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
!     dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)

     dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

!     write(*,*) dVarUp, dVarDown, dVarLimited(i)

  end do

  i = 0
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  i = nAlts+1
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  do i=1,nAlts+1
     VarLeft(i)  = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i)
     VarRight(i) = Var(i)   - 0.5*dVarLimited(i)   * dAlt_F(i)
  end do

end subroutine calc_facevalues_alts


subroutine calc_all_fluxes_hydro(DtIn, RhoS, PressureS, &
            RhoSFlux, MomentumSFlux, &
            EnergyFlux, MomentumFlux, &
            RadDist)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax, VertVel, Gamma_1D, &
                          LogRho, Vel_GD, MeanMajorMass_1d, &
                          CellVol1D, Area_P12, Area_M12, &
                          Temp, Altitude_G, dAlt_F


  use ModPlanet, only : nSpecies, nIonsAdvect, Mass, RBody, &
                        iN2_, cSpecies
  use ModGITM, only : Dt,iUp_, iEast_, iNorth_
  use ModConstants, only : Boltzmanns_Constant

  implicit none

! Passed from Vertical_Solver In
  real, intent(in) :: DtIn
  real, intent(in) :: RhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) :: PressureS(-1:nAlts+2,1:nSpecies)
  real, intent(out):: RhoSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: MomentumSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: EnergyFlux(1:nAlts)
  real, intent(out):: MomentumFlux(1:nAlts,1:3)
  real, intent(in) :: RadDist(-1:nAlts+2)

  real, dimension(1:nAlts,1:nSpecies) :: RhoSLeft_M12, RhoSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: RhoSLeft_P12, RhoSRight_P12

  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogRhoS
  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogPS

  real, dimension(1:nAlts,1:nSpecies) :: LogRhoSLeft_M12, LogRhoSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogRhoSLeft_P12, LogRhoSRight_P12

  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: RhoILeft_M12, RhoIRight_M12
  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: RhoILeft_P12, RhoIRight_P12
  real, dimension(-1:nAlts+2,1:nIonsAdvect) :: LogRhoI
  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: LogRhoILeft_M12, LogRhoIRight_M12
  real, dimension( 1:nAlts  ,1:nIonsAdvect) :: LogRhoILeft_P12, LogRhoIRight_P12


  real, dimension(1:nAlts,1:nSpecies) :: PressureSLeft_M12, PressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: PressureSLeft_P12, PressureSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: LogPressureSLeft_M12, LogPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogPressureSLeft_P12, LogPressureSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: VelLeft_M12, VelRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: VelLeft_P12, VelRight_P12

  real, dimension(1:nAlts,3) :: VelGDLeft_M12, VelGDRight_M12
  real, dimension(1:nAlts,3) :: VelGDLeft_P12, VelGDRight_P12

  real, dimension(1:nAlts,3) :: IVelLeft_M12, IVelRight_M12
  real, dimension(1:nAlts,3) :: IVelLeft_P12, IVelRight_P12

  real, dimension(1:nAlts) :: PLeft_M12, PRight_M12
  real, dimension(1:nAlts) :: PLeft_P12, PRight_P12

  real, dimension(1:nAlts) :: GammaLeft_M12, GammaRight_M12
  real, dimension(1:nAlts) :: GammaLeft_P12, GammaRight_P12

  real, dimension(1:nAlts) :: ELeft_M12, ERight_M12
  real, dimension(1:nAlts) :: ELeft_P12, ERight_P12

  real, dimension(1:nAlts) :: RhoLeft_M12, RhoRight_M12
  real, dimension(1:nAlts) :: RhoLeft_P12, RhoRight_P12

  real, dimension(1:nAlts) :: CSLeft_M12, CSRight_M12
  real, dimension(1:nAlts) :: CSLeft_P12, CSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: RhoSFlux_M12, RhoSFlux_P12
  real, dimension(1:nAlts,1:nSpecies) :: MomentumSFlux_M12, MomentumSFlux_P12
  real, dimension(1:nAlts,3) :: Momentum_M12, Momentum_P12
  real, dimension(1:nAlts,1:nIonsAdvect) :: RhoIFlux_M12, RhoIFlux_P12

  real, dimension(1:nAlts) :: EnergyFlux_M12, EnergyFlux_P12

  real :: SubCs
  integer :: iSpecies, iAlt, iDim
  !------------------------------------------------------------

  ! ==================== AUSM Flux Variables

  real, dimension( 1:nAlts,1:nSpecies) :: NumericalVelocity_P12, &
                                          NumericalVelocity_M12
  real, dimension( 1:nAlts,1:nSpecies) :: NumericalPressure_P12, &
                                          NumericalPressure_M12   

  real, dimension( 1:nAlts) :: BulkNumericalVelocity_P12, &
                               BulkNumericalVelocity_M12    
  real, dimension( 1:nAlts) :: BulkIVel_P12, BulkIVel_M12    
  real, dimension( 1:nAlts) :: MeanCS_P12, MeanCS_M12    

  real, dimension( 1:nAlts, 1:nSpecies) :: MeanPressureS_P12,&
                                           MeanPressureS_M12    
  real, dimension( 1:nAlts, 1:nSpecies) :: MeanRhoS_P12, &
                                           MeanRhoS_M12    

  real :: Kp(1:nSpecies), Ku(1:nAlts,1:nSpecies)
  real :: LiouKp, LiouKu

  real :: LiouKpS(1:nAlts,1:nSpecies), LiouKuS(1:nAlts,1:nSpecies)
  real :: MaxKpS(1:nAlts,1:nSpecies), MaxKuS(1:nAlts,1:nSpecies)
  real :: MinKpS(1:nAlts,1:nSpecies), MinKuS(1:nAlts,1:nSpecies)
  real :: KpWidth, KpAltMidPoint
  integer :: AltIndex

  real, dimension( 1:nAlts) :: LeftRadius, RightRadius    
  real, dimension( 1:nAlts) :: AreaFunction_P12, AreaFunction_M12    
  real, dimension( 1:nAlts) :: LocalCellVolume

!!!!!! ====================  JMB:  Liou Stuff ===============================
  real, dimension(1:nAlts) :: LiouCSLeft_M12, LiouCSRight_M12
  real, dimension(1:nAlts) :: LiouCSLeft_P12, LiouCSRight_P12

  real, dimension(1:nAlts) :: LiouEnthalpyLeft_M12, LiouEnthalpyRight_M12
  real, dimension(1:nAlts) :: LiouEnthalpyLeft_P12, LiouEnthalpyRight_P12

  real, dimension(1:nAlts) :: InterfaceCS_M12, InterfaceCS_P12

  real, dimension(1:nAlts,1:nSpecies) :: MLeft_M12, MRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: MLeft_P12, MRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: M2Bar_M12, M2Bar_P12
  real, dimension(1:nAlts,1:nSpecies) :: M2Zero_M12, M2Zero_P12

  real, dimension(1:nAlts,1:nSpecies) :: MZero_M12, MZero_P12
  real:: MInf, LiouBeta
  real, dimension(1:nAlts,1:nSpecies):: ModifiedZeta

  real, dimension(1:nAlts,1:nSpecies) :: FA_M12, FA_P12

!!! Polynomial Mach Functions

!!! First Order Polynomial
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Left_M12, MF1N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Right_M12, MF1N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Left_P12, MF1N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Right_P12, MF1N_Right_P12

!!!! 2nd Order Polynomial
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Left_M12, MF2N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Right_M12, MF2N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Left_P12, MF2N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Right_P12, MF2N_Right_P12

!!!! 4th Order Polynomial
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Left_M12, MF4N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Right_M12, MF4N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Left_P12, MF4N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Right_P12, MF4N_Right_P12

!!!! Pressure Mach Number
  real, dimension(1:nAlts,1:nSpecies) :: MPress_M12, MPress_P12

!!!! Interface Mach Number
  real, dimension(1:nAlts,1:nSpecies) :: InterfaceMach_M12, InterfaceMach_P12

  real, dimension(1:nAlts,1:nSpecies) :: LiouNumericalVelocity_M12, &
                                         LiouNumericalVelocity_P12

!!! Extract Mean Atmosphere Values

  real :: LogKpS(1:nAlts,1:nSpecies), LogKuS(1:nAlts,1:nSpecies)
  real :: LogMinKpS(1:nAlts,1:nSpecies), LogMinKuS(1:nAlts,1:nSpecies)
  real :: LogMaxKpS(1:nAlts,1:nSpecies), LogMaxKuS(1:nAlts,1:nSpecies)
!
  real :: ZVar, VL, VR, CBar, ML, MR
!
  real :: LiouKpS_P12(1:nAlts,1:nSpecies), LiouKuS_P12(1:nAlts,1:nSpecies)
  real :: LiouKpS_M12(1:nAlts,1:nSpecies), LiouKuS_M12(1:nAlts,1:nSpecies)
  real :: LogKpS_P12(1:nAlts,1:nSpecies), LogKuS_P12(1:nAlts,1:nSpecies)
  real :: LogKpS_M12(1:nAlts,1:nSpecies), LogKuS_M12(1:nAlts,1:nSpecies)
  real :: Alt_P12, Alt_M12
  real :: MachScaling

  MInf = 1.0e-19
  LiouBeta = 1.0/8.0

  LogRhoS(-1:nAlts+2,1:nSpecies) = alog(RhoS(-1:nAlts+2,1:nSpecies))
  LogPS(-1:nAlts+2,1:nSpecies) = alog(PressureS(-1:nAlts+2,1:nSpecies))

  Kp(1:nSpecies) = 0.10             !! Ullrich et al. [2011]
  Ku(1:nAlts,1:nSpecies) = 2.0             !! Ullrich et al. [2011]

!!! 1-D Settings
  MinKuS(1:nAlts,1:nSpecies) = 0.25
  MaxKuS(1:nAlts,1:nSpecies) = 0.50

! These Settings Looked Okay
! 03/01/2021:  This worked well
!  MinKpS(1:nAlts,1:nSpecies) = 1.0e-09
!  MaxKpS(1:nAlts,1:nSpecies) = 1.0e-02
! This looked pretty good
!  MinKpS(1:nAlts,1:nSpecies) = 1.0e-06
!  MaxKpS(1:nAlts,1:nSpecies) = 1.0e-02

  ! The Original Settings
!  MinKuS(1:nAlts,1:nSpecies) = 2.00
!  MaxKuS(1:nAlts,1:nSpecies) = 2.00
!  MinKpS(1:nAlts,1:nSpecies) = 1.0e-04
!  MaxKpS(1:nAlts,1:nSpecies) = 1.0e-03

  !! The Original Settings
  !MinKuS(1:nAlts,1:nSpecies) = 2.00
  !MaxKuS(1:nAlts,1:nSpecies) = 3.00
  !MinKpS(1:nAlts,1:nSpecies) = 1.0e-04
  !MaxKpS(1:nAlts,1:nSpecies) = 0.25
  ! The Original Settings
  MinKuS(1:nAlts,1:nSpecies) = 2.00
  MaxKuS(1:nAlts,1:nSpecies) = 2.00
  MinKpS(1:nAlts,1:nSpecies) = 1.0e-09
  MaxKpS(1:nAlts,1:nSpecies) = 1.00e-03



  !!KpWidth = 100.0e+03
  KpWidth = 10.0e+03
  KpAltMidPoint = Altitude_G(1)

  LogMaxKuS = alog(MaxKuS)
  LogMinKuS = alog(MinKuS)

  LogMaxKpS = alog(MaxKpS)
  LogMinKpS = alog(MinKpS)

  do iAlt = 1, nAlts
     do iSpecies = 1, nSpecies
     
        ! Alt_P12 = 0.5*(Altitude_G(iAlt) + Altitude_G(iAlt+1))
        ! Alt_M12 = 0.5*(Altitude_G(iAlt) + Altitude_G(iAlt-1))
        ! LogKpS_P12(iAlt,iSpecies) = LogMinKpS(iAlt,iSpecies) +  &
        !       0.5*(LogMaxKpS(iAlt,iSpecies) - LogMinKpS(iAlt,iSpecies))*&
        !       ( 1.0 + tanh(  (Alt_P12- KpAltMidPoint)/KpWidth ) )
!
!         LogKpS_M12(iAlt,iSpecies) = LogMinKpS(iAlt,iSpecies) +  &
!               0.5*(LogMaxKpS(iAlt,iSpecies) - LogMinKpS(iAlt,iSpecies))*&
!               ( 1.0 + tanh(  (Alt_M12 - KpAltMidPoint)/KpWidth ) )
!
!
!         LogKuS_P12(iAlt,iSpecies) = LogMinKuS(iAlt,iSpecies) +  &
!               0.5*(LogMaxKuS(iAlt,iSpecies) - LogMinKuS(iAlt,iSpecies))*&
!               ( 1.0 + tanh(  (Alt_P12- KpAltMidPoint)/KpWidth ) )
!
!         LogKuS_M12(iAlt,iSpecies) = LogMinKuS(iAlt,iSpecies) +  &
!               0.5*(LogMaxKuS(iAlt,iSpecies) - LogMinKuS(iAlt,iSpecies))*&
!               ( 1.0 + tanh(  (Alt_M12 - KpAltMidPoint)/KpWidth ) )
!

     
         Alt_P12 = 0.5*(Altitude_G(iAlt) + Altitude_G(iAlt+1))
         Alt_M12 = 0.5*(Altitude_G(iAlt) + Altitude_G(iAlt-1))

         LogKpS_P12(iAlt,iSpecies) = LogMinKpS(iAlt,iSpecies) +  &
               (LogMaxKpS(iAlt,iSpecies) - LogMinKpS(iAlt,iSpecies))*&
               (tanh(  (Alt_P12- KpAltMidPoint)/KpWidth ) )

         LogKpS_M12(iAlt,iSpecies) = LogMinKpS(iAlt,iSpecies) +  &
               (LogMaxKpS(iAlt,iSpecies) - LogMinKpS(iAlt,iSpecies))*&
               (tanh(  (Alt_M12 - KpAltMidPoint)/KpWidth ) )

!
         LogKuS_P12(iAlt,iSpecies) = LogMinKuS(iAlt,iSpecies) +  &
               (LogMaxKuS(iAlt,iSpecies) - LogMinKuS(iAlt,iSpecies))*&
               (tanh(  (Alt_P12- KpAltMidPoint)/KpWidth ) )

         LogKuS_M12(iAlt,iSpecies) = LogMinKuS(iAlt,iSpecies) +  &
               (LogMaxKuS(iAlt,iSpecies) - LogMinKuS(iAlt,iSpecies))*&
               (tanh(  (Alt_M12 - KpAltMidPoint)/KpWidth ) )
!
         !MachScaling = min(1.0, abs(0.5*(MLeft_P12(iAlt,iSpecies) + &
         !                               MRight_P12(iAlt,iSpecies))))
         !MachScaling = 1.0
         !LiouKpS_P12(iAlt,iSpecies) = &
         !      (MachScaling**2.0)*exp(LogKpS_P12(iAlt,iSpecies))
         !LiouKpS_P12(iAlt,iSpecies) = &
         !      (MachScaling)*exp(LogKpS_P12(iAlt,iSpecies))

         !MachScaling = min(1.0, abs(0.5*(MLeft_M12(iAlt,iSpecies) + &
         !                               MRight_M12(iAlt,iSpecies))))

         !MachScaling = 1.0
         !LiouKpS_M12(iAlt,iSpecies) = &
         !      (MachScaling**2.0)*exp(LogKpS_M12(iAlt,iSpecies))
!         LiouKpS_M12(iAlt,iSpecies) = &
!               (MachScaling**2.0)*exp(LogKpS_M12(iAlt,iSpecies))

         LiouKpS_P12(iAlt,iSpecies) = exp(LogKpS_P12(iAlt,iSpecies))
         LiouKpS_M12(iAlt,iSpecies) = exp(LogKpS_M12(iAlt,iSpecies))

         LiouKuS_P12(iAlt,iSpecies) = exp(LogKuS_P12(iAlt,iSpecies))
         LiouKuS_M12(iAlt,iSpecies) = exp(LogKuS_M12(iAlt,iSpecies))

     enddo 
  enddo 

!!!! Grab the left and right states of the Variables
!!!!  on boh Interfaces (P12 = +1/2 and M12 = -1/2)
    do iSpecies = 1, nSpecies
          !! Calculate the Left and Right Faces of the RhoS 
           call calc_kt_facevalues(LogRhoS(-1:nAlts+2,iSpecies), &
                           LogRhoSLeft_M12( 1:nAlts  ,iSpecies), &
                          LogRhoSRight_M12( 1:nAlts  ,iSpecies), &
                           LogRhoSLeft_P12( 1:nAlts  ,iSpecies), &
                          LogRhoSRight_P12( 1:nAlts  ,iSpecies) )

           RhoSLeft_M12(:,iSpecies) = exp( LogRhoSLeft_M12(:,iSpecies)) 
          RhoSRight_M12(:,iSpecies) = exp(LogRhoSRight_M12(:,iSpecies)) 

           RhoSLeft_P12(:,iSpecies) = exp( LogRhoSLeft_P12(:,iSpecies)) 
          RhoSRight_P12(:,iSpecies) = exp(LogRhoSRight_P12(:,iSpecies)) 

    enddo 

!    write(*,*) '==================== CALC_HYDRO_FLUXES ================'
!    write(*,*) '======= M12 FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, RhoS(i-1), Left_M12, Right_M12, RhoS(i) = ', &
!                 iAlt, RhoS(iAlt-1,iSpecies), RhoSLeft_M12(iAlt,iSpecies), &
!          RhoSRight_M12(iAlt,iSpecies), RhoS(iAlt,iSpecies)  
!    enddo 
!    enddo 
!
!    write(*,*) '======= P12 FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, RhoS(i-1), Left_P12, Right_P12, RhoS(i) = ', &
!                 iAlt, RhoS(iAlt,iSpecies), RhoSLeft_P12(iAlt,iSpecies), &
!          RhoSRight_P12(iAlt,iSpecies), RhoS(iAlt+1,iSpecies)  
!    enddo 
!    enddo 
!
!stop


   do iSpecies = 1, nSpecies
         !! Calculate the Left and Right Faces of the PressureS
      call calc_kt_facevalues(LogPS(:,iSpecies), &
          LogPressureSLeft_M12(:,iSpecies), LogPressureSRight_M12(:,iSpecies), &
          LogPressureSLeft_P12(:,iSpecies), LogPressureSRight_P12(:,iSpecies) )

      PressureSLeft_M12(:,iSpecies) = exp( LogPressureSLeft_M12(:,iSpecies)) 
      PressureSRight_M12(:,iSpecies) = exp(LogPressureSRight_M12(:,iSpecies)) 

      PressureSLeft_P12(:,iSpecies) = exp( LogPressureSLeft_P12(:,iSpecies)) 
      PressureSRight_P12(:,iSpecies) = exp(LogPressureSRight_P12(:,iSpecies)) 

   enddo 

   do iSpecies = 1, nSpecies
         !! Calculate the Left and Right Faces of the Var (Rho) 
     call calc_kt_facevalues(VertVel(:,iSpecies), &
                         VelLeft_M12(:,iSpecies), VelRight_M12(:,iSpecies), &
                         VelLeft_P12(:,iSpecies), VelRight_P12(:,iSpecies) )
   enddo 

!    write(*,*) '======= M12 VELOCITY FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, Vel(i-1), VelLeft_M12, VelRight_M12, Vel(i) = ', &
!                 iAlt, VertVel(iAlt-1,iSpecies), VelLeft_M12(iAlt,iSpecies), &
!          VelRight_M12(iAlt,iSpecies), VertVel(iAlt,iSpecies)  
!    enddo 
!    enddo 
!
!    write(*,*) '======= P12 VELOCITY FACEVALUES ========'
!    do iSpecies = 1, nSpecies
!       write(*,*) '  SPECIES TYPE :    ',iSpecies, cSpecies(iSpecies)
!    do iAlt = 1, nAlts
!        write(*,*) 'iAlt, Vel(i-1), VelLeft_P12, VelRight_P12, Vel(i) = ', &
!                 iAlt, VertVel(iAlt-1,iSpecies), VelLeft_P12(iAlt,iSpecies), &
!          VelRight_P12(iAlt,iSpecies), VertVel(iAlt,iSpecies)  
!    enddo 
!    enddo 


   RhoLeft_M12(:) = 0.0
   RhoRight_M12(:) = 0.0

   RhoLeft_P12(:) = 0.0
   RhoRight_P12(:) = 0.0

   do iSpecies = 1, nSpecies
      RhoLeft_M12(1:nAlts) = RhoLeft_M12(1:nAlts) + &
                            RhoSLeft_M12(1:nAlts,iSpecies)
      RhoRight_M12(1:nAlts) = RhoRight_M12(1:nAlts) + &
                             RhoSRight_M12(1:nAlts,iSpecies)

      RhoLeft_P12(1:nAlts) = RhoLeft_P12(1:nAlts) + &
                            RhoSLeft_P12(1:nAlts,iSpecies)
      RhoRight_P12(1:nAlts) = RhoRight_P12(1:nAlts) + &
                             RhoSRight_P12(1:nAlts,iSpecies)
   enddo 

   PLeft_M12(:) = 0.0
   PRight_M12(:) = 0.0
 
   PLeft_P12(:) = 0.0
   PRight_P12(:) = 0.0

   do iSpecies = 1, nSpecies
      PLeft_M12(1:nAlts) = PLeft_M12(1:nAlts) + &
                   PressureSLeft_M12(1:nAlts,iSpecies)
      PRight_M12(1:nAlts) = PRight_M12(1:nAlts) + &
                    PressureSRight_M12(1:nAlts,iSpecies)

      PLeft_P12(1:nAlts) = PLeft_P12(1:nAlts) + &
                   PressureSLeft_P12(1:nAlts,iSpecies)
      PRight_P12(1:nAlts) = PRight_P12(1:nAlts) + &
                    PressureSRight_P12(1:nAlts,iSpecies)
   enddo 


   do iDim = 1, 3
      !! Calculate the Left and Right Faces of the Var (Rho) 
       call calc_kt_facevalues(Vel_GD(:,iDim), &
                        VelGDLeft_M12(:,iDim), &
                       VelGDRight_M12(:,iDim), &
                     VelGDLeft_P12(:,iDim), &
                       VelGDRight_P12(:,iDim) )
   enddo 

   VelGDLeft_M12(:,iUp_) = 0.0
   VelGDRight_M12(:,iUp_) = 0.0

   VelGDLeft_P12(:,iUp_) = 0.0
   VelGDRight_P12(:,iUp_) = 0.0

   do iSpecies = 1, nSpecies

        VelGDLeft_M12(1:nAlts,iUp_) = &
                 VelGDLeft_M12(1:nAlts,iUp_) + &
                 RhoSLeft_M12(1:nAlts,iSpecies)*&
                 VelLeft_M12(1:nAlts,iSpecies)/RhoLeft_M12(1:nAlts)

        VelGDRight_M12(1:nAlts,iUp_) = &
                 VelGDRight_M12(1:nAlts,iUp_) + &
                 RhoSRight_M12(1:nAlts,iSpecies)*&
                 VelRight_M12(1:nAlts,iSpecies)/RhoRight_M12(1:nAlts)

        VelGDLeft_P12(1:nAlts,iUp_) = &
                 VelGDLeft_P12(1:nAlts,iUp_) + &
                 RhoSLeft_P12(1:nAlts,iSpecies)*&
                 VelLeft_P12(1:nAlts,iSpecies)/RhoLeft_P12(1:nAlts)

        VelGDRight_P12(1:nAlts,iUp_) = &
                 VelGDRight_P12(1:nAlts,iUp_) + &
                 RhoSRight_P12(1:nAlts,iSpecies)*&
                 VelRight_P12(1:nAlts,iSpecies)/RhoRight_P12(1:nAlts)
   enddo 

        !! Calculate the Left and Right Faces of the Pressure 
   call calc_kt_facevalues(Gamma_1d(:), GammaLeft_M12(:), GammaRight_M12(:), &
                                        GammaLeft_P12(:), GammaRight_P12(:) )

    do iAlt = 1, nAlts
!!!! ============= Bulk Values
        ELeft_M12(iAlt) = &
            ( 1.0/(GammaLeft_M12(iAlt) - 1.0))*PLeft_M12(iAlt) + &
              0.5*RhoLeft_M12(iAlt)*&
             (VelGDLeft_M12(iAlt,iUp_)**2.0 + VelGDLeft_M12(iAlt,iEast_)**2.0 + &
              VelGDLeft_M12(iAlt,iNorth_)**2.0)

        ERight_M12(iAlt) = &
            ( 1.0/(GammaRight_M12(iAlt) - 1.0))*PRight_M12(iAlt) + &
              0.5*RhoRight_M12(iAlt)*&
             (VelGDRight_M12(iAlt,iUp_)**2.0 + VelGDRight_M12(iAlt,iEast_)**2.0 + &
              VelGDRight_M12(iAlt,iNorth_)**2.0)

        ELeft_P12(iAlt) = &
            ( 1.0/(GammaLeft_P12(iAlt) - 1.0))*PLeft_P12(iAlt) + &
              0.5*RhoLeft_P12(iAlt)* &
             (VelGDLeft_P12(iAlt,iUp_)**2.0 + VelGDLeft_P12(iAlt,iEast_)**2.0 + &
              VelGDLeft_P12(iAlt,iNorth_)**2.0)

        ERight_P12(iAlt) = &
            ( 1.0/(GammaRight_P12(iAlt) - 1.0))*PRight_P12(iAlt) + &
              0.5*RhoRight_P12(iAlt)* &
             (VelGDRight_P12(iAlt,iUp_)**2.0 + VelGDRight_P12(iAlt,iEast_)**2.0 + &
              VelGDRight_P12(iAlt,iNorth_)**2.0)

    enddo 



!!!! Liou et al. [2006] suggest using the Enthalpy for the 
!!!! numerical speed of sound.
!!!! Calculate the Enthalpy at the cell faces here.
!    write(*,*) 'ENTHALPY CALC:=============='
    do iAlt = 1, nAlts 

!      write(*,*) 'Gamma(i-1), GammaLeft_M12, GammaRight_M12, Gamma(i) =', &
!            Gamma_1d(iAlt-1), GammaLeft_M12(iAlt), GammaRight_M12(iAlt), &
!            Gamma_1d(iAlt)

       LiouEnthalpyLeft_M12(iAlt) = &
         0.5*(VelGDLeft_M12(iAlt,iUp_)**2.0 + &
              VelGDLeft_M12(iAlt,iEast_)**2.0 + &
              VelGDLeft_M12(iAlt,iNorth_)**2.0) + &
            (GammaLeft_M12(iAlt)/(GammaLeft_M12(iAlt) - 1.0))*&
            PLeft_M12(iAlt)/RhoLeft_M12(iAlt)

       LiouEnthalpyRight_M12(iAlt) = &
            0.5*(VelGDRight_M12(iAlt,iUp_)**2.0 + &
                 VelGDRight_M12(iAlt,iEast_)**2.0 + &
                 VelGDRight_M12(iAlt,iNorth_)**2.0) + &
            (GammaRight_M12(iAlt)/(GammaRight_M12(iAlt) - 1.0))*&
            PRight_M12(iAlt)/RhoRight_M12(iAlt)


       LiouEnthalpyLeft_P12(iAlt) = &
            0.5*(VelGDLeft_P12(iAlt,iUp_)**2.0 + &
                 VelGDLeft_P12(iAlt,iEast_)**2.0 + &
                 VelGDLeft_P12(iAlt,iNorth_)**2.0) + &
            (GammaLeft_P12(iAlt)/(GammaLeft_P12(iAlt) - 1.0))*&
            PLeft_P12(iAlt)/RhoLeft_P12(iAlt)


       LiouEnthalpyRight_P12(iAlt) = &
            0.5*(VelGDRight_P12(iAlt,iUp_)**2.0 + &
                 VelGDRight_P12(iAlt,iEast_)**2.0 + &
                 VelGDRight_P12(iAlt,iNorth_)**2.0) + &
            (GammaRight_P12(iAlt)/(GammaRight_P12(iAlt) - 1.0))*&
            PRight_P12(iAlt)/RhoRight_P12(iAlt)

!      write(*,*) 'Gamma(i), GammaLeft_P12, GammaRight_P12, Gamma(i+1) =', &
!            Gamma_1d(iAlt), GammaLeft_P12(iAlt), GammaRight_P12(iAlt), &
!            Gamma_1d(iAlt+1)

!      write(*,*) ' PLeft_M12(iAlt), PRight_M12, PLeft_P12, PRight_P12 =', &
!            PLeft_M12(iAlt), PRight_M12(iAlt), PLeft_P12(iAlt), &
!            PRight_P12(iAlt)
    enddo 

!!! Liou Methodology

   do iAlt = 1, nAlts

       SubCs = & 
         sqrt(2.0*( (GammaLeft_M12(iAlt) - 1.0 )/(GammaLeft_M12(iAlt) + 1.0)) *&
         LiouEnthalpyLeft_M12(iAlt) )

         LiouCSLeft_M12(iAlt) = & 
             (SubCs**2.0)/max(SubCs, VelGDLeft_M12(iAlt,iUp_))

       SubCs = &
       sqrt(2.0*( (GammaRight_M12(iAlt) - 1.0 )/(GammaRight_M12(iAlt) + 1.0)) *&
        LiouEnthalpyRight_M12(iAlt) )

       LiouCSRight_M12(iAlt) = &
           (SubCs**2.0)/max(SubCs, -1.0*VelGDRight_M12(iAlt,iUp_))

       InterfaceCS_M12(iAlt) = min(LiouCSLeft_M12(iAlt), LiouCSRight_M12(iAlt))

      SubCs = &
         sqrt(2.0*( (GammaLeft_P12(iAlt) - 1.0 )/(GammaLeft_P12(iAlt) + 1.0)) *&
            LiouEnthalpyLeft_P12(iAlt) )
      LiouCSLeft_P12(iAlt) = &
            (SubCs**2.0)/max(SubCs, VelGDLeft_P12(iAlt,iUp_))

      SubCs = &
      sqrt(2.0*( (GammaRight_P12(iAlt) - 1.0 )/(GammaRight_P12(iAlt) + 1.0)) *&
         LiouEnthalpyRight_P12(iAlt) )

         LiouCSRight_P12(iAlt) = &
              (SubCs**2.0)/max(SubCs, -1.0*VelGDRight_P12(iAlt,iUp_))

      InterfaceCS_P12(iAlt) = min(LiouCSLeft_P12(iAlt), LiouCSRight_P12(iAlt))


!      write(*,*) 'iAlt, InterfaceCS_M12(iAlt), InterfaceCS_P12(iAlt) =',&
!                  iAlt, InterfaceCS_M12(iAlt), InterfaceCS_P12(iAlt)
   enddo 

   ! Next, we modify the Velocities (only velocities) normal to the interface
   ! According to the local Mach Number (by species)

   do iSpecies = 1, nSpecies
      do iAlt = 1, nAlts
       !real :: ZVar, VL, CL, VR, CR
         CBar = InterfaceCS_P12(iAlt)
         VL =  VelLeft_P12(iAlt,iSpecies)
         VR = VelRight_P12(iAlt,iSpecies)
         ML = VL/CBar
         MR = VR/CBar
         ZVar = min(1.0, max(ML,MR))
         ! New Left and Right States
          VelLeft_P12(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VL - VR)
         VelRight_P12(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VR - VL)

         CBar = InterfaceCS_M12(iAlt)
         VL =  VelLeft_M12(iAlt,iSpecies)
         VR = VelRight_M12(iAlt,iSpecies)
         ML = VL/CBar
         MR = VR/CBar
         ZVar = min(1.0, max(ML,MR))
         ! New Left and Right States
          VelLeft_M12(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VL - VR)
         VelRight_M12(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VR - VL)

         
      enddo 
   enddo 


!stop
  
  MeanPressureS_P12(1:nAlts,1:nSpecies) = &
        0.5*(PressureSLeft_P12(1:nAlts,1:nSpecies) + &
            PressureSRight_P12(1:nAlts,1:nSpecies))

  MeanPressureS_M12(1:nAlts,1:nSpecies) = &
        0.5*(PressureSLeft_M12(1:nAlts,1:nSpecies) + &
            PressureSRight_M12(1:nAlts,1:nSpecies))

  MeanRhoS_P12(1:nAlts,1:nSpecies) = &
       0.5*(RhoSLeft_P12(1:nAlts,1:nSpecies) + RhoSRight_P12(1:nAlts,1:nSpecies))

  MeanRhoS_M12(1:nAlts,1:nSpecies) = &
       0.5*(RhoSLeft_M12(1:nAlts,1:nSpecies) + RhoSRight_M12(1:nAlts,1:nSpecies))

  do iAlt = 1, nAlts 
     MeanCS_P12(iAlt) = InterfaceCS_P12(iAlt)
     MeanCS_M12(iAlt) = InterfaceCS_M12(iAlt)
  enddo 

!!! Next, define local mach numbers at the interfaces

   do iAlt = 1, nAlts 
     do iSpecies = 1, nSpecies

        MLeft_M12(iAlt,iSpecies) = &
           VelLeft_M12(iAlt,iSpecies)/MeanCS_M12(iAlt)

        MRight_M12(iAlt,iSpecies) = &
           VelRight_M12(iAlt,iSpecies)/MeanCS_M12(iAlt)

        MLeft_P12(iAlt,iSpecies) = &
           VelLeft_P12(iAlt,iSpecies)/MeanCS_P12(iAlt)

        MRight_P12(iAlt,iSpecies) = &
           VelRight_P12(iAlt,iSpecies)/MeanCS_P12(iAlt)
    
     enddo 
   enddo 

     do iSpecies = 1, nSpecies
       M2Bar_M12(1:nAlts,iSpecies) = &
        0.5*(MLeft_M12(1:nAlts,iSpecies)**2.0 + &
             MRight_M12(1:nAlts,iSpecies)**2.0 )

        M2Bar_P12(1:nAlts,iSpecies) = &
           0.5*(MLeft_P12(1:nAlts,iSpecies)**2.0 + &
               MRight_P12(1:nAlts,iSpecies)**2.0 )
     enddo 

     do iSpecies = 1, nSpecies
      do iAlt = 1, nAlts 

        M2Zero_M12(iAlt,iSpecies) = &
              min(1.0, max(M2Bar_M12(iAlt,iSpecies), MInf)) 
        M2Zero_P12(iAlt,iSpecies) = &
              min(1.0, max(M2Bar_P12(iAlt,iSpecies), MInf)) 

        MZero_M12(iAlt,iSpecies) = sqrt(M2Zero_M12(iAlt,iSpecies))
        MZero_P12(iAlt,iSpecies) = sqrt(M2Zero_P12(iAlt,iSpecies))

      enddo 
     enddo 

     do iSpecies = 1, nSpecies
      do iAlt = 1, nAlts 

        FA_M12(iAlt,iSpecies) = &
               MZero_M12(iAlt,iSpecies)*(2.0 - MZero_M12(iAlt,iSpecies))
        FA_P12(iAlt,iSpecies) = &
               MZero_P12(iAlt,iSpecies)*(2.0 - MZero_P12(iAlt,iSpecies))

      enddo 
     enddo 

   do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts

       MF1P_Left_M12(iAlt,iSpecies) = &
             0.5*(MLeft_M12(iAlt,iSpecies) + abs(MLeft_M12(iAlt,iSpecies)) )
       MF1N_Left_M12(iAlt,iSpecies) = &
             0.5*(MLeft_M12(iAlt,iSpecies) - abs(MLeft_M12(iAlt,iSpecies)) )

       MF1P_Right_M12(iAlt,iSpecies) = &
             0.5*(MRight_M12(iAlt,iSpecies) + abs(MRight_M12(iAlt,iSpecies)) )
       MF1N_Right_M12(iAlt,iSpecies) = &
             0.5*(MRight_M12(iAlt,iSpecies) - abs(MRight_M12(iAlt,iSpecies)) )

       MF1P_Left_P12(iAlt,iSpecies) = &
             0.5*(MLeft_P12(iAlt,iSpecies) + abs(MLeft_P12(iAlt,iSpecies)) )
       MF1N_Left_P12(iAlt,iSpecies) = &
             0.5*(MLeft_P12(iAlt,iSpecies) - abs(MLeft_P12(iAlt,iSpecies)) )

       MF1P_Right_P12(iAlt,iSpecies) = &
             0.5*(MRight_P12(iAlt,iSpecies) + abs(MRight_P12(iAlt,iSpecies)) )
       MF1N_Right_P12(iAlt,iSpecies) = &
             0.5*(MRight_P12(iAlt,iSpecies) - abs(MRight_P12(iAlt,iSpecies)) )
     enddo 
  enddo 

   do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts

      MF2P_Left_M12(iAlt,iSpecies) =  0.25*(MLeft_M12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Left_M12(iAlt,iSpecies) = -0.25*(MLeft_M12(iAlt,iSpecies) - 1.0)**2.0

      MF2P_Right_M12(iAlt,iSpecies) =  0.25*(MRight_M12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Right_M12(iAlt,iSpecies) = -0.25*(MRight_M12(iAlt,iSpecies) - 1.0)**2.0

      MF2P_Left_P12(iAlt,iSpecies) =  0.25*(MLeft_P12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Left_P12(iAlt,iSpecies) = -0.25*(MLeft_P12(iAlt,iSpecies) - 1.0)**2.0

      MF2P_Right_P12(iAlt,iSpecies) =  0.25*(MRight_P12(iAlt,iSpecies) + 1.0)**2.0
      MF2N_Right_P12(iAlt,iSpecies) = -0.25*(MRight_P12(iAlt,iSpecies) - 1.0)**2.0

     enddo 
  enddo 


!!! Begin 4th Order Mach Number Polynomials
  do iSpecies = 1, nSpecies
    do iAlt = 1, nAlts 
      if ( abs(MLeft_M12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4P_Left_M12(iAlt,iSpecies) = MF1P_Left_M12(iAlt,iSpecies)
      else
         MF4P_Left_M12(iAlt,iSpecies) = MF2P_Left_M12(iAlt,iSpecies)*&
                   (1.0 - 16.0*LiouBeta*MF2N_Left_M12(iAlt,iSpecies))
      endif 

      if ( abs(MRight_M12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4N_Right_M12(iAlt,iSpecies) = MF1N_Right_M12(iAlt,iSpecies)
      else
         MF4N_Right_M12(iAlt,iSpecies) = MF2N_Right_M12(iAlt,iSpecies)*&
                    (1.0 + 16.0*LiouBeta*MF2P_Right_M12(iAlt,iSpecies))
      endif 

      if ( abs(MLeft_P12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4P_Left_P12(iAlt,iSpecies) = MF1P_Left_P12(iAlt,iSpecies)
      else
         MF4P_Left_P12(iAlt,iSpecies) = MF2P_Left_P12(iAlt,iSpecies)*&
                   (1.0 - 16.0*LiouBeta*MF2N_Left_P12(iAlt,iSpecies))
      endif 

      if ( abs(MRight_P12(iAlt,iSpecies)) .ge. 1.0) then 
         MF4N_Right_P12(iAlt,iSpecies) = MF1N_Right_P12(iAlt,iSpecies)
      else
         MF4N_Right_P12(iAlt,iSpecies) = MF2N_Right_P12(iAlt,iSpecies)*&
                    (1.0 + 16.0*LiouBeta*MF2P_Right_P12(iAlt,iSpecies))
      endif 


    enddo 
  enddo 

!!! Begin 2nd Order Mach Number Polynomials
   do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts 
 
        ModifiedZeta(iAlt,iSpecies) = 1.0
 
        MPress_P12(iAlt,iSpecies) = &
           LiouKpS_P12(iAlt,iSpecies)*max( (1.0 - M2Bar_P12(iAlt,iSpecies)), 0.0)*&
           (PressureSRight_P12(iAlt, iSpecies) - PressureSLeft_P12(iAlt,iSpecies) )/&
           ( MeanRhoS_P12(iAlt,iSpecies)*MeanCS_P12(iAlt)*&
           (FA_P12(iAlt,iSpecies)*MeanCS_P12(iAlt) + &
           ModifiedZeta(iAlt,iSpecies)*dAlt_F(iAlt+1)/DtIn)  ) 
 
        MPress_M12(iAlt,iSpecies) = &
          LiouKpS_M12(iAlt,iSpecies)*max( (1.0 - M2Bar_M12(iAlt,iSpecies)), 0.0)*&
          (PressureSRight_M12(iAlt, iSpecies) - PressureSLeft_M12(iAlt,iSpecies) )/&
          ( MeanRhoS_M12(iAlt,iSpecies)*MeanCS_M12(iAlt)*&
          (FA_M12(iAlt,iSpecies)*MeanCS_M12(iAlt) + &
           ModifiedZeta(iAlt,iSpecies)*dAlt_F(iAlt)/DtIn)  ) 
 
      enddo ! iSpecies
    enddo ! iAlt


   do iAlt = 1, nAlts 
      do iSpecies = 1, nSpecies

          InterfaceMach_M12(iAlt,iSpecies) =  &
                 MF4P_Left_M12(iAlt,iSpecies) + MF4N_Right_M12(iAlt,iSpecies) &
                    - MPress_M12(iAlt,iSpecies)

          LiouNumericalVelocity_M12(iAlt,iSpecies) = &
                   MeanCS_M12(iAlt)*InterfaceMach_M12(iAlt,iSpecies) 

          InterfaceMach_P12(iAlt,iSpecies) =  &
                 MF4P_Left_P12(iAlt,iSpecies) + MF4N_Right_P12(iAlt,iSpecies) &
                    - MPress_P12(iAlt,iSpecies)
 
          LiouNumericalVelocity_P12(iAlt,iSpecies) = &
                   MeanCS_P12(iAlt)*InterfaceMach_P12(iAlt,iSpecies) 
      enddo ! iSpecies 
   enddo  ! iAlt

      NumericalVelocity_M12(1:nAlts,1:nSpecies) = &
  LiouNumericalVelocity_M12(1:nAlts,1:nSpecies)

      NumericalVelocity_P12(1:nAlts,1:nSpecies) = &
  LiouNumericalVelocity_P12(1:nAlts,1:nSpecies)

!! NUMERICAL PRESSURE
  do iSpecies = 1, nSpecies
       do iAlt = 1, nAlts

             NumericalPressure_P12(iAlt,iSpecies) = &
               (MeanPressureS_P12(iAlt,iSpecies) )&
                - 0.5*LiouKuS_P12(iAlt,iSpecies)*MeanCS_P12(iAlt)*&
               ( (RhoSRight_P12(iAlt,iSpecies) )*VelRight_P12(iAlt,iSpecies) - &
                  (RhoSLeft_P12(iAlt,iSpecies) )* VelLeft_P12(iAlt,iSpecies) )

             NumericalPressure_M12(iAlt,iSpecies) = &
               (MeanPressureS_M12(iAlt,iSpecies) )&
                - 0.5*LiouKuS_M12(iAlt,iSpecies)*MeanCS_M12(iAlt)*&
               ( (RhoSRight_M12(iAlt,iSpecies) )*VelRight_M12(iAlt,iSpecies) - &
                  (RhoSLeft_M12(iAlt,iSpecies) )* VelLeft_M12(iAlt,iSpecies) )
           
       enddo !!iAlt = 1, nAlts
  enddo !!! nSpecies


  do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts 
        if (NumericalVelocity_P12(iAlt,iSpecies) .ge. 0.0) then 
           RhoSFlux_P12(iAlt,iSpecies) = &
               (RhoSLeft_P12(iAlt,iSpecies)*NumericalVelocity_P12(iAlt,iSpecies)) 

           MomentumSFlux_P12(iAlt,iSpecies) = &
               (RhoSLeft_P12(iAlt,iSpecies)*VelLeft_P12(iAlt,iSpecies)*&
               NumericalVelocity_P12(iAlt,iSpecies)) 
        else
           RhoSFlux_P12(iAlt,iSpecies) = &
               (RhoSRight_P12(iAlt,iSpecies)*NumericalVelocity_P12(iAlt,iSpecies)) 

           MomentumSFlux_P12(iAlt,iSpecies) = &
               (RhoSRight_P12(iAlt,iSpecies)*VelRight_P12(iAlt,iSpecies)*&
               NumericalVelocity_P12(iAlt,iSpecies)) 
        endif
        if (NumericalVelocity_M12(iAlt,iSpecies) .ge. 0.0) then 
           RhoSFlux_M12(iAlt,iSpecies) = &
               (RhoSLeft_M12(iAlt,iSpecies)*NumericalVelocity_M12(iAlt,iSpecies)) 

           MomentumSFlux_M12(iAlt,iSpecies) = &
               (RhoSLeft_M12(iAlt,iSpecies)*VelLeft_M12(iAlt,iSpecies)*&
               NumericalVelocity_M12(iAlt,iSpecies)) 
        else
           RhoSFlux_M12(iAlt,iSpecies) = &
               (RhoSRight_M12(iAlt,iSpecies)*NumericalVelocity_M12(iAlt,iSpecies)) 

           MomentumSFlux_M12(iAlt,iSpecies) = &
               (RhoSRight_M12(iAlt,iSpecies)*VelRight_M12(iAlt,iSpecies)*&
               NumericalVelocity_M12(iAlt,iSpecies)) 
        endif
     enddo 
  enddo 


    BulkNumericalVelocity_P12(1:nAlts) = 0.0
    BulkNumericalVelocity_M12(1:nAlts) = 0.0

    do iSpecies = 1, nSpecies
        BulkNumericalVelocity_P12(1:nAlts) = &
        BulkNumericalVelocity_P12(1:nAlts) + &
           ( RhoSLeft_P12(1:nAlts,iSpecies) + RhoSRight_P12(1:nAlts,iSpecies) )*&
             NumericalVelocity_P12(1:nAlts,iSpecies)/&
             (RhoLeft_P12(1:nAlts) + RhoRight_P12(1:nAlts) )

        BulkNumericalVelocity_M12(1:nAlts) = &
        BulkNumericalVelocity_M12(1:nAlts) + &
           ( RhoSLeft_M12(1:nAlts,iSpecies) + RhoSRight_M12(1:nAlts,iSpecies) )*&
             NumericalVelocity_M12(1:nAlts,iSpecies)/&
             (RhoLeft_M12(1:nAlts) + RhoRight_M12(1:nAlts) )
    enddo 


    do iAlt = 1, nAlts 
       if (BulkNumericalVelocity_P12(iAlt) .ge. 0.0) then 
           EnergyFlux_P12(iAlt) = &
             ( ELeft_P12(iAlt) + PLeft_P12(iAlt) ) &
               *BulkNumericalVelocity_P12(iAlt) 

           Momentum_P12(iAlt,1:3) = &
                RhoLeft_P12(iAlt)*VelGDLeft_P12(iAlt,1:3)*&
                BulkNumericalVelocity_P12(iAlt) 
       else
           EnergyFlux_P12(iAlt) = &
             ( ERight_P12(iAlt) + PRight_P12(iAlt) ) &
               *BulkNumericalVelocity_P12(iAlt) 

           Momentum_P12(iAlt,1:3) = &
                RhoRight_P12(iAlt)*VelGDRight_P12(iAlt,1:3)*&
                BulkNumericalVelocity_P12(iAlt) 
       endif
       if (BulkNumericalVelocity_M12(iAlt) .ge. 0.0) then 
           EnergyFlux_M12(iAlt) = &
             ( ELeft_M12(iAlt) + PLeft_M12(iAlt) ) &
               *BulkNumericalVelocity_M12(iAlt) 

           Momentum_M12(iAlt,1:3) = &
                RhoLeft_M12(iAlt)*VelGDLeft_M12(iAlt,1:3)*&
                BulkNumericalVelocity_M12(iAlt) 
       else
           EnergyFlux_M12(iAlt) = &
             ( ERight_M12(iAlt) + PRight_M12(iAlt) ) &
               *BulkNumericalVelocity_M12(iAlt) 

           Momentum_M12(iAlt,1:3) = &
                RhoRight_M12(iAlt)*VelGDRight_M12(iAlt,1:3)*&
                BulkNumericalVelocity_M12(iAlt) 
       endif
    enddo !! End iAlt Loop 

    do iAlt = 1, nAlts
          LeftRadius(iAlt) = 0.5*(RadDist(iAlt) + RadDist(iAlt-1))
          RightRadius(iAlt) = 0.5*(RadDist(iAlt) + RadDist(iAlt+1))
    enddo 

    do iAlt = 1, nAlts
          !AreaFunction_P12(iAlt) = RightRadius(iAlt)**2.0
          !AreaFunction_M12(iAlt) = LeftRadius(iAlt)**2.0
          !LocalCellVolume(iAlt) = &
          !   (1.0/3.0)*(RightRadius(iAlt)**3.0 - LeftRadius(iAlt)**3.0)

          AreaFunction_P12(iAlt) = Area_P12(iAlt)
          AreaFunction_M12(iAlt) = Area_M12(iAlt)
          LocalCellVolume(iAlt) = CellVol1D(iAlt)
    enddo 

    do iSpecies = 1, nSpecies
       RhoSFlux(1:nAlts,iSpecies) = &
            ( (AreaFunction_P12(1:nAlts))*RhoSFlux_P12(1:nAlts,iSpecies) - &
              (AreaFunction_M12(1:nAlts))*RhoSFlux_M12(1:nAlts,iSpecies) )/&
              LocalCellVolume(1:nAlts)

       MomentumSFlux(1:nAlts,iSpecies) = &
            ( (AreaFunction_P12(1:nAlts))*MomentumSFlux_P12(1:nAlts,iSpecies) - &
              (AreaFunction_M12(1:nAlts))*MomentumSFlux_M12(1:nAlts,iSpecies) )/&
              LocalCellVolume(1:nAlts)

       MomentumSFlux(1:nAlts,iSpecies) = &
       MomentumSFlux(1:nAlts,iSpecies) + &
           (NumericalPressure_P12(1:nAlts,iSpecies) - &
            NumericalPressure_M12(1:nAlts,iSpecies))/dAlt_C(1:nAlts)
    enddo 

    EnergyFlux(1:nAlts) = &
         ( (AreaFunction_P12(1:nAlts))*EnergyFlux_P12(1:nAlts) - &
           (AreaFunction_M12(1:nAlts))*EnergyFlux_M12(1:nAlts))/&
           LocalCellVolume(1:nAlts)

     do iDim = 1, 3
       MomentumFlux(1:nAlts,iDim) = &
            ( AreaFunction_P12(1:nAlts)*Momentum_P12(1:nAlts,iDim) - &
              AreaFunction_M12(1:nAlts)*Momentum_M12(1:nAlts,iDim))/&
              LocalCellVolume(1:nAlts)
     enddo 

end subroutine calc_all_fluxes_hydro


 subroutine calc_kt_facevalues(Var, VarLeft_M12, VarRight_M12, &
                                    VarLeft_P12, VarRight_P12)
 
   use ModVertical, only: &
       dAlt_F, dAlt_C, InvDAlt_F, Altitude_G, &
       Mesh_ULP120, Mesh_ULP121, Mesh_ULP122, &
       Mesh_CLP120, Mesh_CLP121, Mesh_CLP122, &
       Mesh_URM120, Mesh_URM121, Mesh_URM122, &
       Mesh_CRM120, Mesh_CRM121, Mesh_CRM122, &
       UB_MeshCoefs, LB_MeshCoefs, &
       Mesh_IS0, Mesh_IS1, Mesh_IS2
     
   use ModSizeGITM, only: nAlts
   use ModLimiterGitm
 
   implicit none
   
   real, intent(in) :: Var(-1:nAlts+2)
   real, intent(out):: VarLeft_P12(1:nAlts), VarRight_P12(1:nAlts)
   real, intent(out):: VarLeft_M12(1:nAlts), VarRight_M12(1:nAlts)
 
   real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)
   real :: rp, rn, rp2, rn2
   real :: Rn3, Rp3, R
   real :: Meshkm2, Meshkm1, Meshk, Meshkp1, Meshkp2

   real :: UL_P120, UL_P121, UL_P122
   real :: UR_M120, UR_M121, UR_M122
   real :: CL_P120, CL_P121, CL_P122
   real :: CR_M120, CR_M121, CR_M122
   real :: X_P12, X_P32, X_P52
   real :: X_M12, X_M32, X_M52
   real :: IS0, IS1, IS2
   real :: IS0Z, IS1Z, IS2Z
   real :: TRis
   real :: Alpha_P0, Alpha_P1, Alpha_P2
   real :: Alpha_M0, Alpha_M1, Alpha_M2
   real :: W_P0, W_P1, W_P2
   real :: W_M0, W_M1, W_M2
   integer :: i, k,iAlt
   real :: RISk
   real :: WENOEpsilon
   !------------ USE THESE FOR THE FAST VERSION --------------
   real :: C_P120, C_P121, C_P122
   real :: C_M120, C_M121, C_M122
   real :: U_P12, U_M12
   real :: U_P32, U_M32
   real :: U_BarP12, U_BarM12
   real :: U_BarP32, U_BarM32
   real :: QL_P120, QL_P121, QL_P122
   real :: QR_M120, QR_M121, QR_M122
   real :: hi, hip1
   real :: Alpha_M12, Alpha_M32, Alpha_P12, Alpha_P32
   !-------
   real :: DenominatorL, DenominatorR

!!! Use for 4-th Order Forward Differences
!!! Need a 5-point Stencil
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4

!!! Use for 4-th Order Backward Differences
!!! Need a 5-point Stencil
  real :: hm1, hm2, hm3, hm4
  real :: MeshHm1, MeshHm2, MeshHm3, MeshHm4
  real :: MeshCoefm0, MeshCoefm1, &
          MeshCoefm2, MeshCoefm3, &
          MeshCoefm4
  real :: dVar
  real :: LocalVar(-2:nAlts+3)  ! Need this for extrapolation
  real :: LocalVarLeft(0:nAlts), LocalVarRight(0:nAlts)
  real :: LocalVarLeft_M12(1:nAlts), LocalVarRight_M12(1:nAlts)
  real :: LocalVarLeft_P12(1:nAlts), LocalVarRight_P12(1:nAlts)
  real :: Alt_M2, Alt_M3
  real :: Alt_P2, Alt_P3
  ! WENOZ  - Factors
  real :: Tau5Z


  LocalVar(-1:nAlts+2) = Var(-1:nAlts+2)
   ! Use WENO Reconstruction
   !WENOEpsilon = 1.0e-6
   WENOEpsilon = 1.0e-2
   do iAlt=1,nAlts  
      UL_P120 = &
        Mesh_ULP120(iAlt,1)*   & 
                LocalVar(iAlt+1)  + &
        Mesh_ULP120(iAlt,2)*   &
                (LocalVar(iAlt) - LocalVar(iAlt+1)) - &
        Mesh_ULP120(iAlt,3)*   &
                (LocalVar(iAlt+2) - LocalVar(iAlt+1))
!          
      UL_P121 = &
        Mesh_ULP121(iAlt,1)*   & 
                   LocalVar(iAlt  ) + &
        Mesh_ULP121(iAlt,2)*   &
               (LocalVar(iAlt+1) - LocalVar(iAlt  )) - &
        Mesh_ULP121(iAlt,3)*   &
                (LocalVar(iAlt-1) - LocalVar(iAlt  )) 

      UL_P122 = &
        Mesh_ULP122(iAlt,1)*   & 
                   LocalVar(iAlt-1) + &
        Mesh_ULP122(iAlt,2)*   &
               (LocalVar(iAlt-2) - LocalVar(iAlt-1)) + &
        Mesh_ULP122(iAlt,3)*   &
               (LocalVar(iAlt  ) - LocalVar(iAlt-1)) 
       
      UR_M120 = &
        Mesh_URM120(iAlt,1)*   & 
                   LocalVar(iAlt+1) + &
        Mesh_URM120(iAlt,2)*   &
                (LocalVar(iAlt  ) - LocalVar(iAlt+1)) + & 
        Mesh_URM120(iAlt,3)*   &
                (LocalVar(iAlt+2) - LocalVar(iAlt+1)) 

      UR_M121 = &
        Mesh_URM121(iAlt,1)*   & 
                   LocalVar(iAlt  ) + &
        Mesh_URM121(iAlt,2)*   &
                  (LocalVar(iAlt-1) - LocalVar(iAlt  )) - &
        Mesh_URM121(iAlt,3)*   &
                  (LocalVar(iAlt+1) - LocalVar(iAlt  )) 
       
      UR_M122 = &
        Mesh_URM122(iAlt,1)*   & 
                   LocalVar(iAlt-1) + &
        Mesh_URM122(iAlt,2)*   &
                  (LocalVar(iAlt) - LocalVar(iAlt-1)) - &
        Mesh_URM122(iAlt,3)*   &
                  (LocalVar(iAlt-2) - LocalVar(iAlt-1)) 

      CL_P120 = Mesh_CLP120(iAlt)
      CL_P121 = Mesh_CLP121(iAlt)
      CL_P122 = Mesh_CLP122(iAlt)

      CR_M120 = Mesh_CRM120(iAlt)
      CR_M121 = Mesh_CRM121(iAlt)
      CR_M122 = Mesh_CRM122(iAlt)

      IS0 = Mesh_IS0(iAlt,1)*&
            (LocalVar(iAlt+2) - LocalVar(iAlt+1))**2.0 + &
            Mesh_IS0(iAlt,2)*&
            ((LocalVar(iAlt+2) - LocalVar(iAlt+1))*(LocalVar(iAlt) - LocalVar(iAlt+1))) + &
            Mesh_IS0(iAlt,3)*&
            ((LocalVar(iAlt  ) - LocalVar(iAlt+1))**2.0)
            
      IS1 = Mesh_IS1(iAlt,1)*&
            (LocalVar(iAlt-1) - LocalVar(iAlt  ))**2.0 + &
            Mesh_IS1(iAlt,2)*&
            ((LocalVar(iAlt+1) - LocalVar(iAlt  ))*(LocalVar(iAlt-1) - LocalVar(iAlt))) + &
            Mesh_IS1(iAlt,3)*&
            ((LocalVar(iAlt+1) - LocalVar(iAlt  ))**2.0)

      IS2 = Mesh_IS2(iAlt,1)*&
            (LocalVar(iAlt-2) - LocalVar(iAlt-1))**2.0 + &
            Mesh_IS2(iAlt,2)*&
            ((LocalVar(iAlt  ) - LocalVar(iAlt-1))*(LocalVar(iAlt-2) - LocalVar(iAlt-1))) + &
            Mesh_IS2(iAlt,3)*&
            ((LocalVar(iAlt  ) - LocalVar(iAlt-1))**2.0)

      Tau5Z = abs(IS0 - IS2)

      Alpha_P0 = CL_P120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
      Alpha_P1 = CL_P121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
      Alpha_P2 = CL_P122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

      Alpha_M0 = CR_M120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
      Alpha_M1 = CR_M121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
      Alpha_M2 = CR_M122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

      DenominatorL = Alpha_P0 + Alpha_P1 + Alpha_P2
      DenominatorR = Alpha_M0 + Alpha_M1 + Alpha_M2

      W_P0 = Alpha_P0/DenominatorL
      W_P1 = Alpha_P1/DenominatorL
      W_P2 = Alpha_P2/DenominatorL

      W_M0 = Alpha_M0/DenominatorR
      W_M1 = Alpha_M1/DenominatorR
      W_M2 = Alpha_M2/DenominatorR

       LocalVarLeft_P12(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122
      LocalVarRight_M12(iAlt  ) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122
   enddo !i=1,nAlts

   do iAlt = 2, nAlts
      LocalVarLeft_M12(iAlt) = LocalVarLeft_P12(iAlt-1)
   enddo !iAlt = 1, nAlts
   do iAlt = 1, nAlts-1
      LocalVarRight_P12(iAlt) = LocalVarRight_M12(iAlt+1)
   enddo !iAlt = 1, nAlts

   iAlt = -2
   dVar  = LB_MeshCoefs(1,1)*Var(iAlt+1) + &  
           LB_MeshCoefs(1,2)*Var(iAlt+2) + &  
           LB_MeshCoefs(1,3)*Var(iAlt+3) + &  
           LB_MeshCoefs(1,4)*Var(iAlt+4) + &  
           LB_MeshCoefs(1,5)*Var(iAlt+5)      
   LocalVar(iAlt) = Var(iAlt+1) - dAlt_F(iAlt+1)*dVar 


   iAlt = nAlts + 3
   dVar  = UB_MeshCoefs(3,1)*Var(iAlt-1) + &  
           UB_MeshCoefs(3,2)*Var(iAlt-2) + &  
           UB_MeshCoefs(3,3)*Var(iAlt-3) + &  
           UB_MeshCoefs(3,4)*Var(iAlt-4) + &  
           UB_MeshCoefs(3,5)*Var(iAlt-5)      
   LocalVar(iAlt) = LocalVar(iAlt-1) + dVar*dAlt_F(iAlt-1)


   ! Now, the extra ghost cells (-2, and nAlts + 3) are filled, so we can use
   ! The WENOZ-RL(5) for the VarLeft(0), VarRight(nAlts)
   iAlt = 0
   UL_P120 = &
        Mesh_ULP120(iAlt,1)*   & 
                LocalVar(iAlt+1)  + &
        Mesh_ULP120(iAlt,2)*   &
                (LocalVar(iAlt) - LocalVar(iAlt+1)) - &
        Mesh_ULP120(iAlt,3)*   &
                (LocalVar(iAlt+2) - LocalVar(iAlt+1))
!          
   UL_P121 = &
        Mesh_ULP121(iAlt,1)*   & 
                   LocalVar(iAlt  ) + &
        Mesh_ULP121(iAlt,2)*   &
               (LocalVar(iAlt+1) - LocalVar(iAlt  )) - &
        Mesh_ULP121(iAlt,3)*   &
                (LocalVar(iAlt-1) - LocalVar(iAlt  )) 

   UL_P122 = &
        Mesh_ULP122(iAlt,1)*   & 
                   LocalVar(iAlt-1) + &
        Mesh_ULP122(iAlt,2)*   &
               (LocalVar(iAlt-2) - LocalVar(iAlt-1)) + &
        Mesh_ULP122(iAlt,3)*   &
               (LocalVar(iAlt  ) - LocalVar(iAlt-1)) 
       
   CL_P120 = Mesh_CLP120(iAlt)
   CL_P121 = Mesh_CLP121(iAlt)
   CL_P122 = Mesh_CLP122(iAlt)

   IS0 = Mesh_IS0(iAlt,1)*&
         (LocalVar(iAlt+2) - LocalVar(iAlt+1))**2.0 + &
         Mesh_IS0(iAlt,2)*&
         ((LocalVar(iAlt+2) - LocalVar(iAlt+1))*(LocalVar(iAlt) - LocalVar(iAlt+1))) + &
         Mesh_IS0(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt+1))**2.0)
            
   IS1 = Mesh_IS1(iAlt,1)*&
         (LocalVar(iAlt-1) - LocalVar(iAlt  ))**2.0 + &
         Mesh_IS1(iAlt,2)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))*(LocalVar(iAlt-1) - LocalVar(iAlt))) + &
         Mesh_IS1(iAlt,3)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))**2.0)

   IS2 = Mesh_IS2(iAlt,1)*&
         (LocalVar(iAlt-2) - LocalVar(iAlt-1))**2.0 + &
         Mesh_IS2(iAlt,2)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))*(LocalVar(iAlt-2) - LocalVar(iAlt-1))) + &
         Mesh_IS2(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))**2.0)

   Tau5Z = abs(IS0 - IS2)

   Alpha_P0 = CL_P120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
   Alpha_P1 = CL_P121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
   Alpha_P2 = CL_P122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

   DenominatorL = Alpha_P0 + Alpha_P1 + Alpha_P2

   W_P0 = Alpha_P0/DenominatorL
   W_P1 = Alpha_P1/DenominatorL
   W_P2 = Alpha_P2/DenominatorL

   LocalVarLeft_M12(iAlt+1) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122


   iAlt = nAlts+1
   UR_M120 = &
     Mesh_URM120(iAlt,1)*   & 
                LocalVar(iAlt+1) + &
     Mesh_URM120(iAlt,2)*   &
             (LocalVar(iAlt  ) - LocalVar(iAlt+1)) + & 
     Mesh_URM120(iAlt,3)*   &
             (LocalVar(iAlt+2) - LocalVar(iAlt+1)) 

   UR_M121 = &
     Mesh_URM121(iAlt,1)*   & 
                LocalVar(iAlt  ) + &
     Mesh_URM121(iAlt,2)*   &
               (LocalVar(iAlt-1) - LocalVar(iAlt  )) - &
     Mesh_URM121(iAlt,3)*   &
               (LocalVar(iAlt+1) - LocalVar(iAlt  )) 
       
   UR_M122 = &
     Mesh_URM122(iAlt,1)*   & 
                LocalVar(iAlt-1) + &
     Mesh_URM122(iAlt,2)*   &
               (LocalVar(iAlt) - LocalVar(iAlt-1)) - &
     Mesh_URM122(iAlt,3)*   &
               (LocalVar(iAlt-2) - LocalVar(iAlt-1)) 

   CR_M120 = Mesh_CRM120(iAlt)
   CR_M121 = Mesh_CRM121(iAlt)
   CR_M122 = Mesh_CRM122(iAlt)

   IS0 = Mesh_IS0(iAlt,1)*&
         (LocalVar(iAlt+2) - LocalVar(iAlt+1))**2.0 + &
         Mesh_IS0(iAlt,2)*&
         ((LocalVar(iAlt+2) - LocalVar(iAlt+1))*(LocalVar(iAlt) - LocalVar(iAlt+1))) + &
         Mesh_IS0(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt+1))**2.0)
            
   IS1 = Mesh_IS1(iAlt,1)*&
         (LocalVar(iAlt-1) - LocalVar(iAlt  ))**2.0 + &
         Mesh_IS1(iAlt,2)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))*(LocalVar(iAlt-1) - LocalVar(iAlt))) + &
         Mesh_IS1(iAlt,3)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))**2.0)

   IS2 = Mesh_IS2(iAlt,1)*&
         (LocalVar(iAlt-2) - LocalVar(iAlt-1))**2.0 + &
         Mesh_IS2(iAlt,2)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))*(LocalVar(iAlt-2) - LocalVar(iAlt-1))) + &
         Mesh_IS2(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))**2.0)

   Tau5Z = abs(IS0 - IS2)

   Alpha_M0 = CR_M120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
   Alpha_M1 = CR_M121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
   Alpha_M2 = CR_M122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

   DenominatorR = Alpha_M0 + Alpha_M1 + Alpha_M2

   W_M0 = Alpha_M0/DenominatorR
   W_M1 = Alpha_M1/DenominatorR
   W_M2 = Alpha_M2/DenominatorR

   LocalVarRight_P12(iAlt-1) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122

   do iAlt = 1, nAlts
       VarLeft_P12(iAlt) = LocalVarLeft_P12(iAlt)
      VarRight_P12(iAlt) = LocalVarRight_P12(iAlt)

       VarLeft_M12(iAlt) =  LocalVarLeft_M12(iAlt)
      VarRight_M12(iAlt) = LocalVarRight_M12(iAlt)
   enddo 

 end subroutine calc_kt_facevalues


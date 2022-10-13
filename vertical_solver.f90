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
            DivRhoSFlux, DivMomentumSFlux, &
            DivEnergyFlux, DivMomentumFlux, RadDist) 

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax, VertVel, Gamma_1D, &
                          Vel_GD, MeanMajorMass_1d, &
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
  real, intent(out):: DivRhoSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: DivMomentumSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: DivEnergyFlux(1:nAlts)
  real, intent(out):: DivMomentumFlux(1:nAlts,1:3)
  real, intent(in) :: RadDist(-1:nAlts+2)

  integer :: iSpecies, iAlt, iDim
  real :: SubCs

  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogRhoS
  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogPS

  ! Local Fluxes through the Faces
  real, dimension(0:nAlts,1:nSpecies) :: RhoSFlux
  real, dimension(0:nAlts,1:nSpecies) :: MomentumSFlux
  real, dimension(0:nAlts) :: EnergyFlux
  real, dimension(0:nAlts,3) :: MomentumFlux
  real, dimension(0:nAlts) :: AreaFunction

  real, dimension( 1:nAlts) :: LocalCellVolume
  real:: MInf, LiouBeta

!

  real :: MachScaling
  real :: ZVar, VL, VR, CBar, ML, MR

  ! Updated Kp/Ku Variables
  ! 0:nAlts:  Also, no Species component
  real :: KpWidth, KpAltMidPoint
  integer :: AltIndex
  real :: LiouKpS(0:nAlts,1:nSpecies), LiouKuS(0:nAlts,1:nSpecies)
  real :: Kp, Ku

  !                  Pressure, VerticalVelocity
  real, dimension(0:nAlts,1:nSpecies) :: LogRhoSLeft, LogRhoSRight
  real, dimension(0:nAlts,1:nSpecies) :: RhoSLeft, RhoSRight

  real, dimension(0:nAlts,1:nSpecies) :: LogPressureSLeft, LogPressureSRight
  real, dimension(0:nAlts,1:nSpecies) :: PressureSLeft, PressureSRight
  real, dimension(0:nAlts,1:nSpecies) :: VelLeft, VelRight

  ! Mean Variables @ Interface (mean of Left and Right Vars)
  real, dimension(0:nAlts,1:nSpecies) :: MeanPressureS
  real, dimension(0:nAlts,1:nSpecies) :: MeanRhoS
  ! Vertical Velocity Variables
  ! Bulk Rho Variables: Rho, Pressure, Velocity, Gamma, Energy
  real, dimension(0:nAlts) :: RhoLeft, RhoRight
  real, dimension(0:nAlts) :: PLeft, PRight
  real, dimension(0:nAlts,3) :: VelGDLeft, VelGDRight
  real, dimension(0:nAlts) :: GammaLeft, GammaRight
  real, dimension(0:nAlts) :: ELeft, ERight

  ! AUSM-solver specific Variables
  ! Full Enthalpy Variable at Interface needed for num. speed of sound
  real, dimension(0:nAlts) :: LiouEnthalpyLeft, LiouEnthalpyRight
  real, dimension(0:nAlts) :: LiouCSLeft, LiouCSRight

  ! Mach Number Variables
  real, dimension(0:nAlts,1:nSpecies) :: MLeft, MRight
  real, dimension(0:nAlts) :: MeanCS
  real, dimension(0:nAlts,1:nSpecies) :: M2Bar
  real, dimension(0:nAlts,1:nSpecies) :: M2Zero
  real, dimension(0:nAlts,1:nSpecies) :: MZero
  real, dimension(0:nAlts,1:nSpecies) :: FA

  ! Mach Number Polynomials 
  ! First Order Polynomial
  real, dimension(0:nAlts,1:nSpecies) :: MF1P_Left
  real, dimension(0:nAlts,1:nSpecies) :: MF1N_Right
  ! 2nd Order Polynomial
  real, dimension(0:nAlts,1:nSpecies) :: MF2P_Left, MF2N_Left
  real, dimension(0:nAlts,1:nSpecies) :: MF2P_Right, MF2N_Right
  ! 4th Order Polynomial
  real, dimension(0:nAlts,1:nSpecies) :: MF4P_Left, MF4N_Left
  real, dimension(0:nAlts,1:nSpecies) :: MF4P_Right, MF4N_Right
  ! MPress: Numerical Pressure
  real, dimension(0:nAlts,1:nSpecies) :: MPress
  ! Entropy Correction Factor
  real, dimension(0:nAlts,1:nSpecies):: ModifiedZeta
  ! Interface Mach Number
  real, dimension(0:nAlts,1:nSpecies) :: InterfaceMach
  ! Numerical Velocity
  real, dimension(0:nAlts,1:nSpecies) :: NumericalVelocity
  ! AUSM Numerical Pressure 
  real, dimension(0:nAlts,1:nSpecies) :: NumericalPressure
  ! Numerical Velocity
  real, dimension(0:nAlts) :: BulkNumericalVelocity

!  ! BEGIN SETTING UP VARIABLES
  ! Basic Geometry of the Grid
  !   |------------------|--------- X ---------|------------------|
  !   |-----(i-1)--------|--------- i ---------|------(i+1)-------|
  !                   (i-1/2)               (i+1/2)
  !                    L | R                 L | R
  !            UL_M12(i) | UR_M12(i)  UL_P12(i) | UR_P12(i)
  !                     i ranges from 1 -> nAlts
  !
  ! Current Variable Setup focuses only on the +1/2 Values from 0 to nAlts
  !         |           UL(i-1)| UR(i-1)       UL(i)| UR(i)
  !                           i ranges from 0 -> nAlts
  ! U(i)  |---0---|---- 1 ------|  ...(i-1) -- | -- (i) -- | ... (nAlts) -- | -- (nAlts+1) 
  ! UI(i) |--   UI(0) -----   UI(1) ...      UI(i-1)         ...        UI(nAlts) 
  ! L/R   | [UL(0)|UR(0)] [UL(1)|UR(1)] [UL(i-1)|UR(i-1)]    ...  [UL(nAlts)|UR(nAlts)]

  ! Set global AUSM Variables: MInf is our free stream Mach Number
  ! Not sure what to specify, but chose MInf = 1.0e-19
  MInf = 1.0e-19
  LiouBeta = 1.0/8.0

  ! Take Logarithm of Exponentially Varying Quantities
  LogRhoS(-1:nAlts+2,1:nSpecies) = alog(     RhoS(-1:nAlts+2,1:nSpecies))
  LogPS(  -1:nAlts+2,1:nSpecies) = alog(PressureS(-1:nAlts+2,1:nSpecies))

  do iSpecies = 1, nSpecies
     !! Calculate the Left and Right Faces of the RhoS 
     call calc_kt_facevalues(LogRhoS(:,iSpecies), &
                       LogRhoSLeft( :  ,iSpecies), &
                      LogRhoSRight( :  ,iSpecies))
          RhoSLeft(:,iSpecies) = exp(LogRhoSLeft(:,iSpecies))
         RhoSRight(:,iSpecies) = exp(LogRhoSRight(:,iSpecies))
  enddo 

  do iSpecies = 1, nSpecies
     call calc_kt_facevalues(LogPS(-1:nAlts+2,iSpecies), &
                   LogPressureSLeft( :  ,iSpecies), &
                  LogPressureSRight( :  ,iSpecies))
          PressureSLeft(:,iSpecies) = exp(LogPressureSLeft( :,iSpecies))
         PressureSRight(:,iSpecies) = exp(LogPressureSRight(:,iSpecies))
  enddo 

  do iSpecies = 1, nSpecies
     call calc_kt_facevalues(VertVel(:,iSpecies), &
                                 VelLeft(:,iSpecies), &
                                VelRight(:,iSpecies))
  enddo 

  RhoLeft(:) = 0.0
  RhoRight(:) = 0.0
  do iSpecies = 1, nSpecies
     RhoLeft(: ) = RhoLeft( :) + RhoSLeft(:,iSpecies)
     RhoRight(:) = RhoRight(:) + RhoSRight(:,iSpecies)
  enddo 

  PLeft(:) = 0.0
  PRight(:) = 0.0
  do iSpecies = 1, nSpecies
     PLeft(: ) = PLeft(:)  + PressureSLeft(:,iSpecies)
     PRight(:) = PRight(:) + PressureSRight(:,iSpecies)
  enddo 

  do iDim = 1, 3
     call calc_kt_facevalues(Vel_GD(:,iDim), &
                               VelGDLeft(:,iDim), &
                             VelGDRight(:,iDim))
  enddo 

  VelGDLeft(:,iUp_) = 0.0
  VelGDRight(:,iUp_) = 0.0

  do iSpecies = 1, nSpecies
     VelGDLeft(:,iUp_) = VelGDLeft(:,iUp_) + &
           RhoSLeft(:,iSpecies)*VelLeft(:,iSpecies)/&
            RhoLeft(:)
 
     VelGDRight(:,iUp_) = VelGDRight(:,iUp_) + &
           RhoSRight(:,iSpecies)*VelRight(:,iSpecies)/&
            RhoRight(:)
  enddo 

  call calc_kt_facevalues(Gamma_1d(:), GammaLeft(:), GammaRight(:))
    ! Fill the old variables with our new ones

  do iAlt = 0, nAlts
      ELeft(iAlt) = &
          ( 1.0/(GammaLeft(iAlt) - 1.0))*PLeft(iAlt) + &
            0.5*RhoLeft(iAlt)* &
           (VelGDLeft(iAlt,iUp_)**2.0 + VelGDLeft(iAlt,iEast_)**2.0 + &
              VelGDLeft(iAlt,iNorth_)**2.0)
  
      ERight(iAlt) = &
            ( 1.0/(GammaRight(iAlt) - 1.0))*PRight(iAlt) + &
            0.5*RhoRight(iAlt)* &
             (VelGDRight(iAlt,iUp_)**2.0 + VelGDRight(iAlt,iEast_)**2.0 + &
              VelGDRight(iAlt,iNorth_)**2.0)
  enddo 

  do iAlt = 0, nAlts 
       LiouEnthalpyLeft(iAlt) = &
            0.5*(VelGDLeft(iAlt,iUp_)**2.0 + &
                 VelGDLeft(iAlt,iEast_)**2.0 + &
                 VelGDLeft(iAlt,iNorth_)**2.0) + &
            (GammaLeft(iAlt)/(GammaLeft(iAlt) - 1.0))*&
            PLeft(iAlt)/RhoLeft(iAlt)
  
 
       LiouEnthalpyRight(iAlt) = &
            0.5*(VelGDRight(iAlt,iUp_)**2.0 + &
                 VelGDRight(iAlt,iEast_)**2.0 + &
                 VelGDRight(iAlt,iNorth_)**2.0) + &
            (GammaRight(iAlt)/(GammaRight(iAlt) - 1.0))*&
            PRight(iAlt)/RhoRight(iAlt)
  enddo 

  do iAlt = 0, nAlts
      SubCs = sqrt(2.0*( (GammaLeft(iAlt) - 1.0 )/(GammaLeft(iAlt) + 1.0)) *&
           LiouEnthalpyLeft(iAlt) )
      LiouCSLeft(iAlt) = (SubCs**2.0)/max(SubCs, VelGDLeft(iAlt,iUp_))
 
      SubCs = sqrt(2.0*( (GammaRight(iAlt) - 1.0 )/(GammaRight(iAlt) + 1.0)) *&
           LiouEnthalpyRight(iAlt) )
      LiouCSRight(iAlt) = (SubCs**2.0)/max(SubCs, -1.0*VelGDRight(iAlt,iUp_))
      MeanCS(iAlt) = min(LiouCSLeft(iAlt), LiouCSRight(iAlt))
  enddo 

  do iSpecies = 1, nSpecies
     do iAlt = 0, nAlts
       !real :: ZVar, VL, CL, VR, CR
        CBar = MeanCS(iAlt)
        VL =  VelLeft(iAlt,iSpecies)
        VR = VelRight(iAlt,iSpecies)
        ML = VL/CBar
        MR = VR/CBar
        ZVar = min(1.0, max(ML,MR))
        ! New Left and Right States
         VelLeft(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VL - VR)
        VelRight(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VR - VL)
     enddo 
  enddo 

  
  MeanPressureS(:,1:nSpecies) = &
        0.5*(PressureSLeft(:,1:nSpecies) + &
            PressureSRight(:,1:nSpecies))
  MeanRhoS(:,1:nSpecies) = &
       0.5*(RhoSLeft(:,1:nSpecies) + RhoSRight(:,1:nSpecies))

  ! Local Mach Number Variables
  do iAlt = 0, nAlts 
     do iSpecies = 1, nSpecies
         MLeft(iAlt,iSpecies) = &
            VelLeft(iAlt,iSpecies)/MeanCS(iAlt)

         MRight(iAlt,iSpecies) = &
            VelRight(iAlt,iSpecies)/MeanCS(iAlt)

        M2Bar(iAlt,iSpecies) = &
           0.5*(MLeft(iAlt,iSpecies)**2.0 + &
               MRight(iAlt,iSpecies)**2.0 )

        M2Zero(iAlt,iSpecies) = &
              min(1.0, max(M2Bar(iAlt,iSpecies), MInf)) 
        MZero(iAlt,iSpecies) = sqrt(M2Zero(iAlt,iSpecies))

        FA(iAlt,iSpecies) = &
               MZero(iAlt,iSpecies)*(2.0 - MZero(iAlt,iSpecies))
       ! Fourth Order Mach Number Polynomial
       if ( abs(MLeft(iAlt,iSpecies)) .lt. 1.0) then 
          MF4P_Left(iAlt,iSpecies) = &
                         ( &
                                0.25*(MLeft(iAlt,iSpecies)      + 1.0)**2.0 + &
                            LiouBeta*(MLeft(iAlt,iSpecies)**2.0 - 1.0)**2.0 & 
                         )
       else
          MF4P_Left(iAlt,iSpecies) = &
             0.5*(MLeft(iAlt,iSpecies) + abs(MLeft(iAlt,iSpecies)) )
       endif 

       ! Fourth Order Mach Number Polynomial
       if ( abs(MRight(iAlt,iSpecies)) .lt. 1.0) then 
          MF4N_Right(iAlt,iSpecies) = &
                         ( &
                               -0.25*(MRight(iAlt,iSpecies)      - 1.0)**2.0 - &
                            LiouBeta*(MRight(iAlt,iSpecies)**2.0 - 1.0)**2.0 & 
                         )
       else
          MF4N_Right(iAlt,iSpecies) = &
             0.5*(MRight(iAlt,iSpecies) - abs(MRight(iAlt,iSpecies)) )
       endif 

     enddo 
  enddo 
  Kp = 0.10             !! Ullrich et al. [2011]
  Ku = 1.0             !! Ullrich et al. [2011]
  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
      ! Note that this enforces a pressure-type low-mach number variation
      MachScaling = MZero(iAlt,iSpecies)**3.0
      LiouKpS(iAlt,iSpecies) = (MachScaling)*Kp
      LiouKuS(iAlt,iSpecies) = (MachScaling)*Ku
     enddo 
  enddo 
  ! Numerical AUSM Fluxes
  do iSpecies = 1, nSpecies
     do iAlt = 0, nAlts 
       ModifiedZeta(iAlt,iSpecies) = 1.0
           MPress(iAlt,iSpecies) = &
              LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar(iAlt,iSpecies)), 0.0)*&
              (PressureSRight(iAlt, iSpecies) - PressureSLeft(iAlt,iSpecies) )/&
              ( MeanRhoS(iAlt,iSpecies)*MeanCS(iAlt)*&
              (FA(iAlt,iSpecies)*MeanCS(iAlt) + &
              ModifiedZeta(iAlt,iSpecies)*dAlt_F(iAlt+1)/DtIn)  ) 
      enddo ! iSpecies
    enddo ! iAlt

    do iSpecies = 1, nSpecies
        InterfaceMach(:,iSpecies) =  &
               MF4P_Left(:,iSpecies) + MF4N_Right(:,iSpecies) &
                - MPress(:,iSpecies)
        NumericalVelocity(:,iSpecies) = &
                 MeanCS(:)*InterfaceMach(:,iSpecies) 
    enddo ! iSpecies 

 !! NUMERICAL PRESSURE
   do iSpecies = 1, nSpecies
      NumericalPressure(:,iSpecies) = &
            (MeanPressureS(:,iSpecies) )&
              - 0.5*LiouKuS(:,iSpecies)*MeanCS(:)*&
              ( (RhoSRight(:,iSpecies) )*VelRight(:,iSpecies) - &
                 (RhoSLeft(:,iSpecies) )* VelLeft(:,iSpecies) )
   enddo !iSpecies = 1, nSpecies
   do iSpecies = 1, nSpecies
      do iAlt = 0, nAlts 
          ! Rewrite without the need for if statments
          RhoSFlux(iAlt,iSpecies) = &
                 0.5*RhoSLeft(iAlt,iSpecies)* &
                 (NumericalVelocity(iAlt,iSpecies) + & 
              abs(NumericalVelocity(iAlt,iSpecies) ) ) & 
                 +  &
                 0.5*RhoSRight(iAlt,iSpecies)* &
                 (NumericalVelocity(iAlt,iSpecies) - & 
              abs(NumericalVelocity(iAlt,iSpecies) ) ) 

             MomentumSFlux(iAlt,iSpecies) = &
                 0.5*RhoSLeft(iAlt,iSpecies)*VelLeft(iAlt,iSpecies)*&
                 (NumericalVelocity(iAlt,iSpecies) + & 
              abs(NumericalVelocity(iAlt,iSpecies) ) ) & 
                 +  &
                 0.5*RhoSRight(iAlt,iSpecies)*VelRight(iAlt,iSpecies)*&
                 (NumericalVelocity(iAlt,iSpecies) - & 
              abs(NumericalVelocity(iAlt,iSpecies) ) ) 
      enddo !iAlt = 0, nAlts 
   enddo !iSpecies = 1, nSpecies

   BulkNumericalVelocity(:) = 0.0
   do iSpecies = 1, nSpecies
      BulkNumericalVelocity(:) = &
      BulkNumericalVelocity(:) + &
         ( RhoSLeft(:,iSpecies) + RhoSRight(:,iSpecies))*&
            NumericalVelocity(:,iSpecies)/&
            (RhoLeft(:) + RhoRight(:) )
   enddo !iSpecies = 1, nSpecies
   do iAlt = 0, nAlts 
      EnergyFlux(iAlt) = &
         0.5*( ELeft(iAlt) + PLeft(iAlt) ) &
           *( BulkNumericalVelocity(iAlt)  + &
         abs( BulkNumericalVelocity(iAlt)) ) &
             + & 
         0.5*( ERight(iAlt) + PRight(iAlt) ) &
           *( BulkNumericalVelocity(iAlt)  - &
         abs( BulkNumericalVelocity(iAlt)) ) 

       MomentumFlux(iAlt,1:3) = &
            0.5*RhoLeft(iAlt)*VelGDLeft(iAlt,1:3)*&
            (   BulkNumericalVelocity(iAlt) + & 
            abs(BulkNumericalVelocity(iAlt)) ) &
             + &
            0.5*RhoRight(iAlt)*VelGDRight(iAlt,1:3)*&
            (   BulkNumericalVelocity(iAlt) - & 
            abs(BulkNumericalVelocity(iAlt)) ) 
   enddo !iAlt = 0, nAlts 

    do iAlt = 1, nAlts
          AreaFunction(iAlt) = Area_P12(iAlt)
          LocalCellVolume(iAlt) = CellVol1D(iAlt)
    enddo 
    ! The first Cell left hand side
    AreaFunction(0) = Area_M12(1)

    do iAlt = 1, nAlts
       do iSpecies = 1, nSpecies
         DivRhoSFlux(iAlt,iSpecies) = &
              ( (AreaFunction(iAlt  ))*RhoSFlux(iAlt  ,iSpecies) - &
                (AreaFunction(iAlt-1))*RhoSFlux(iAlt-1,iSpecies) )/&
                LocalCellVolume(iAlt)
         DivMomentumSFlux(iAlt,iSpecies) = &
              ( (AreaFunction(iAlt  ))*RhoSFlux(iAlt  ,iSpecies) - &
                (AreaFunction(iAlt-1))*RhoSFlux(iAlt-1,iSpecies) )/&
                LocalCellVolume(iAlt)
         DivMomentumSFlux(iAlt,iSpecies) = &
         DivMomentumSFlux(iAlt,iSpecies) + &
               (NumericalPressure(iAlt  ,iSpecies) - &
                NumericalPressure(iAlt-1,iSpecies))/dAlt_C(iAlt)
       enddo !iSpecies = 1, nSpecies
       DivEnergyFlux(iAlt) = &
            ( (AreaFunction(iAlt  ))*EnergyFlux(iAlt  ) - &
              (AreaFunction(iAlt-1))*EnergyFlux(iAlt-1))/&
              LocalCellVolume(iAlt)
       do iDim = 1, 3
         DivMomentumFlux(iAlt,iDim) = &
              ( AreaFunction(iAlt  )*MomentumFlux(iAlt  ,iDim) - &
                AreaFunction(iAlt-1)*MomentumFlux(iAlt-1,iDim))/&
              LocalCellVolume(iAlt)
       enddo !iDim

     enddo ! iAlt

end subroutine calc_all_fluxes_hydro


 subroutine calc_kt_facevalues(Var, VarLeft, VarRight)
 
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
   real, intent(out):: VarLeft(0:nAlts), VarRight(0:nAlts)
 
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
   !-------
   real :: DenominatorL, DenominatorR

!!! Use for 4-th Order Backward Differences
!!! Need a 5-point Stencil
  real :: dVar
  real :: LocalVar(-2:nAlts+3)  ! Need this for extrapolation
  real :: LocalVarLeft(0:nAlts), LocalVarRight(0:nAlts)
!  real :: LocalVarLeft_M12(0:nAlts+1), LocalVarRight_M12(0:nAlts+1)
!  real :: LocalVarLeft_P12(0:nAlts+1), LocalVarRight_P12(0:nAlts+1)
  ! WENOZ  - Factors
  real :: Tau5Z
  LocalVar(-1:nAlts+2) = Var(-1:nAlts+2)
  ! =\ 
  ! ==\
  ! Extend Local Vars upward and downward
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

      ! For a given "i" we calculate the inner faces of the cell
      ! That is we interpolate to the right-face of the -1/2 edge
      ! That is we interpolate to the left-face  of the +1/2 edge
      ! Let's define
      !     U(i)          = LocalVar(i)
      !     UL_P12/M12(i) = LocalVarLeft_P12(i)/LocalVarLeft_M12(i)
      !     UR_P12/M12(i) = LocalVarRight_P12(i)/LocalVarRight_M12(i)
      !      (i-1)        |               (i)                |       (i+1)
      !                   | UR(i)  <----- U(i) ----->  UL(i) |
      !  Old Notation
      !                   | UR_M12(i)  <--U(i) --> UL_P12(i) |
      !   Old Scheme (i) ranged from i = 1, nAlts
      !
      !  New Notation
      !                   | UR(i-1)    <--U(i) --> UL(i)     |
      !   New Scheme (i) ranges from i = 0, nAlts
      !   In New Scheme (Special Cases):  i = 0:  

!       LocalVarLeft_P12(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122
!      LocalVarRight_M12(iAlt  ) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122

       LocalVarLeft(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122
      LocalVarRight(iAlt-1) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122
   enddo !i=1,nAlts

   ! Special Cases
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

   LocalVarLeft(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122


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

   ! This isi LocalVarRight_P12(nAlts)
   LocalVarRight(iAlt-1) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122

   do iAlt = 0, nAlts
       VarLeft(iAlt) = LocalVarLeft(iAlt)
      VarRight(iAlt) = LocalVarRight(iAlt)
   enddo 

 end subroutine calc_kt_facevalues




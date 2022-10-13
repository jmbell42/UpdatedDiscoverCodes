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
  use ModInputs, only: UseBarriers, iDebugLevel
  implicit none
  !-----------------------------------------------------------
  ! JMB 07-10-2017
  ! The code below implements the classic explicit 
  ! Runge-Kutta 4th order solver.
  ! This is based upon Ullrich and Jablonowski [2012], Journal of Comp. Physics
  !  			"MCore: a non-hydrostatic atmospheric dynamical core ..."
  ! RK-4 Update and naming of the stages are taken from Equations (66) - (69)
  ! q(1) -> Stage1, the original state 
  ! q(2) -> Stage2 
  ! q(3) -> Stage3 
  ! q(4) -> Stage4 
  ! q(*) -> Final
  !
  ! q(2) = q(1) + (Dt/2.0)*H(q(1))
  ! q(3) = q(1) + (Dt/2.0)*H(q(2))
  ! q(4) = q(1) + (Dt    )*H(q(3))
  ! q(*) = (1.0/3.0)*( -q(1) + q(2) + 2*q(3) + q(4) + (Dt/2.0)*H(q(4)))
  !
  integer :: iError, iAlt, iSpecies, iDir
  ! These are the intermediate updated states.  This is overwritten several times
  real :: UpdatedLogNS(-1:nAlts+2,1:nSpecies)
  real :: UpdatedLogINS(-1:nAlts+2,1:nIons)
  real :: UpdatedLogRho(-1:nAlts+2)
  real :: UpdatedVel_GD(-1:nAlts+2,1:3)
  real :: UpdatedTemp(-1:nAlts+2)
  real :: UpdatedVS(-1:nAlts+2,1:nSpecies)
  ! These are the final RK-4 Updates to be passed out of vertical_solver.
  real :: FinalLogNS(-1:nAlts+2,1:nSpecies)
  real :: FinalLogINS(-1:nAlts+2,1:nIons)
  real :: FinalLogRho(-1:nAlts+2)
  real :: FinalVel_GD(-1:nAlts+2,1:3)
  real :: FinalTemp(-1:nAlts+2)
  real :: FinalVS(-1:nAlts+2,1:nSpecies)

!!! RStage-4 Coefficients
! For RK-4 We need 4 separate evaluations
! Hence Stage1 - Stage4
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

  ! We sub-divide the Original Dt into smaller segments
  ! Each segment is DtIn long, and is passed to the vertical_solver
  real :: DtIn

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 6) write(*,*) "=======> vertical bcs 1", iproc
  ! Step 1, Fill in Ghost Cells
  call set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp,LogINS,IVel,VertVel)

  ! Set q(1) -> Stage 1
  !
   Stage1LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage1LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage1LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage1Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage1Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage1VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel

  ! Calculate Dt/2.0*H(q(1))
  DtIn = 0.5*Dt  !!! 1st stage needs a 1/2 Dt step
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
  ! \
  ! Calculate q(2) 
  ! q(2) = q(1) + (Dt/2.0)*H(q(1))
       Stage2LogNS(-1:nAlts+2,1:nSpecies)      = &
       Stage1LogNS(-1:nAlts+2,1:nSpecies)      + &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      Stage2LogINS(-1:nAlts+2,1:nIons)  = &
      Stage1LogINS(-1:nAlts+2,1:nIons)  + &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      Stage2LogRho(-1:nAlts+2)                = &
      Stage1LogRho(-1:nAlts+2)                + &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      Stage2Vel_GD(-1:nAlts+2,1:3)            = & 
      Stage1Vel_GD(-1:nAlts+2,1:3)            + & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        Stage2Temp(-1:nAlts+2)                  = &
        Stage1Temp(-1:nAlts+2)                  + &
     (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          Stage2VS(-1:nAlts+2,1:nSpecies)      = &
          Stage1VS(-1:nAlts+2,1:nSpecies)      + &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

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

!!! Set Upper Boundary
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                        UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

   ! Fill the Variables properly and pass to _1stage
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

  ! Calculate (Dt/2.0)*H(q(2) )
  DtIn = 0.5*Dt
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! Calculate q(3)
  ! q(3) = q(1) + (Dt/2.0)*H(q(2))
       Stage3LogNS(-1:nAlts+2,1:nSpecies)      = &
       Stage1LogNS(-1:nAlts+2,1:nSpecies)      + &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      Stage3LogINS(-1:nAlts+2,1:nIons)  = &
      Stage1LogINS(-1:nAlts+2,1:nIons)  + &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      Stage3LogRho(-1:nAlts+2)                = &
      Stage1LogRho(-1:nAlts+2)                + &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      Stage3Vel_GD(-1:nAlts+2,1:3)            = & 
      Stage1Vel_GD(-1:nAlts+2,1:3)            + & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        Stage3Temp(-1:nAlts+2)                  = &
        Stage1Temp(-1:nAlts+2)                  + &
      (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          Stage3VS(-1:nAlts+2,1:nSpecies)      = &
          Stage1VS(-1:nAlts+2,1:nSpecies)      + &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

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
  ! Calculate Dt*H(q(3))
  DtIn = Dt
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! Calculate q(4) = q(1) + Dt*H(q(3))
       Stage4LogNS(-1:nAlts+2,1:nSpecies)      = &
       Stage1LogNS(-1:nAlts+2,1:nSpecies)      + &
     (NewLogNS(-1:nAlts+2,1:nSpecies) - LogNS(-1:nAlts+2,1:nSpecies))

      Stage4LogINS(-1:nAlts+2,1:nIons)  = &
      Stage1LogINS(-1:nAlts+2,1:nIons)  + &
    (NewLogINS(-1:nAlts+2,1:nIons) - LogINS(-1:nAlts+2,1:nIons))

      Stage4LogRho(-1:nAlts+2)                = &
      Stage1LogRho(-1:nAlts+2)                + &
    (NewLogRho(-1:nAlts+2) - LogRho(-1:nAlts+2))

      Stage4Vel_GD(-1:nAlts+2,1:3)            = & 
      Stage1Vel_GD(-1:nAlts+2,1:3)            + & 
    (NewVel_GD(-1:nAlts+2,1:3) - Vel_GD(-1:nAlts+2,1:3))

        Stage4Temp(-1:nAlts+2)                  = &
        Stage1Temp(-1:nAlts+2)                  + &
      (NewTemp(-1:nAlts+2) -  Temp(-1:nAlts+2))

          Stage4VS(-1:nAlts+2,1:nSpecies)      = &
          Stage1VS(-1:nAlts+2,1:nSpecies)      + &
   (NewVertVel(-1:nAlts+2,1:nSpecies) - VertVel(-1:nAlts+2,1:nSpecies))

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

  ! Update the Boundaries
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
  
!! Calculate (Dt/2.0)*H(q(4))
  DtIn = 0.5*Dt
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
  
  ! This section ensures that our lower boundary conditions are maintained
  ! and not overwritten.
  ! Note that Stage1 is our original state vector
  FinalLogNS = Stage1LogNS
  FinalLogINS = Stage1LogINS
  FinalLogRho = Stage1LogRho
  FinalVel_GD = Stage1Vel_GD
  FinalTemp   = Stage1Temp
  FinalVS     = Stage1VS

  ! Calculate the final RK-4 State q(*)
  ! q(*) = (1.0/3.0)*( -q(1) + q(2) + 2*q(3) + q(4) + (Dt/2.0)*H(q(4)))
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

  ! Final Boundary Condition Update
  call set_vertical_bcs(FinalLogRho, FinalLogNS, FinalVel_GD, &
                          FinalTemp, FinalLogINS, IVel, FinalVS)

  ! Pass updated fields back to advance_vertical.f90
   LogNS = FinalLogNS
  LogINS = FinalLogINS
  LogRho = FinalLogRho
  Vel_GD = FinalVel_GD
    Temp = FinalTemp
 VertVel = FinalVS

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d

!=============================================================================
subroutine advance_vertical_1stage_ausm( DtIn, &
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
       KappaTemp_1d, ViscCoefS_1d, &
       ChemSources_1d, &
       Centrifugal, &
       Gravity_G, Altitude_G,Cv_1D, dAlt_F
  use ModTime
  use ModInputs
  use ModConstants
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
  real, intent(inout) :: NewVertVel(-1:nAlts+2,nSpecies)

  real :: NS(-1:nAlts+2,nSpecies)
  real :: Rho(-1:nAlts+2)

  real, dimension(1:nAlts)    :: GradTemp, GradTempKoM, &
       DiffTemp, GradTmp, DiffTmp, DiviVel
  real, dimension(1:nAlts,3) :: GradiVel_CD, DiffiVel_CD

  real, dimension(1:nAlts,nIons) :: GradLogINS, DiffLogINS
  real :: NewSumRho, NewLogSumRho

  integer :: iAlt, iSpecies,  iDim

  real, dimension(-1:nAlts+2)    :: NT
  real, dimension(-1:nAlts+2)    :: Press
  real, dimension(1:nAlts,nSpecies)    :: EddyDiffusionVel

  real :: nVel(1:nAlts,1:nSpecies)
  integer :: nFilter, iFilter
  real :: LowFilter

  ! Parameters Used for the Sponge
  ! This Sponge is useful to dampen out spurious modes
  ! oscillating between the bottom and top of the model.
  integer :: nAltsSponge = 12
  real :: kSP, NuSP, AmpSP

  ! JMB:  Adding Eddy Diffusion Variables here
  ! Note:  These are used in the calc_neutral_friction
  !--------------------------------------------------------------------------
  !! Eddy Diffusion Variables
  real, dimension(1:nAlts,nSpecies)    :: GradLogConS
  real, dimension(-1:nAlts+2,nSpecies) :: ConS, LogConS
  real, dimension(1:nAlts,nSpecies)    :: EddyCoefRatio_1d
  !--------------------------------------------------------------------------
  ! 4th Order Gradients on a Non-Uniform Mesh (5-point Stencil)
  ! Used for calculating the d(ln[Chi])/dr -> Log of the concentration gradient
  !--------------------------------------------------------------------------
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4
  ! ----------------------------------------------------
  ! JMB:  AUSM Variables
  ! ----------------------------------------------------
  real ::    RhoS(-1:nAlts+2,1:nSpecies),&
          NewRhoS(-1:nAlts+2,1:nSpecies),&
      AUSMRhoSFluxes(1:nAlts,1:nSpecies)

  real :: InvScaleHeight, MeanGravity, MeanTemp, MeanMass

  real ::   HydroNS(-1:nAlts+2,1:nSpecies),&
          HydroRhoS(-1:nAlts+2,1:nSpecies), &
            HydroNT(-1:nAlts+2),&
     HydroPressureS(-1:nAlts+2,1:nSpecies),&
      DeviationRhoS(-1:nAlts+2,1:nSpecies),&
 DeviationRhoSRatio(-1:nAlts+2,1:nSpecies),&
      HydroPressure(-1:nAlts+2), &
           HydroRho(-1:nAlts+2)

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

  logical :: SubtractHydrostatic(-1:nAlts+2,1:nSpecies)   
  real :: RadialDistance_C(-1:nAlts+2)
  real :: EffectiveGravity(-1:nAlts+2)

  ! JMB:  Use these as Limiters on Winds for an initial startup
  real :: TimeFactor, Vel0, DeltaV, VelocityCap

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  
  ! Dynamic Velocity Cap.
  ! Vel0 is the initial velocity cap at time = 0.0
  ! The VelocityCap variable will vary from Vel0 -> MaximumVertVel
  ! Over time with an exponential timescale given by
  ! TimeFactor = exp(-tSimulation/alpha), alpha is set by the user
  Vel0 = 100.0 ! initial velocity in m/s
  TimeFactor = exp(-tSimulation/43200.0)
  DeltaV = MaximumVerticalVelocity - Vel0
  VelocityCap = Vel0 + DeltaV*(1.0 - TimeFactor)
  ! NOTE:  If you want VelocityCap to be MaximumVelocity, then just
  ! set Vel0 = MaximumVerticalVelocity and TimeFactor to 1.0

  ! Effective Gravity is the net gravity + Centrifugal term
  do iAlt = -1, nAlts + 2
     EffectiveGravity(iAlt) = &
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

  ! Pressure for the Total energy term
  do iAlt = -1, nAlts + 2
    Press(iAlt) = NT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
  enddo
  ! Grad Temp used by the thermal diffusion term
  ! retain the old method for doing this here
  call calc_rusanov_alts_ausm(Temp   ,GradTemp,    DiffTemp)

  ! Add geometrical correction to gradient and obtain divergence
  DiviVel = GradiVel_CD(:,iUp_) + 2*iVel(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)

  do iSpecies=1,nIons-1
     call calc_rusanov_alts_ausm(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo

  ! Step 1:  Calculate Ln(rho_{s}/Rho)
  ! You can use either rho_{s}/rho or n_{s}/n.
  ! rho_{s}/rho has the added benefit of including variations in 
  ! the mean major mass in the eddy diffusion term
  do iSpecies = 1, nSpecies
    LogConS(-1:nAlts+2,iSpecies) = &
         alog(Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)/Rho(-1:nAlts+2))
    !LogConS(-1:nAlts+2,iSpecies) = &
    !     alog(NS(-1:nAlts+2,iSpecies)/NT(-1:nAlts+2))
  enddo 

  ! Use 4th order non-uniform gradient
  ! Need a 5-point stencil
  do iAlt = 1, nAlts
     h1 = dAlt_F(iAlt-1)
     h2 = dAlt_F(iAlt+0)
     h3 = dAlt_F(iAlt+1)
     h4 = dAlt_F(iAlt+2)

     MeshH2 = h2 + h1
     MeshH3 = h3 + h2 + h1
     MeshH4 = h4 + h3 + h2 + h1

     MeshCoef0 = (h2*h3*(h3+h4))/(h1*MeshH2*MeshH3*MeshH4)
     MeshCoef1 = -1.0*(MeshH2*h3*(h3 + h4))/(h1*h2*(h2+h3)*(h2+h3+h4))
     MeshCoef3 = MeshH2*h2*(h4 + h3)/(MeshH3*(h2+h3)*h3*h4) 
     MeshCoef4 = -1.0*MeshH2*h2*h3/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)

     MeshCoef2 = (h2*h3*(h3+h4) + &
                  MeshH2*h3*(h3+h4) - &
                  MeshH2*h2*(h3+h4) - &
                  MeshH2*h2*h3)/&
                  (MeshH2*h2*h3*(h3+h4))

     do iSpecies = 1, nSpecies
        GradLogConS(iAlt,iSpecies) =  &
          MeshCoef0*LogConS(iAlt-2,iSpecies)&
       +  MeshCoef1*LogConS(iAlt-1,iSpecies)&
       +  MeshCoef2*LogConS(iAlt  ,iSpecies)&
       +  MeshCoef3*LogConS(iAlt+1,iSpecies)&
       +  MeshCoef4*LogConS(iAlt+2,iSpecies)
     enddo  ! iSpecies Loop
  enddo ! iAlt Loop

  ! JMB: 07/20/2017:
  ! AUSM SOLVER DESCRIPTION:
  !    For the AUSM solver as described in Ullrich and Jablonowski [2012],
  !    we need to define a background hydrostatic state
  !    to subtract off.
  !    In essence we need a term such that Grad(P_{s}) - p_{s}*g ~ 0.0
  !    The residual after the hydrostatic background is subtracted
  !    drives the non-hydrostatic dynamics
  ! 
  ! NOTE:  We subtract a hydrostatic state for each species
  ! because the fundamental force balance is between Grad(P_{S}) and
  ! Rho_{S}*g.  You could subtract Grad(P_{total}) - Rho_total, but this
  ! will not allow the species' momentum balance to be accurately
  ! calculated when the species is the major background species.   
  !
  ! Define the Hydrostatic State:  (1) HydroNS -each species
  !		                   (2) HydroNT --Total background
  HydroNS(-1:nAlts+2,1:nSpecies) = NS(-1:nAlts+2,1:nSpecies)
  HydroNT(-1:nAlts+2)            = NT(-1:nAlts+2)
  ! For Hydro NT, we assume NT(0) is fixed and HydroNT(0) must be the
  ! same value
  ! Next, we use the standard hydrostatic solution to integrate up.
  ! n(z) = n(0)*T(0)/T(z)*exp(-delta_z/Ha)
  do iAlt = 1, nAlts + 2
    MeanMass    = 0.5*(MeanMajorMass_1d(iAlt-1) + MeanMajorMass_1d(iAlt))
    MeanGravity = 0.5*( EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1) )
    MeanTemp    = 0.5*( Temp(iAlt) + Temp(iAlt-1) )
    InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                     (Boltzmanns_Constant*MeanTemp)
    HydroNT(iAlt) = HydroNT(iAlt-1)*(Temp(iAlt-1)/Temp(iAlt))*&
          exp (-1.0*dAlt_F(iAlt)*InvScaleHeight)
  enddo 
   ! Next, we integrate down to (-1) cell from the fixed NT(0)
  iAlt = -1
  MeanMass = 0.5*( MeanMajorMass_1d(-1) + MeanMajorMass_1d(0))
  MeanGravity = 0.5*( EffectiveGravity(iAlt+1) + EffectiveGravity(iAlt) )
  MeanTemp = 0.5*( Temp(iAlt+1) + Temp(iAlt) )
  InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                  (Boltzmanns_Constant*MeanTemp)
  HydroNT(iAlt) = HydroNT(iAlt+1)*(Temp(iAlt+1)/Temp(iAlt))*&
                  exp(dAlt_F(iAlt)*InvScaleHeight)

  ! Calculate a Background Atmosphere for each species
  ! Note: This really only works for the major species
  ! which is fine, since minor species will have large difference
  ! Grad(P) - rho*g terms and are not subject to the same level
  ! of sensitivity.
  !
  ! Apply same method as above.  Start from cell 0 and integrate up
  do iSpecies =  1, nSpecies
     do iAlt = 1, nAlts + 2
     MeanMass = Mass(iSpecies)
     MeanGravity = 0.5*( EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1) )
     MeanTemp = 0.5*( Temp(iAlt) + Temp(iAlt-1) )
     InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                      (Boltzmanns_Constant*MeanTemp)
     HydroNS(iAlt,iSpecies) = &
           HydroNS(iAlt-1,iSpecies)*(Temp(iAlt-1)/Temp(iAlt))*&
           exp (-1.0*dAlt_F(iAlt)*InvScaleHeight)
     enddo 
  enddo 
  ! Extend the densities downward with the mean mass of the atmosphere
  ! Note:  This assumes that the lower boundary is in the 
  ! turbulent regime.  This can be altered to the species mass
  ! if necessary.
  do iSpecies =  1, nSpecies
     iAlt = -1
     MeanMass = 0.5*( MeanMajorMass_1d(-1) + MeanMajorMass_1d(0))
     MeanGravity = 0.5*( EffectiveGravity(iAlt+1) + EffectiveGravity(iAlt) )
     MeanTemp = 0.5*( Temp(iAlt+1) + Temp(iAlt) )
     InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                      (Boltzmanns_Constant*MeanTemp)
     HydroNS(iAlt,iSpecies) = &
               HydroNS(iAlt+1,iSpecies)*(Temp(iAlt+1)/Temp(iAlt))*&
               exp(dAlt_F(iAlt)*InvScaleHeight)
  enddo 

  ! Calculate the Hydrostatic Pressure background.
  ! (1) for each species
  ! (2) for the bulk atmosphere
  do iAlt = -1, nAlts + 2
     do iSpecies =  1, nSpecies
       HydroPressureS(iAlt,iSpecies) = &
           HydroNS(iAlt,iSpecies)*Boltzmanns_Constant*Temp(iAlt)
       HydroRhoS(iAlt,iSpecies) = Mass(iSpecies)*HydroNS(iAlt,iSpecies)
     enddo 
       HydroPressure(iAlt) = HydroNT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
            HydroRho(iAlt) = HydroNT(iAlt)*MeanMajorMass_1d(iAlt)
  enddo 

  ! Initialize the AUSM+-up Variables
  !
  ! Main Prognostic variables are:
  ! (1) RhoS      -- m_{s}*n_{s}; Species Mass densities
  ! (2) MomentumS -- rho_{s}*V_{s; Species Vertical Momentum
  ! (3) Momentum  -- rho_{total}*U_{theta,phi}; Bulk Horizontal Momentum
  ! (4) Total Energy -- P/(gamma -1) + 0.5*rho*V^2.0

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
  ! Deviation Rho is used in the acceleration terms below:
  ! Quantifies the amount that rho_{s} and HydroRho_{s}
  ! Differs. 
  DeviationRhoS(-1:nAlts+2,1:nSpecies) = &
           RhoS(-1:nAlts+2, 1:nSpecies) - &
      HydroRhoS(-1:nAlts+2,1:nSpecies)

  ! Need to calculate the deviation from the local density
  ! from the background density.  
  ! Reason:  For some species that are not at all hydrostatic
  !          The AUSM+-up method of subtracting off the 
  !          hydrostatic background state fails miserably.
  !
  ! DeviationRhoSRatio is what we use to determine if
  ! we should subtract off the background state
  DeviationRhoSRatio(-1:nAlts+2,1:nSpecies) = & 
            abs(RhoS(-1:nAlts+2,1:nSpecies) - &
           HydroRhoS(-1:nAlts+2,1:nSpecies))/&
                RhoS(-1:nAlts+2,1:nSpecies)

  !
  ! Subtract Hydrostatic:  Logical Argument 
  !
  ! Logical check to see if we need to subtract off the
  ! background state.
  !
  ! Why Altitude Dependent?  Some species (i.e. CH4 at Titan) can
  ! exist in a nearly hydrostatic state deep in the atmosphere
  ! but get stripped upward by non-thermal processes, meaning that it's
  ! very much NOT hydrostatic.  This allows us to adapt the model
  ! to a situation where a species may be partially hydrostatic in 
  ! one regime, but not hydrostatic in another.
  do iAlt = -1, nAlts + 2
    do iSpecies = 1, nSpecies
      if ( abs(DeviationRhoSRatio(iAlt,iSpecies)) .gt. 1.0) then
          SubtractHydrostatic(iAlt,iSpecies) = .false.
      else
          SubtractHydrostatic(iAlt,iSpecies) = .true.
      endif 
    enddo 
  enddo 

  ! Initalize updated state variables.
  NewRho = Rho
  NewPress = Press
  NewTotalEnergy = TotalEnergy

  ! calc_all_fluxes_hydro takes the current state variables and returns
  ! the divergence of the fluxes for:
  ! (1) RhoS, MomentumS, Momentum(horizontal), and Energy
  call calc_all_fluxes_hydro(DtIn, RhoS, PressureS, HydroPressureS, HydroRhoS, &
        HydroPressure,  HydroRho, AUSMRhoSFluxes,AUSMMomentumSFluxes, &
        AUSMTotalEnergyFluxes, AUSMMomentumFluxes, RadialDistance_C,  &
        SubtractHydrostatic)

  ! Optional Sponge Layer Variables
  ! Taken from Dowling et al. [1999], EPIC Model description
  AmpSP = (1.0/(10.0*DtIn))
  kSP = nAltsSponge + 1
  do iAlt = 1,nAlts

     !Update the Mass Densities for each species
     do iSpecies=1,nSpecies
        NewRhoS(iAlt,iSpecies) = RhoS(iAlt,iSpecies) &
               - DtIn*(AUSMRhoSFluxes(iAlt,iSpecies))
        NewLogNS(iAlt,iSpecies) = alog( NewRhoS(iAlt,iSpecies)/Mass(iSpecies) )
     enddo

     ! Update Ion Densities
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
        endif
     enddo
  enddo !iAlt = 1,nAlts


  ! Update the non-log species densities (NewNS)
  ! Update the non-log total densities (NewNT)
  ! Update the non-log total mass densities (NewRho)
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

  ! Update Vertical Winds
  do iAlt = 1,nAlts
     ! Optional Sponge Layer Specifiction
     if (iAlt >= (nAlts - nAltsSponge)) then
        NuSP = AmpSP*(1.0 - cos( pi*(kSP - (nAlts - iAlt))/kSP) )
     else
        NuSP = 0.0
     endif

     if (UseDamping) then
        VertTau(iAlt) = &
             15 - (1 - exp(-1.0*altitude_G(ialt)/1000.0/40.0))*5.0
     endif

     do iSpecies=1,nSpecies

        ! If we are not subtracting off hydrostatic background
        ! Mainly for minor photochemical species or 
        ! light species H, H2, He which are escaping and are highly
        ! non-hydrostatic over much of the computational domain.
        if (.not. SubtractHydrostatic(iAlt,iSpecies)) then
           NewMomentumS(iAlt,iSpecies) = MomentumS(iAlt,iSpecies) - &
                 DtIn*(AUSMMomentumSFluxes(iAlt,iSpecies)) + &
                 DtIn*RhoS(iAlt,iSpecies)*EffectiveGravity(iAlt) 
        else
        ! More standard case for major species
        ! Background has been subtracted.
        ! Use only DeviationRhoS for acceleration due to gravity
           NewMomentumS(iAlt,iSpecies) = MomentumS(iAlt,iSpecies) - &
                 DtIn*(AUSMMomentumSFluxes(iAlt,iSpecies)) + &
                 DtIn*DeviationRhoS(iAlt,iSpecies)*EffectiveGravity(iAlt) 
        endif 
        ! In conservative form, we must add the chemical source/loss
        ! term as an additional momentum source/loss term
        ! ChemSources_1d -> net loss (#/m^3/s) from calc_sources
        NewMomentumS(iAlt,iSpecies) = NewMomentumS(iAlt,iSpecies) + &
              DtIn*ChemSources_1d(iAlt,iSpecies)*&
              VertVel(iAlt,iSpecies)*Mass(iSpecies)

        ! Updated Vertical Winds from the updated momentum
        NewVertVel(iAlt,iSpecies) = &
           NewMomentumS(iAlt,iSpecies)/NewRhoS(iAlt,iSpecies) 

        ! Add Explicit Spherical Curvature Terms here
        NewVertVel(iAlt,iSpecies) = &
           NewVertVel(iAlt,iSpecies) + DtIn*&
            (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
            * InvRadialDistance_C(iAlt) 
        ! Add Coriolis Acceleration here
        if (UseCoriolis) then
           NewVertVel(iAlt,ispecies) = NewVertVel(iAlt,ispecies) + DtIn * ( &
                Coriolis * Vel_GD(iAlt,iEast_))
        endif

        ! Thermal Diffusion Effects (For Light Species H2, H, and He) 
        ! ThermalDiffCoefS is set in calc_rates
        ! Note:  ThermalDiffCoefS is < 0.0 for light species
        ! Note:  ThermalDiffCoefS is > 0.0 for heavy species
        ! This forces light species into hot zones and heavy species into cold zones
        NewVertVel(iAlt,iSpecies) = NewVertVel(iAlt,iSpecies) - &
          DtIn*(ThermalDiffCoefS(iSpecies)*Boltzmanns_Constant*&
             GradTemp(iAlt))/Mass(iSpecies)
     enddo ! iSpecies
  enddo ! iAlt

  ! Neutral Neutral Friction.
  ! We need this here to ensure that minor species adopt the mean
  ! background gas velocity in the collisional regime.
  ! Use a matrix inversion method to solve stiff term. 
  ! Note:  This is only linear in time (Backward Euler)
  !        but this is embedded in our RK-4, so we expect
   !       overall 4th order accuracy
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


  ! Cap the Winds to reasonable values
  NewVel_GD(-1:nAlts+2,iUp_) = 0.0
  do iAlt = -1, nAlts+2
     do iSpecies=1,nSpecies
        NewVertVel(iAlt, iSpecies) = max(-VelocityCap, &
             NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( VelocityCap, &
             NewVertVel(iAlt, iSpecies))
        ! Calculate bulk vertical winds
        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies)*NewRhoS(iAlt,iSpecies)/NewRho(iAlt)
     enddo
  enddo


  ! Update bulk horizontal momentum transport due to vertical winds 
  ! and energy/temperatures
  do iAlt = 1, nAlts
      ! Horizontal (East/North) momentum updates
      NewMomentum(iAlt,iEast_) = Momentum(iAlt,iEast_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iEast_)) 
 
      NewVel_GD(iAlt,iEast_) = NewMomentum(iAlt,iEast_)/NewRho(iAlt)

      NewMomentum(iAlt,iNorth_) = Momentum(iAlt,iNorth_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iNorth_)) 
 
      NewVel_GD(iAlt,iNorth_) = NewMomentum(iAlt,iNorth_)/NewRho(iAlt)

     ! Energy Update
     NewTotalEnergy(iAlt)   = TotalEnergy(iAlt) - &
         DtIn*AUSMTotalEnergyFluxes(iAlt) + &
         DtIn*Rho(iAlt)*Vel_GD(iAlt,iUp_)*EffectiveGravity(iAlt) 
 
     ! Extract the new pressure from the new Energy
     NewPress(iAlt) = &
        (NewTotalEnergy(iAlt) - &
         0.5*NewRho(iAlt)*(NewVel_GD(iAlt,iUp_)**2.0 + &
                           NewVel_GD(iAlt,iNorth_)**2.0 + &
                           NewVel_GD(iAlt,iEast_ )**2.0 ) )*&
                           (Gamma_1d(iAlt) - 1.0)
 
     ! Extract the New Temp from the new Pressure
     NewTemp(iAlt) = NewPress(iAlt)/(Boltzmanns_Constant*NewNT(iAlt))
  enddo ! iAlt

  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

end subroutine advance_vertical_1stage_ausm

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov_alts_ausm(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

  call calc_facevalues_alts_ausm(Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt_C(1:nAlts)

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt_C(1:nAlts)

end subroutine calc_rusanov_alts_ausm

!\
! ------------------------------------------------------------
! calc_facevalues_alts_ausm
! ------------------------------------------------------------
!/

subroutine calc_facevalues_alts_ausm(Var, VarLeft, VarRight)

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

end subroutine calc_facevalues_alts_ausm


! calc_all_fluxes_hydro:
! Purpose is to calculate the conservative fluxes of the prognostic variables
! in the radial/vertical direction
subroutine calc_all_fluxes_hydro(DtIn, RhoS, PressureS, HydroPressureS, HydroRhoS, &
            HydroPressure, HydroRho, RhoSFlux, MomentumSFlux, &
            EnergyFlux, MomentumFlux, &
            RadDist, SubtractHydrostatic)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax, VertVel, Gamma_1D, &
                          LogRho, Vel_GD, MeanMajorMass_1d, &
                          Temp, Altitude_G
  use ModPlanet, only : nSpecies, nIonsAdvect, Mass, RBody, &
                        iN2_, cSpecies
  use ModGITM, only : Dt,iUp_, iEast_, iNorth_
  use ModConstants, only : Boltzmanns_Constant

  implicit none

! Passed from Vertical_Solver In
  real, intent(in) :: DtIn
  real, intent(in) :: RhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) :: PressureS(-1:nAlts+2,1:nSpecies)
  real, intent(in) :: HydroRhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) :: HydroPressureS(-1:nAlts+2,1:nSpecies)
  real, intent(in) :: HydroRho(-1:nAlts+2) 
  real, intent(in) :: HydroPressure(-1:nAlts+2)
  real, intent(in) :: RadDist(-1:nAlts+2)
  logical, intent(in) :: SubtractHydrostatic(-1:nAlts+2,1:nSpecies)
  ! Conservative Fluxes for output
  real, intent(out):: RhoSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: MomentumSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: EnergyFlux(1:nAlts)
  real, intent(out):: MomentumFlux(1:nAlts,1:3)

  ! Left and Right Face values for our prognostic variables and some 
  ! derived variables.
  ! M12 -> denotes the values at the -1/2 face
  ! P12 -> denotes the values at the +1/2 face
  ! Left and Right refer to which side you are on:
  !
  ! Example is RhoSLeft_M12 is the RhoS at the -1/2 on the Left side

  ! Note: We use LogRho and LogP for our interpolation to the facevalues,
  ! since they are exponential with altitude.
  real, dimension(1:nAlts,1:nSpecies) :: RhoSLeft_M12, RhoSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: RhoSLeft_P12, RhoSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: PressureSLeft_M12, PressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: PressureSLeft_P12, PressureSRight_P12
  ! Log Variables for RhoS, Pressure_{s}
  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogRhoS
  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogPS

  real, dimension(1:nAlts,1:nSpecies) :: LogRhoSLeft_M12, LogRhoSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogRhoSLeft_P12, LogRhoSRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: LogPressureSLeft_M12, LogPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogPressureSLeft_P12, LogPressureSRight_P12
  ! Vertical Velocity for each species at the faces
  real, dimension(1:nAlts,1:nSpecies) :: VelLeft_M12, VelRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: VelLeft_P12, VelRight_P12
  ! Horizontal Velocity at the faces
  real, dimension(1:nAlts,3) :: VelGDLeft_M12, VelGDRight_M12
  real, dimension(1:nAlts,3) :: VelGDLeft_P12, VelGDRight_P12
  ! Total Pressure at the faces
  real, dimension(1:nAlts) :: PLeft_M12, PRight_M12
  real, dimension(1:nAlts) :: PLeft_P12, PRight_P12
  ! Gamma is needed for the Energy
  real, dimension(1:nAlts) :: GammaLeft_M12, GammaRight_M12
  real, dimension(1:nAlts) :: GammaLeft_P12, GammaRight_P12
  ! Total Energy @ faces
  real, dimension(1:nAlts) :: ELeft_M12, ERight_M12
  real, dimension(1:nAlts) :: ELeft_P12, ERight_P12
  ! Total Rho @ faces
  real, dimension(1:nAlts) :: RhoLeft_M12, RhoRight_M12
  real, dimension(1:nAlts) :: RhoLeft_P12, RhoRight_P12
  ! Numerical Speed of Sound at the faces (Liou et al. 2006, AUSM+-up)
  real, dimension(1:nAlts) :: CSLeft_M12, CSRight_M12
  real, dimension(1:nAlts) :: CSLeft_P12, CSRight_P12
  ! Fluxes at the faces:  These are needed to calculate the divergence of 
  ! the fluxes
  real, dimension(1:nAlts,1:nSpecies) :: RhoSFlux_M12, RhoSFlux_P12
  real, dimension(1:nAlts,1:nSpecies) :: MomentumSFlux_M12, MomentumSFlux_P12
  real, dimension(1:nAlts,3) :: Momentum_M12, Momentum_P12
  real, dimension(1:nAlts) :: EnergyFlux_M12, EnergyFlux_P12
  ! Hydrostatic Variables:
  !   These are used to subtract off from the numerical fluxes
  real, dimension(-1:nAlts+2,1:nSpecies) :: LogHydroPressureS
  real, dimension(1:nAlts,1:nSpecies) :: LogHydroPressureSLeft_M12, &
                                         LogHydroPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: LogHydroPressureSLeft_P12, &
                                         LogHydroPressureSRight_P12
  real, dimension(1:nAlts,1:nSpecies) :: HydroPressureSLeft_M12, &
                                        HydroPressureSRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: HydroPressureSLeft_P12, &
                                         HydroPressureSRight_P12
  ! It is often useful to use the simple mean of the Left and Right
  ! facevalues for the pressure in the energy fluxes.
  real, dimension(1:nAlts,1:nSpecies) :: MeanHydroPressureS_M12
  real, dimension(1:nAlts,1:nSpecies) :: MeanHydroPressureS_P12
  ! SubCs is simply used to calculate the numerical sound speed
  real :: SubCs
  integer :: iSpecies, iAlt, iDim
  !------------------------------------------------------------
  ! Numerical Velocities from Liou et al. [2006] on the faces
  ! Numerical Pressures from Liou et al. [2006] on the faces
  ! These are used to calculate the actual fluxes at the interface
  ! They are a combination of physical and numerical winds
  real, dimension( 1:nAlts,1:nSpecies) :: NumericalVelocity_P12, &
                                          NumericalVelocity_M12
  real, dimension( 1:nAlts,1:nSpecies) :: NumericalPressure_P12, &
                                          NumericalPressure_M12   
  ! BulkNumericalVelocity is the mean between the left and right 
  ! numerical winds
  real, dimension( 1:nAlts) :: BulkNumericalVelocity_P12, &
                               BulkNumericalVelocity_M12    

  ! Mean Values:  Simple average of the left and right face
  ! values.
  real, dimension( 1:nAlts, 1:nSpecies) :: MeanPressureS_P12,&
                                           MeanPressureS_M12    
  real, dimension( 1:nAlts, 1:nSpecies) :: MeanRhoS_P12, &
                                           MeanRhoS_M12    
  ! Numerical Flux constants
  ! LiouKp, LiouKpS are empircal functions that should be
  ! between 0.0 and 1.0.
  !
  ! We have found that having an altitude-varying LiouKpS works best.
  ! mainly because high densities at low altitudes can give 
  ! very high numerical fluxes relative to the higher altitudes where
  ! there are lower densities
  real :: Ku(1:nAlts,1:nSpecies)
  real :: LiouKp, LiouKu
  real :: LiouKpS(1:nAlts,1:nSpecies), LiouKuS(1:nAlts,1:nSpecies)
  real :: MaxKpS(1:nAlts,1:nSpecies), MaxKuS(1:nAlts,1:nSpecies)
  real :: MinKpS(1:nAlts,1:nSpecies), MinKuS(1:nAlts,1:nSpecies)
  real :: KpWidth, KpAltMidPoint
  integer :: AltIndex
  real :: LogKpS(1:nAlts,1:nSpecies), LogKuS(1:nAlts,1:nSpecies)
  real :: LogMinKpS(1:nAlts,1:nSpecies), LogMinKuS(1:nAlts,1:nSpecies)
  real :: LogMaxKpS(1:nAlts,1:nSpecies), LogMaxKuS(1:nAlts,1:nSpecies)

  ! Variables to calculate the cell face areas and volumes
  real, dimension( 1:nAlts) :: LeftRadius, RightRadius    
  real, dimension( 1:nAlts) :: AreaFunction_P12, AreaFunction_M12    
  real, dimension( 1:nAlts) :: LocalCellVolume

  ! Liou et al and the AUSM Method
  ! require the mach number on the cell faces
  ! Liou also recommends replacing the physical speed of sound
  ! with a "numerical sound speed" to resolve low mach-number
  ! flows.
  ! These varibles are how we get to this numerical speed of sound
  real, dimension(1:nAlts) :: LiouCSLeft_M12, LiouCSRight_M12
  real, dimension(1:nAlts) :: LiouCSLeft_P12, LiouCSRight_P12

  real, dimension(1:nAlts) :: LiouEnthalpyLeft_M12, LiouEnthalpyRight_M12
  real, dimension(1:nAlts) :: LiouEnthalpyLeft_P12, LiouEnthalpyRight_P12

  real, dimension(1:nAlts) :: InterfaceCS_M12, InterfaceCS_P12
  ! Mach Numbers used by Liou et al. [2006]
  real, dimension(1:nAlts,1:nSpecies) :: MLeft_M12, MRight_M12
  real, dimension(1:nAlts,1:nSpecies) :: MLeft_P12, MRight_P12

  real, dimension(1:nAlts,1:nSpecies) :: M2Bar_M12, M2Bar_P12
  real, dimension(1:nAlts,1:nSpecies) :: M2Zero_M12, M2Zero_P12

  real, dimension(1:nAlts,1:nSpecies) :: MZero_M12, MZero_P12
  real:: MInf, LiouBeta
  real, dimension(1:nAlts,1:nSpecies):: ModifiedZeta

  real, dimension(1:nAlts,1:nSpecies) :: FA_M12, FA_P12

  ! Liou et al. [2006] requires functions of the cell face
  ! Mach Numbers.  there are 1st order, 2nd order and 4th order
  ! functions, defined below
  ! First Order Mach number polynomials
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Left_M12, MF1N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Right_M12, MF1N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Left_P12, MF1N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF1P_Right_P12, MF1N_Right_P12

  ! Second Order Mach number polynomials
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Left_M12, MF2N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Right_M12, MF2N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Left_P12, MF2N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF2P_Right_P12, MF2N_Right_P12

  ! Fourth Order Mach number polynomials
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Left_M12, MF4N_Left_M12
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Right_M12, MF4N_Right_M12

  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Left_P12, MF4N_Left_P12
  real, dimension(1:nAlts,1:nSpecies) :: MF4P_Right_P12, MF4N_Right_P12

  real, dimension(1:nAlts,1:nSpecies) :: MPress_M12, MPress_P12
  real, dimension(1:nAlts,1:nSpecies) :: InterfaceMach_M12, InterfaceMach_P12

  !MInf and LiouBeta are constants
  ! MInf is supposed to the be the "free stream" mach number.
  ! Minf cannot be zero without killing the code.
  ! A small but finite value is chosen.
  ! Liou Beta = 1.0/8.0 was suggested as optimal in Liou et al. [2006]
  MInf = 1.0e-19
  LiouBeta = 1.0/8.0

  ! Calculate Logarithm of exponential prognostic variables
  LogRhoS(-1:nAlts+2,1:nSpecies) = alog(RhoS(-1:nAlts+2,1:nSpecies))
  LogPS(-1:nAlts+2,1:nSpecies) = alog(PressureS(-1:nAlts+2,1:nSpecies))
  LogHydroPressureS(-1:nAlts+2,1:nSpecies) = &
                            alog(HydroPressureS(-1:nAlts+2,1:nSpecies))

  Ku(1:nAlts,1:nSpecies) = 2.0             !! Ullrich et al. [2011]

  ! Next, an altitude variable Ku and Kp can be defined.
  ! I have found that an exponentially varying Kp and Ku work best
  ! You can use a constant value by specifying 
  ! MinKp and MaxKp to be the same value
  MinKuS(1:nAlts,1:nSpecies) = 0.25
  MaxKuS(1:nAlts,1:nSpecies) = 0.50

!  MinKpS(1:nAlts,1:nSpecies) = 1.0e-19
!  MaxKpS(1:nAlts,1:nSpecies) = 1.0e-09

  MinKpS(1:nAlts,1:nSpecies) = 0.5
  MaxKpS(1:nAlts,1:nSpecies) = 0.5

  LogMaxKuS = alog(MaxKuS)
  LogMinKuS = alog(MinKuS)

  LogMaxKpS = alog(MaxKpS)
  LogMinKpS = alog(MinKpS)

  ! JMB:  07/10/2017
  ! Below I specify a tanh function for the Kp, Ku
  ! the user is free to change this. 
  ! the tanh was arrived at heuristically
  KpWidth = 15.0e+03
  KpAltMidPoint = 160.0e+03
  do iAlt = 1, nAlts
     do iSpecies = 1, nSpecies
         LogKpS(iAlt,iSpecies) = LogMinKpS(iAlt,iSpecies) +  &
               0.5*(LogMaxKpS(iAlt,iSpecies) - LogMinKpS(iAlt,iSpecies))*&
               ( 1.0 + tanh(  (Altitude_G(iAlt) - KpAltMidPoint)/KpWidth ) )

         LogKuS(iAlt,iSpecies) = LogMinKuS(iAlt,iSpecies) +  &
               0.5*(LogMaxKuS(iAlt,iSpecies) - LogMinKuS(iAlt,iSpecies))*&
               ( 1.0 + tanh(  (Altitude_G(iAlt) - KpAltMidPoint)/KpWidth ) )

        LiouKpS(iAlt,iSpecies) = exp(LogKpS(iAlt,iSpecies))
             Ku(iAlt,iSpecies) = exp(LogKuS(iAlt,iSpecies))
     enddo 
  enddo 
  ! NOTE:  LiouKpS and Ku will determine how much numerical
  ! fluxes are included in the final solution

   ! Calculate the left and right states of the Variables
   ! on both Interfaces (P12 = +1/2, and M12 = -1/2)
    do iSpecies = 1, nSpecies
       call calc_kt_facevalues(LogRhoS(-1:nAlts+2,iSpecies), &
                       LogRhoSLeft_M12( 1:nAlts  ,iSpecies), &
                      LogRhoSRight_M12( 1:nAlts  ,iSpecies), &
                       LogRhoSLeft_P12( 1:nAlts  ,iSpecies), &
                      LogRhoSRight_P12( 1:nAlts  ,iSpecies) )
       ! Using our LogRhoS, get the RhoS on Interfaces
        RhoSLeft_M12(:,iSpecies) = exp( LogRhoSLeft_M12(:,iSpecies)) 
       RhoSRight_M12(:,iSpecies) = exp(LogRhoSRight_M12(:,iSpecies)) 

        RhoSLeft_P12(:,iSpecies) = exp( LogRhoSLeft_P12(:,iSpecies)) 
       RhoSRight_P12(:,iSpecies) = exp(LogRhoSRight_P12(:,iSpecies)) 
    enddo 

    ! Get the Interface values for the Species Hydro Pressures
    do iSpecies = 1, nSpecies
       call calc_kt_facevalues(LogHydroPressureS(-1:nAlts+2,iSpecies), &
                          LogHydroPressureSLeft_M12(1:nAlts,iSpecies), &
                         LogHydroPressureSRight_M12(1:nAlts,iSpecies), &
                          LogHydroPressureSLeft_P12(1:nAlts,iSpecies), &
                         LogHydroPressureSRight_P12(1:nAlts,iSpecies) )

        HydroPressureSLeft_M12(:,iSpecies) = &
              exp( LogHydroPressureSLeft_M12(:,iSpecies)) 
       HydroPressureSRight_M12(:,iSpecies) = &
              exp(LogHydroPressureSRight_M12(:,iSpecies)) 

        HydroPressureSLeft_P12(:,iSpecies) = &
              exp( LogHydroPressureSLeft_P12(:,iSpecies)) 
       HydroPressureSRight_P12(:,iSpecies) = &
              exp(LogHydroPressureSRight_P12(:,iSpecies)) 
    enddo 

    ! Get the Interface values for the Species Pressures
    do iSpecies = 1, nSpecies
       call calc_kt_facevalues(LogPS(:,iSpecies), &
          LogPressureSLeft_M12(:,iSpecies), LogPressureSRight_M12(:,iSpecies), &
          LogPressureSLeft_P12(:,iSpecies), LogPressureSRight_P12(:,iSpecies) )

       PressureSLeft_M12(:,iSpecies) = exp( LogPressureSLeft_M12(:,iSpecies)) 
       PressureSRight_M12(:,iSpecies) = exp(LogPressureSRight_M12(:,iSpecies)) 

       PressureSLeft_P12(:,iSpecies) = exp( LogPressureSLeft_P12(:,iSpecies)) 
       PressureSRight_P12(:,iSpecies) = exp(LogPressureSRight_P12(:,iSpecies)) 

    enddo 

    ! Get the Interface values for the Vertical Winds
    do iSpecies = 1, nSpecies
       call calc_kt_facevalues(VertVel(:,iSpecies), &
                         VelLeft_M12(:,iSpecies), VelRight_M12(:,iSpecies), &
                         VelLeft_P12(:,iSpecies), VelRight_P12(:,iSpecies) )
    enddo 

    ! Bulk Rho Values
    ! Note: we assume that Rho = sum(rho_{s}) on the interfaces
    ! No need to re-calculate face values separately
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
 
    ! Total Pressure at Interfaces
    ! Derived from the species pressures
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

    ! Bulk Horizontal Winds
    do iDim = 1, 3
        call calc_kt_facevalues(Vel_GD(:,iDim), &
                         VelGDLeft_M12(:,iDim), &
                        VelGDRight_M12(:,iDim), &
                      VelGDLeft_P12(:,iDim), &
                        VelGDRight_P12(:,iDim) )
    enddo 

    ! Derive from the Vertical Winds at Interfaces
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
    ! Gamma Interface values
    call calc_kt_facevalues(Gamma_1d(:), GammaLeft_M12(:), GammaRight_M12(:), &
                                         GammaLeft_P12(:), GammaRight_P12(:) )
    ! Total Energy at interfaces
    do iAlt = 1, nAlts
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

    ! Liou et al. [2006] suggest using the Enthalpy for the 
    ! numerical speed of sound.
    ! Calculate the Enthalpy at the cell faces here.
    do iAlt = 1, nAlts 
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
    enddo 
    ! Use Enthalpy to calculate the Sound Speeds
    ! Based upon Liou et al. [2006], AUSM+-up.
    ! Found that the interface speeds that produced the most stable results 
    ! for all speeds were defined below
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
    enddo 

    MeanPressureS_P12(1:nAlts,1:nSpecies) = &
          0.5*(PressureSLeft_P12(1:nAlts,1:nSpecies) + &
              PressureSRight_P12(1:nAlts,1:nSpecies))

    MeanPressureS_M12(1:nAlts,1:nSpecies) = &
          0.5*(PressureSLeft_M12(1:nAlts,1:nSpecies) + &
              PressureSRight_M12(1:nAlts,1:nSpecies))

    MeanHydroPressureS_P12(1:nAlts,1:nSpecies) = &
          0.5*(HydroPressureSLeft_P12(1:nAlts,1:nSpecies) + &
              HydroPressureSRight_P12(1:nAlts,1:nSpecies))

    MeanHydroPressureS_M12(1:nAlts,1:nSpecies) = &
          0.5*(HydroPressureSLeft_M12(1:nAlts,1:nSpecies) + &
              HydroPressureSRight_M12(1:nAlts,1:nSpecies))

    MeanRhoS_P12(1:nAlts,1:nSpecies) = &
          0.5*(RhoSLeft_P12(1:nAlts,1:nSpecies) + RhoSRight_P12(1:nAlts,1:nSpecies))

    MeanRhoS_M12(1:nAlts,1:nSpecies) = &
        0.5*(RhoSLeft_M12(1:nAlts,1:nSpecies) + RhoSRight_M12(1:nAlts,1:nSpecies))


    ! Calculate the Interface Mach Numbers
    ! Used in the Polynomials to calculate the numerical fluxes
    do iAlt = 1, nAlts 
      do iSpecies = 1, nSpecies
         MLeft_M12(iAlt,iSpecies) = &
            VelLeft_M12(iAlt,iSpecies)/InterfaceCS_M12(iAlt)

         MRight_M12(iAlt,iSpecies) = &
            VelRight_M12(iAlt,iSpecies)/InterfaceCS_M12(iAlt)

         MLeft_P12(iAlt,iSpecies) = &
            VelLeft_P12(iAlt,iSpecies)/InterfaceCS_P12(iAlt)

         MRight_P12(iAlt,iSpecies) = &
            VelRight_P12(iAlt,iSpecies)/InterfaceCS_P12(iAlt)
    
      enddo 
    enddo 

    ! M2Bar is the Mean of the square Mach number 
    do iSpecies = 1, nSpecies
      M2Bar_M12(1:nAlts,iSpecies) = &
       0.5*(MLeft_M12(1:nAlts,iSpecies)**2.0 + &
            MRight_M12(1:nAlts,iSpecies)**2.0 )

       M2Bar_P12(1:nAlts,iSpecies) = &
          0.5*(MLeft_P12(1:nAlts,iSpecies)**2.0 + &
              MRight_P12(1:nAlts,iSpecies)**2.0 )
    enddo 
    ! M2Zero:  is a scaled mach number that is limited by MInf and 1.0
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
    ! FA is a scaling function 
    do iSpecies = 1, nSpecies
     do iAlt = 1, nAlts 

       FA_M12(iAlt,iSpecies) = &
              MZero_M12(iAlt,iSpecies)*(2.0 - MZero_M12(iAlt,iSpecies))
       FA_P12(iAlt,iSpecies) = &
              MZero_P12(iAlt,iSpecies)*(2.0 - MZero_P12(iAlt,iSpecies))

     enddo 
    enddo 
    ! First Order Mach number polynomials
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
    ! 2nd Order Mach Number Polynomials
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

    ! 4th Order Mach Number Polynomials
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

  
    !JMB 07/10/2017
    !    FINALLY, Calculate the Full Numerical Fluxes for the AUSM+-up Update
    !    Two Versions:  (1) when hydrostatic background is subtracted
    !                   (2) when the hydrostatic background is NOT subtracted
    ! Note:  We employ an entropy fix suggested by Moguen et al. [2009].
    ! 		This shows up as a ModifiedZeta term
    !		proportional to L/tau -> where L is a characterisitc length
    !		tau is a characteristic time.
    !  
    !	     I use dz/Dt as the characteristic cell length and time
    do iSpecies = 1, nSpecies
      do iAlt = 1, nAlts 

         ModifiedZeta(iAlt,iSpecies) = 1.0
         ! Subtract hydrostatic background
         MPress_M12(iAlt,iSpecies) =      &
           LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_M12(iAlt,iSpecies)), 0.0)*&
           ((PressureSRight_M12(iAlt, iSpecies) - &
           HydroPressureSRight_M12(iAlt,iSpecies) ) - &
           (PressureSLeft_M12(iAlt, iSpecies) - &
            HydroPressureSLeft_M12(iAlt,iSpecies) ) )/ &
           ( MeanRhoS_M12(iAlt,iSpecies)*InterfaceCS_M12(iAlt)*&
           (FA_M12(iAlt,iSpecies)*InterfaceCS_M12(iAlt) + &
            ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 

         MPress_P12(iAlt,iSpecies) = &
            LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_P12(iAlt,iSpecies)), 0.0)*&
            (  (PressureSRight_P12(iAlt, iSpecies) - &
            HydroPressureSRight_P12(iAlt,iSpecies) ) - &
            (PressureSLeft_P12(iAlt, iSpecies) - &
            HydroPressureSLeft_P12(iAlt,iSpecies) ) )/ &
            ( MeanRhoS_P12(iAlt,iSpecies)*InterfaceCS_P12(iAlt)*&
            (FA_P12(iAlt,iSpecies)*InterfaceCS_P12(iAlt) + &
             ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 

       ! When the background is not subtracted
       if (.not. SubtractHydrostatic(iAlt,iSpecies) ) then
         MPress_P12(iAlt,iSpecies) = &
            LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_P12(iAlt,iSpecies)), 0.0)*&
            (PressureSRight_P12(iAlt, iSpecies) - PressureSLeft_P12(iAlt,iSpecies) )/&
            ( MeanRhoS_P12(iAlt,iSpecies)*InterfaceCS_P12(iAlt)*&
            (FA_P12(iAlt,iSpecies)*InterfaceCS_P12(iAlt) + &
            ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 

         MPress_M12(iAlt,iSpecies) = &
           LiouKpS(iAlt,iSpecies)*max( (1.0 - M2Bar_M12(iAlt,iSpecies)), 0.0)*&
           (PressureSRight_M12(iAlt, iSpecies) - PressureSLeft_M12(iAlt,iSpecies) )/&
           ( MeanRhoS_M12(iAlt,iSpecies)*InterfaceCS_M12(iAlt)*&
           (FA_M12(iAlt,iSpecies)*InterfaceCS_M12(iAlt) + &
            ModifiedZeta(iAlt,iSpecies)*dAlt_C(iAlt)/DtIn)  ) 
       endif 
      enddo ! iSpecies
    enddo ! iAlt

    ! Use the MPress to update the Interface Mach Numbers and the Numerical 
    ! Velocities
    do iAlt = 1, nAlts 
       do iSpecies = 1, nSpecies
          InterfaceMach_M12(iAlt,iSpecies) =  &
                 MF4P_Left_M12(iAlt,iSpecies) + MF4N_Right_M12(iAlt,iSpecies) &
                    - MPress_M12(iAlt,iSpecies)

          NumericalVelocity_M12(iAlt,iSpecies) = &
                   InterfaceCS_M12(iAlt)*InterfaceMach_M12(iAlt,iSpecies) 

          InterfaceMach_P12(iAlt,iSpecies) =  &
                 MF4P_Left_P12(iAlt,iSpecies) + MF4N_Right_P12(iAlt,iSpecies) &
                    - MPress_P12(iAlt,iSpecies)
 
          NumericalVelocity_P12(iAlt,iSpecies) = &
                   InterfaceCS_P12(iAlt)*InterfaceMach_P12(iAlt,iSpecies) 
       enddo ! iSpecies 
    enddo  ! iAlt

    ! Numerical Pressure calculation
    ! We follow the Ullrich [2011] method of using the mean Pressure
    ! at the interface and we use the modified
    ! version of the numerical Pressure Fluxes
    do iSpecies = 1, nSpecies
       do iAlt = 1, nAlts
          NumericalPressure_P12(iAlt,iSpecies) = &
                  (MeanPressureS_P12(iAlt,iSpecies) - &
              MeanHydroPressureS_P12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*InterfaceCS_P12(iAlt)*&
               ( (RhoSRight_P12(iAlt,iSpecies) )*VelRight_P12(iAlt,iSpecies) - &
                  (RhoSLeft_P12(iAlt,iSpecies) )* VelLeft_P12(iAlt,iSpecies) )

          NumericalPressure_M12(iAlt,iSpecies) = &
               (MeanPressureS_M12(iAlt,iSpecies) - &
           MeanHydroPressureS_M12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*InterfaceCS_M12(iAlt)*&
               ( (RhoSRight_M12(iAlt,iSpecies) )*VelRight_M12(iAlt,iSpecies) - &
                  (RhoSLeft_M12(iAlt,iSpecies) )* VelLeft_M12(iAlt,iSpecies) )
          if (.not. SubtractHydrostatic(iAlt,iSpecies) ) then

             NumericalPressure_P12(iAlt,iSpecies) = &
               (MeanPressureS_P12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*InterfaceCS_P12(iAlt)*&
               ( (RhoSRight_P12(iAlt,iSpecies) )*VelRight_P12(iAlt,iSpecies) - &
                  (RhoSLeft_P12(iAlt,iSpecies) )* VelLeft_P12(iAlt,iSpecies) )

             NumericalPressure_M12(iAlt,iSpecies) = &
               (MeanPressureS_M12(iAlt,iSpecies) )&
                - 0.5*Ku(iAlt,iSpecies)*InterfaceCS_M12(iAlt)*&
               ( (RhoSRight_M12(iAlt,iSpecies) )*VelRight_M12(iAlt,iSpecies) - &
                  (RhoSLeft_M12(iAlt,iSpecies) )* VelLeft_M12(iAlt,iSpecies) )

          endif 
           
       enddo !!iAlt = 1, nAlts
    enddo !!! nSpecies

    ! JMB:  07/10/2017
    ! CALCULATe NUMERICAL FLUXES HERE
    ! Use the Liou et al. [2006] to determine the proper flux to use
    ! based upon the wind velocity at the interface.
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

    ! Use the Bulk Velocity in the Energy Fluxes
    ! and the horizontal momentum fluxes
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

    ! Calculate cell face areas and cell volumes
    do iAlt = 1, nAlts
          LeftRadius(iAlt) = 0.5*(RadDist(iAlt) + RadDist(iAlt-1))
          RightRadius(iAlt) = 0.5*(RadDist(iAlt) + RadDist(iAlt+1))
    enddo 
    do iAlt = 1, nAlts
          AreaFunction_P12(iAlt) = RightRadius(iAlt)**2.0
          AreaFunction_M12(iAlt) = LeftRadius(iAlt)**2.0
          LocalCellVolume(iAlt) = &
             (1.0/3.0)*(RightRadius(iAlt)**3.0 - LeftRadius(iAlt)**3.0)
    enddo 

    ! CALCULATE THE FLUX DIVERGENCE FOR OUTPUTS
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


! Calculate the KT-facevalues based upon 2nd order methods
! See Kurganov-Tadmoor, 2000, Journal of Computational Physics.
!  Vol 160, 241-282
! KT = Kurganov-Tadmoor Central method
! and using non-unform grid methods
subroutine calc_kt_facevalues(Var, VarLeft_M12, VarRight_M12, &
                                   VarLeft_P12, VarRight_P12)
 
    use ModVertical, only: dAlt_F, InvDAlt_F
    use ModSizeGITM, only: nAlts
    use ModLimiterGitm
 
    implicit none
    
    real, intent(in) :: Var(-1:nAlts+2)
    real, intent(out):: VarLeft_P12(1:nAlts), VarRight_P12(1:nAlts)
    real, intent(out):: VarLeft_M12(1:nAlts), VarRight_M12(1:nAlts)
  
    real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)
  
    integer :: i

    ! We use the fact that we have non-uniform grids with the
    ! InvDAlt_F variables
    do i=0,nAlts+1
         dVarUp            = (Var(i+1) - Var(i  ))   * InvDAlt_F(i+1)
         dVarDown          = (Var(i  ) - Var(i-1)) * InvDAlt_F(i)
         dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)
    end do

    ! We use the KT-central method to get the facevalues at either
    ! interface.
    do i=1,nAlts
        VarLeft_M12(i) = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i-1)
       VarRight_M12(i) = Var(i  ) - 0.5*dVarLimited(i  ) * dAlt_F(i  )
      
        VarLeft_P12(i) = Var(i  ) + 0.5*dVarLimited(i  ) * dAlt_F(i  )
       VarRight_P12(i) = Var(i+1) - 0.5*dVarLimited(i+1) * dAlt_F(i+1)
    end do
 end subroutine calc_kt_facevalues



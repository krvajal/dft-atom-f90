
program dftatom
! --------------------------------------------------------------
! Program to solve the radial part of Schroedinger equation 
! for the hidrogen atom
! --------------------------------------------------------------
use quad, only: trapz
use types
use ode, only: verlet
use dft
use io_utils
use mesh, only: linspace
use roots

implicit none

    ! call hidrogen_orbitals
    ! call hidrogen
    ! call helio
    call litio


contains
   subroutine litio


    implicit none

    type(SystemParams) :: system

    print *, "========================================="
    print *, "DFT program for Litium"
    print *, "Calculation using LSDA-XC with PZ81 parametrization"
    print *, "========================================="

    system%rmax = 30.0_dp
    system%NumRPoints = 2000
    system%rpoints = linspace(0.0_dp,system%rmax, system%NumRPoints)

    system%Z = 3
    system%NumParticles = 3
    system%NumUp = 2
    system%NumDown = 1
    system%DeltaStep = system%rpoints(2) - system%rpoints(1)
    system%PotType = LSDA_XC_PZ81
    print *, "Delta Step", system%DeltaStep

    call SelfConsistencyLoop(system)
    print *, "==============================",new_line('a')
    system%PotType = LSDA_XC_CY16
    call SelfConsistencyLoop(system)

       
   end subroutine litio

subroutine helio


    implicit none

    type(SystemParams) :: system

    print *, "========================================="
    print *, "DFT program for Helio"

    print *, "========================================="

    system%rmax = 30.0_dp
    system%NumRPoints = 3000
    system%rpoints = linspace(0.0_dp,system%rmax, system%NumRPoints)

    system%Z = 2
    system%NumParticles = 2
    system%NumUp = 1
    system%NumDown = 1
    system%DeltaStep = system%rpoints(2) - system%rpoints(1)
    system%PotType =  LSDA_XC_PZ81
    print *, "Delta Step", system%DeltaStep

    call SelfConsistencyLoop(system)


       
end subroutine helio


subroutine hidrogen


    implicit none

    type(SystemParams) :: system

    print *, "========================================="
    print *, "DFT program for Hidrogen"

    print *, "========================================="

    system%rmax = 30.0_dp
    system%NumRPoints = 3000
    system%rpoints = linspace(0.0_dp,system%rmax, system%NumRPoints)

    system%Z = 1
    system%NumParticles = 1
    system%NumUp = 1
    system%NumDown = 0
    system%DeltaStep = system%rpoints(2) - system%rpoints(1)
    system%PotType =  LSDA_XC_CY16
    system%ElementName = "Hidrogen"

    print *, "Delta Step", system%DeltaStep

    call SelfConsistencyLoop(system)
       
end subroutine hidrogen


subroutine hidrogen_orbitals()
    implicit none
    ! computes the hidrogen atom orbitals for a given energy
    real(dp),dimension(3) :: energies = [ -0.5_dp, -0.125_dp, -1.0_dp/18.0_dp]
    integer,parameter :: points = 50000
    integer ::  i
    real(dp) :: RadialFunc(points), RadialPoints(points),Orbital(points)


    RadialPoints = linspace(0.0_dp,50.0_dp,points)

    do  i = 1,3
        print *, "Solving for E = ", energies(i)
        RadialPoints(1) = 1
        RadialFunc = -2 * (energies(i) + 1/RadialPoints)
        RadialPoints(1) = 0.0_dp
        call solve_radial(RadialFunc,Orbital,RadialPoints)
        call savetxt("hidrogen_orbital_E_"//str(energies(i),4)//".txt",RadialPoints, Orbital)
        


    enddo
   end subroutine hidrogen_orbitals



end program  dftatom



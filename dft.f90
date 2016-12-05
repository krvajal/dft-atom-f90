module dft
    
    use ode, only: numerov
    use quad, only: trapz
    use roots
    use utils
    use constants
    use lsda
    use types
    use solvers
    implicit none
        
    integer,parameter   :: LSDA_XC_PZ81 = 11
    integer,parameter   :: LSDA_XC_CY16 = 12
    integer,parameter   :: LDA_XC_CY16  = 112 ! without spin polarization
    integer,parameter   :: LDA_XC_PZ81  = 113 ! without spin polarization



    ! =========================================
    ! basic type with all the info 
    ! needed to run the dft program
    ! Running Instructions:
    !     To perform a dft calculation we need to
    !     instantiate this structure and fill all 
    !     fields properly, then call the subroutine
    !     SelfConsistencyLoop 
    ! =========================================

    type SystemParams
        character(len = 100) :: ElementName 
        real(dp) :: rmax
        integer :: NumRPoints
        real(dp),allocatable, dimension(:) :: RPoints
        integer :: Z
        integer :: NumParticles
        integer :: NumUp, NumDown
        real(dp) :: DeltaStep
        integer  :: PotType
        integer  :: MaxIter
    end type SystemParams
        


contains


    subroutine SelfConsistencyLoop(system)
        implicit none
        type(SystemParams) :: system
        integer :: iter,i,l
        logical :: result
        type(KhonShamOrbitals) :: ks_orbitals
        real(dp) :: TotalEnergy,OldTotalEnergy
        real(dp) :: HartreePotencial(system%NumRPoints), HartreeEnergy
        real(dp) :: ExternalPotential(system%NumRPoints)
        real(dp) :: ExchangePotentialUp(system%NumRPoints), ExchangePotentialDown(system%NumRPoints)
        real(dp) :: ExchangeEnergy, CorrelationEnergy
        real(dp) :: CorrelationPotentialUp(system%NumRPoints), CorrelationPotentialDown(system%NumRPoints)
        real(dp) :: PotentialUp(system%NumRPoints), PotentialDown(system%NumRPoints)
        
        real(dp) :: aux(system%NumRPoints)
        real(dp),dimension(system%NumRPoints,2) :: dataout
        ! init orbitals

        ks_orbitals%NumUpOrbitals = system%NumUp
        ks_orbitals%NumDownOrbitals = system%NumDown
       
        allocate(ks_orbitals%PhiUp(system%NumUp,system%NumRPoints))
        call assert(allocated(ks_orbitals%PhiUp))
        allocate(ks_orbitals%PhiDown(ks_orbitals%NumDownOrbitals,system%NumRPoints))
        call assert(allocated(ks_orbitals%PhiDown))
        allocate(ks_orbitals%EigenEnergiesUp(ks_orbitals%NumUpOrbitals))
        call assert(allocated(ks_orbitals%EigenEnergiesUp))
        allocate(ks_orbitals%EigenEnergiesDown(ks_orbitals%NumDownOrbitals))
        call assert(allocated(ks_orbitals%EigenEnergiesDown))
        allocate(ks_orbitals%DensityUp(system%NumRPoints))
        call assert(allocated(ks_orbitals%DensityUp))
        allocate(ks_orbitals%DensityDown(system%NumRPoints))
        call assert(allocated(ks_orbitals%DensityDown))

        ExternalPotential(2:)  = -system%Z/system%RPoints(2:)
        ExternalPotential(1) = 0
        HartreePotencial= 0
        CorrelationPotentialUp = 0
        CorrelationPotentialDown = 0
        ExchangePotentialUp = 0
        ExchangePotentialDown = 0
        OldTotalEnergy = 0

        do iter  = 1, system%MaxIter


            print *,"iter", iter
            PotentialUp = ExternalPotential + HartreePotencial + CorrelationPotentialUp + ExchangePotentialUp
            PotentialDown = ExternalPotential + HartreePotencial + CorrelationPotentialDown + ExchangePotentialDown
            
            
            

        

            call SolveKhonShamOrbitals(system,ks_orbitals, PotentialUp,PotentialDown)
            ! now compute the density
        
           

           
            call ComputeDensity(ks_orbitals%NumUpOrbitals, system%NumRPoints,&
                                system%rpoints, ks_orbitals%PhiUp, ks_orbitals%DensityUp)
            

            call ComputeDensity(ks_orbitals%NumDownOrbitals, system%NumRPoints,&
                                system%rpoints, ks_orbitals%PhiDown, ks_orbitals%DensityDown)

            print *,"computed density"

            ! solve hartree
            call ComputeHartree(system%rpoints,ks_orbitals%DensityUp + ks_orbitals%DensityDown , &
                                HartreePotencial, HartreeEnergy,system%DeltaStep)


            print *, "Computed Hartree Potential"            
            call ComputeExchangeLSDA(system%DeltaStep, system%rpoints,ks_orbitals,&
                                      ExchangePotentialUp, ExchangePotentialDown,ExchangeEnergy)

            if(system%PotType == LSDA_XC_PZ81) then
                print *, "Computing correlation using XC_PZ81"
                call PZ81Correlation(system%DeltaStep, system%rpoints,ks_orbitals,&
                                     CorrelationPotentialUp,  CorrelationPotentialDown, CorrelationEnergy)            
            else if(system%PotType == LSDA_XC_CY16) then
                print *, "Computing correlation using XC_CY16"
                call ChachiyoCorrelation(system%DeltaStep, system%rpoints,ks_orbitals,&
                                     CorrelationPotentialUp,  CorrelationPotentialDown, CorrelationEnergy) 
            else 
                call stop_error("unknown method")
            endif


                    
            ! do i=1,system%NumRPoints
            !     print *, system%rpoints(i), ks_orbitals%PhiUp(1,i),ks_orbitals%PhiDown(1,i),HartreePotencial(i)
            ! enddo
            ! stop


            !compute exchange potential
            print *,ks_orbitals%EigenEnergiesUp
            print *,ks_orbitals%EigenEnergiesDown
            TotalEnergy = 0
            TotalEnergy = TotalEnergy + sum(ks_orbitals%EigenEnergiesUp)
            TotalEnergy = TotalEnergy + sum(ks_orbitals%EigenEnergiesDown)

            TotalEnergy = TotalEnergy - HartreeEnergy
            TotalEnergy = TotalEnergy + ExchangeEnergy
            TotalEnergy = TotalEnergy + CorrelationEnergy   

            aux = system%rpoints**2 * ExchangePotentialUp * ks_orbitals%DensityUp

            TotalEnergy = TotalEnergy - &
                          trapz(aux,system%DeltaStep)*4 *pi

            aux = system%rpoints**2 * ExchangePotentialDown * ks_orbitals%DensityDown

            TotalEnergy = TotalEnergy - &
                          trapz(aux, system%DeltaStep)*4 *pi

             aux = system%rpoints**2 * CorrelationPotentialUp * ks_orbitals%DensityUp
             TotalEnergy = TotalEnergy - &
                          trapz(aux, system%DeltaStep)*4 *pi
            aux = system%rpoints**2 * CorrelationPotentialDown * ks_orbitals%DensityDown
            TotalEnergy = TotalEnergy - &
                          trapz(aux, system%DeltaStep)*4 *pi


            print *, "Total Energy",TotalEnergy


            if (abs(TotalEnergy - OldTotalEnergy) < 1e-5) then

                print *, new_line('a')
                print *,"Convergence achieved"
                print *,"Iterations count",iter
                print *, "EigenEnergies Up", ks_orbitals%EigenEnergiesUp
                print *, "EigenEnergies Down", ks_orbitals%EigenEnergiesDown
                print *,"HartreeEnergy", HartreeEnergy
                print *,"XCEnergy", ExchangeEnergy + CorrelationEnergy
                print *,"==============================="
                print *, "Total Energy", TotalEnergy;

                dataout(:,1) = system%rpoints
                do l = 1, ks_orbitals%NumUpOrbitals
                    dataout(:,2) = ks_orbitals%PhiUp(l,:)    
                    call savetxt(system%ElementName//"_orb_up_"//str(l)//".txt",dataout)    
                enddo
                
                do l = 1, ks_orbitals%NumDownOrbitals
                    dataout(:,2) = ks_orbitals%PhiDown(l,:)    
                    call savetxt(system%ElementName//"_orb_down_"//str(l)//".txt",dataout)    
                enddo

                call savetxt(system%ElementName//"_hartree.txt",dataout)
                dataout(:,2) = ExchangePotentialUp
                call savetxt(system%ElementName//"_exchange_up.txt",dataout)
                dataout(:,2) = ExchangePotentialDown
                call savetxt(system%ElementName//"_exchange_down.txt",dataout)
                dataout(:,2) = CorrelationPotentialUp
                call savetxt(system%ElementName//"_correlation_up.txt",dataout)
                dataout(:,2) = CorrelationPotentialDown
                call savetxt(system%ElementName//"_correlation_down.txt",dataout)
                dataout(:,2) = CorrelationPotentialUp + ExchangePotentialUp
                call savetxt(system%ElementName//"_xc_up.txt",dataout)
                dataout(:,2) = CorrelationPotentialdown + ExchangePotentialdown
                call savetxt(system%ElementName//"_xc_down.txt",dataout)
                exit
            endif
            OldTotalEnergy = TotalEnergy
        enddo


    end subroutine SelfConsistencyLoop







    subroutine ComputeDensity(NumOrbitals, NumPoints, rpoints,orbitals, Density)
        integer,intent(in) :: NumOrbitals,NumPoints
        real(dp),intent(in)  :: rpoints(NumPoints)
        real(dp), intent(in) :: orbitals(NumOrbitals,NumPoints)
        real(dp) :: Density(NumPoints)
        real(dp) :: aux(NumPoints)
        integer  ::  i
        Density = 0
        
        

        do i  = 1, NumOrbitals
            aux  = orbitals(i,:)
            aux(2:) = aux(2:)/rpoints(2:)
            Density = Density + aux**2
        enddo
        Density(1) = 0
    end subroutine ComputeDensity


    subroutine SolveKhonShamOrbitals(system,Orbitals, EffectivePotentialUp, EffectivePotentialDown)
            integer :: i
            type(KhonShamOrbitals),intent(inout) :: Orbitals
            real(dp),intent(in) :: EffectivePotentialUp(:), EffectivePotentialDown(:)
            type(SystemParams) :: system
            real(dp),dimension(size(EffectivePotentialUp)) :: phi

            
            ! solve up
            do i = 1, Orbitals%NumUpOrbitals
                phi = Orbitals%PhiUp(i,:)
                call SolveOrbital(system,EffectivePotentialUp,&
                                  EigenEnergyIndex  = i, &
                                  Orbital = phi, &
                                  Energy  =Orbitals%EigenEnergiesUp(i))
                Orbitals%PhiUp(i,:) = phi
            enddo
            
            !solve down
            do i = 1, Orbitals%NumDownOrbitals
                phi = Orbitals%PhiDown(i,:)

                call SolveOrbital(system,EffectivePotentialDown, &
                                 EigenEnergyIndex = i, &
                                 Orbital = phi, &
                                 Energy = Orbitals%EigenEnergiesDown(i))
                Orbitals%PhiDown(i,:) = phi
            enddo


    end subroutine SolveKhonShamOrbitals


    subroutine SolveOrbital(system,Potential, EigenEnergyIndex, Orbital, Energy)
        implicit none
        type(SystemParams),intent(in) :: system
        integer, intent(in) :: EigenEnergyIndex
        real(dp),intent(in) :: Potential(system%NumRPoints)
        real(dp), intent(inout) :: Orbital(system%NumRPoints)
        real(dp), intent(inout) :: Energy
        real(dp)  :: Func(system%NumRPoints)
        integer :: intervals(2)
        integer :: interv_found
        real(dp) :: RadialFunc(system%NumRPoints)
        integer :: i

        integer,parameter :: steps = 10000
        real(dp),parameter :: E_start= -10.0, E_end=-0.01,&
                              E_step = (E_end - E_start)/(steps-1)

        real(dp) :: E(steps), u0(steps)
        
        

        E = (/ (E_start + i * E_step, i=0,steps-1)/)
        
        !compute u0 for each  energy in array
        call ComputeU0(E,system%RPoints, E_step,Potential, u0)

        ! find one zero
        call bracket_zeros(u0,EigenEnergyIndex, intervals,interv_found)
        if (interv_found < EigenEnergyIndex) stop "not zero found"
        Energy = bisec(UStartFunctor,(/system%RPoints, Potential /), E(intervals(EigenEnergyIndex)),&
                        E(intervals(EigenEnergyIndex)+1),0.000001_dp) 

        RadialFunc = -2 * (Energy - Potential)
        call SolveSchodingerRadial(RadialFunc,Orbital,system%RPoints)

    end subroutine SolveOrbital



    subroutine ComputeU0(E,r,h,Pot,u0)
        implicit none
        real(dp), intent(in) :: E(:),h,Pot(:),r(:)
        real(dp), intent(out) :: u0(:)
        real(dp),dimension(size(r)) :: u
        integer      :: steps, i
        real(dp),dimension(size(r)) :: radial

        steps = size(E)
        
        do i = 1, steps
            ! compute radial funct
            
            radial = -2 * (E(i) - Pot)
            
            call SolveSchodingerRadial(radial,u,r)
            
            u0(i) = u(1)
        end do
        
    end subroutine ComputeU0


    subroutine ComputeHartree(rpoints,rho, PotencialHartree,EnergyHartree,h)
        implicit none
        real(dp), intent(in)  :: rho(:),rpoints(:)
        real(dp), intent(out) :: PotencialHartree(:), EnergyHartree
        real(dp), intent(in)  :: h    


        
        call SolvePoisson(rpoints,rho,PotencialHartree,h)

        EnergyHartree = 0.5* trapz(rho*PotencialHartree * rpoints**2,h)*4*pi

    end subroutine ComputeHartree

   

end module dft

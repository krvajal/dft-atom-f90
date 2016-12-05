module lsda

use utils
use quad
use constants !pi, 
use types
implicit none

type KhonShamOrbitals
        
        integer :: NumUpOrbitals
        integer :: NumDownOrbitals
        real(dp), allocatable,dimension(:,:) :: PhiUp, PhiDown
        real(dp),allocatable,dimension(:) :: EigenEnergiesUp,EigenEnergiesDown
        
        real(dp),allocatable, dimension(:) :: DensityUp,DensityDown
end type KhonShamOrbitals



contains

    !=============================================
    ! a standard interpolation formula,
    ! first proposed by von Barth and Hedin
    !=============================================
    function interp_function(zeta)

        implicit none
        real(dp) :: zeta(:)
        real(dp),dimension(:),allocatable :: interp_function
        real(dp) :: exponent = 4.0_dp/3.0_dp

        interp_function = (1 + zeta)**exponent + (1 - zeta)**exponent - 2
        interp_function = interp_function/((2.0_dp)**exponent - 2.0_dp )
    end function interp_function


    !====================================
    ! Computes the exchange energy and 
    ! potential in the LSDA approximation
    !====================================

    subroutine ComputeExchangeLSDA(delta, r, orbitals, ExchangePotentialUp,  ExchangePotentialDown, ExchangeEnergy)
        real(dp), intent(in) :: delta
        type(KhonShamOrbitals),intent(in) :: orbitals
        real(dp),intent(in) :: r(:)
        real(dp),intent(out) :: ExchangePotentialUp(:)
        real(dp),intent(out) :: ExchangePotentialDown(:)
        real(dp),intent(out) :: ExchangeEnergy
        real(dp)   :: ExchangeEnergyDensityUnpolarized(size(r)),ExchangeEnergyDensity(size(r))

        real(dp) :: exponent = 1.0_dp/3.0_dp
        real(dp) :: zeta(size(r)),TotalDensity(size(r))

        TotalDensity  = orbitals%DensityUp + orbitals%DensityDown

        zeta = (orbitals%DensityUp - orbitals%DensityDown)/TotalDensity

        zeta(1) = 0.0 ! null the singularity

        ExchangeEnergyDensityUnpolarized = -(3.0_dp/4.0_dp)*(3.0_dp/pi * TotalDensity)**exponent        
        
        ExchangeEnergyDensity =  ExchangeEnergyDensityUnpolarized &
                                 + interp_function(zeta) * ExchangeEnergyDensityUnpolarized*(2.0_dp**exponent - 1.0_dp)
        

        ExchangeEnergy = 4*  pi * trapz(TotalDensity* ExchangeEnergyDensity * r**2, delta) 


        ExchangePotentialUp = - (6.0_dp/pi * orbitals%DensityUp)**exponent

        ExchangePotentialDown = - (6.0_dp/pi * orbitals%DensityDown)**exponent

    end subroutine ComputeExchangeLSDA


    ! ========================================================================
    ! Parametrization of the correlation energy proposed by Chachiyo
    ! the input parameter are self documented
    ! Chachiyo, T. (2016). Communication: Simple and accurate uniform electron gas
    ! correlation energy for the full range of densities. Journal of Chemical Physics,
    ! 145(2), 9–12. http://doi.org/10.1063/1.4958669
    ! =========================================================================

    subroutine ChachiyoCorrelation(delta, r, orbitals, CorrelationPotentialUp,  CorrelationPotentialDown, CorrelationEnergy)
        real(dp), intent(in) :: delta
        type(KhonShamOrbitals),intent(in) :: orbitals
        real(dp),intent(in) :: r(:)
        real(dp),intent(out) :: CorrelationPotentialUp(:)
        real(dp),intent(out) :: CorrelationPotentialDown(:)
        real(dp),intent(out) :: CorrelationEnergy
        real(dp)   :: CorrelationEnergyDensityUnpolarized(size(r))
        real(dp)   :: CorrelationEnergyDensityPolarized(size(r))
        real(dp)   :: CorrelationEnergyDensity(size(r))
        real(dp), dimension(size(r)) :: dfdzeta, drsdn, dec0drs, dec1drs, dzetadnup, dzetadndown
        real(dp),dimension(2)  :: a = [-0.01554534543482745, -0.007772672717413725];
        real(dp),dimension(2)  :: b = [20.4562557, 27.4203609];
        real(dp):: coeff = 1.0_dp/3.0_dp
        real(dp),dimension(size(r)) ::  rs_inv
        real(dp) :: zeta(size(r)),TotalDensity(size(r))
        

        TotalDensity = orbitals%DensityUp + orbitals%DensityDown

        
        TotalDensity(1) = 1
        zeta = (orbitals%DensityUp - orbitals%DensityDown)/TotalDensity
        zeta(1) = 0
        
        
        
        rs_inv = (4*pi*TotalDensity/3.0_dp)**coeff
        rs_inv(1) = 0
        CorrelationEnergyDensityUnpolarized = a(1) * log(1 +  b(1)*rs_inv*(1 + rs_inv )) 
        CorrelationEnergyDensityPolarized = a(2) * log(1 +  b(2)*rs_inv*(1 + rs_inv ))

        
        CorrelationEnergyDensity  = CorrelationEnergyDensityUnpolarized + interp_function(zeta)*&
                                    (CorrelationEnergyDensityPolarized - CorrelationEnergyDensityUnpolarized)

        CorrelationEnergy =  4*pi* trapz(r**2 * TotalDensity*CorrelationEnergyDensity,delta)

       ! $frac{df}{d\zeta}$
        
        dfdzeta = 4*coeff/(2.0_dp**(4.0*coeff) - 2.0_dp)*( (1 + zeta)**coeff - (1 - zeta)**coeff)
        
        dfdzeta(1) = 0

       ! $frac{d r_s}{dn}$
       drsdn =  -coeff * (3.0_dp/(4*pi))**coeff * (TotalDensity **(-4.0*coeff))
       drsdn(1) = 0
        
        

        dec0drs = -a(1)*b(1)*(1 + 2*rs_inv)*rs_inv*rs_inv/(1 + b(1)* rs_inv *(1 + rs_inv))
        dec1drs = -a(2)*b(2)*(1 + 2*rs_inv)*rs_inv*rs_inv/(1 + b(2)* rs_inv *(1 + rs_inv))

        dzetadnup =  2* orbitals%DensityDown /(TotalDensity**2)
        dzetadndown = -2 * orbitals%DensityUp/(TotalDensity**2)



        TotalDensity(1) = 0
        dzetadnup(1) = 0
        dzetadndown(1) = 0

       CorrelationPotentialUp = CorrelationEnergyDensity  + &
                                 TotalDensity*(dec0drs *drsdn + interp_function(zeta)*(dec1drs - dec0drs)*drsdn)
       CorrelationPotentialDown = CorrelationPotentialUp



       CorrelationPotentialUp = CorrelationPotentialUp + &
                                TotalDensity* (CorrelationEnergyDensityPolarized - CorrelationEnergyDensityUnpolarized)* &
                                dfdzeta * dzetadnup 
       CorrelationPotentialDown = CorrelationPotentialDown + &
                                TotalDensity* (CorrelationEnergyDensityPolarized - CorrelationEnergyDensityUnpolarized)* &
                                dfdzeta * dzetadndown


    end subroutine ChachiyoCorrelation

    ! ========================================================================
    ! Parametrization of the correlation energy proposed by
    ! Pardew and Zunger
    ! Perdew, J. P., & Zunger, A. (1981).
    ! Self-interaction correction to density-functional approximations for
    ! many-electron systems.
    ! Physical Review B, 23(10), 5048–5079. http://doi.org/10.1103/PhysRevB.23.5048
    ! =========================================================================
    subroutine PZ81Correlation(delta, r, orbitals, CorrelationPotentialUp,  CorrelationPotentialDown, CorrelationEnergy)
        implicit none
        
        real(dp), intent(in) :: delta

        type(KhonShamOrbitals),intent(in) :: orbitals
        real(dp),intent(in) :: r(:)
        real(dp),intent(out) :: CorrelationPotentialUp(:)
        real(dp),intent(out) :: CorrelationPotentialDown(:)
        real(dp),intent(out) :: CorrelationEnergy
        real(dp)   :: CorrelationEnergyDensityU(size(r)), PotentialU(size(r))
        real(dp)   :: CorrelationEnergyDensityP(size(r)), PotentialP(size(r))
        real(dp)   :: CorrelationEnergyDensity(size(r))
        real(dp)   :: rs(size(r)),zeta(size(r)), TotalDensity(size(r))
        real(dp),dimension(:),allocatable :: rsl1 ,rsg1
        integer :: midpoint,N
        real(dp) :: coeff = 1.0_dp/3.0_dp

        real(dp), dimension(size(r)) :: dfdzeta

        N = size(r)


        TotalDensity = orbitals%DensityUp + orbitals%DensityDown
        TotalDensity(1) = 1
        zeta = orbitals %DensityUp - orbitals%DensityDown
        
        rs = (3.0_dp/(4*pi*TotalDensity))**(1.0_dp/3.0_dp)
        zeta(1) = 1
        rs(1) = 0.01
        
        call ComputeCorrelationPZ81(rs,.false., CorrelationEnergyDensityU,PotentialU )
        call ComputeCorrelationPZ81(rs,.true., CorrelationEnergyDensityP,PotentialP )

        CorrelationEnergyDensity = CorrelationEnergyDensityU  + interp_function(zeta)*&
                                    (CorrelationEnergyDensityP - CorrelationEnergyDensityU)

        CorrelationPotentialUp = CorrelationEnergyDensityU + interp_function(zeta)*&
                                (PotentialP - PotentialU)
        CorrelationPotentialDown = CorrelationPotentialUp

        dfdzeta = 4*coeff/(2.0_dp**(4.0*coeff) - 2.0_dp)*( (1 + zeta)**coeff - (1 - zeta)**coeff)
        ! zeta(1) = 1
       

        CorrelationPotentialUp = CorrelationPotentialUp + (CorrelationEnergyDensityP - CorrelationEnergyDensityU ) * &
                                 (1  - zeta)*  dfdzeta
        CorrelationPotentialDown = CorrelationPotentialDown + (CorrelationEnergyDensityP - CorrelationEnergyDensityU ) * &
                                 (-1  - zeta)* dfdzeta



        CorrelationEnergy = 4 * pi * trapz(r**2 * TotalDensity * CorrelationEnergyDensity, delta)                         
       
    end subroutine PZ81Correlation


    ! aux function to compute correlation energy and potential

    subroutine ComputeCorrelationPZ81(rs, polarized, EnergyDensity, Potential)
        implicit none
        real(dp),intent(in) :: rs(:)
        logical, intent(in) :: polarized
        real(dp),intent(out) :: EnergyDensity(:)
        real(dp),intent(out) :: Potential(:)
        integer :: i  = 1
        real(dp),dimension(2) :: A     = [  0.0311_dp,  0.01555_dp]
        real(dp),dimension(2) :: B     = [ -0.048_dp , -0.0269_dp]
        real(dp),dimension(2) :: C     = [  0.0020_dp,  0.0007_dp]
        real(dp),dimension(2) :: D     = [ -0.0116_dp, -0.0048_dp]
        real(dp),dimension(2) :: gamma = [ -0.1423_dp, -0.0843_dp]
        real(dp),dimension(2) :: beta1 = [  1.0529_dp,  1.3981_dp]
        real(dp),dimension(2) :: beta2 = [  0.3334_dp,  0.2611_dp]
        real(dp),dimension(:),allocatable :: rsl1 ,rsg1
        integer :: midpoint

        if(polarized) i = 2

        midpoint = minloc(abs(rs - 1.0_dp),1)
        
        if(rs(midpoint ) < 1) midpoint = midpoint + 1
        ! for rs(:midpoint-1)
        rsl1 = rs(:midpoint-1)
        rsg1 = rs(midpoint:)

        ! compute energy rs < 1
        EnergyDensity(:midpoint-1) =  (A(i)  + C(i)*rsl1)*log(rsl1)&
                                                  + B(i)  + D(i)*rsl1
        ! compute energy rs > 1
        EnergyDensity(midpoint:) = gamma(i)/(1 + beta1(i)*sqrt(rsg1) + beta2(i)*rsg1) 

        ! compute potential rs < 1

        Potential(:midpoint-1) = A(i)*log(rsl1) + B(i) - A(i)/3.0_dp &
                                              + 2.0_dp/3 *C(i)*rsl1*log(rsl1) + (2*D(i) - C(i))*rsl1/3                                            
    
        Potential(midpoint:) =  EnergyDensity(midpoint:)  * (1 + 7.0_dp/6.0_dp * beta1(i)*sqrt(rsg1) + &
                                             beta2(1)*4/3.0_dp * rsg1)

        Potential(midpoint:) = Potential(midpoint:)/(1 + beta1(i)*sqrt(rsg1) + beta2(i)*rsg1)


    end subroutine ComputeCorrelationPZ81

end module lsda
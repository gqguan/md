c     This module defines the common variables etc
      module CommonDef
c     Define the double precise number
      real*8, parameter :: ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0
      real*8, parameter :: THREE = 3.0d0, FOUR = 4.0d0, FIVE = 5.0d0
      real*8, parameter :: HALF = 0.5d0
c     Define the constants
      real*8, parameter :: PI = 3.14159265d0
c     Define the grid number
      Integer, parameter :: COM_RadialGridNum = 21
      Integer, parameter :: COM_AxialGridNum = 21
c     The maximum iterative number      
      Integer, parameter :: MaxIterNum = 100
c     Define the data file names
      character*11, parameter :: mod_def_file = 'tmp_MOD.txt'
      character*11, parameter :: stm_def_file = 'tmp_STM.txt'
c     Define the data file formats
      character*13, parameter :: mod_def_file_fmt = '(A30, 1E12.5)'
      character*13, parameter :: stm_def_file_fmt = '(A30, 2E12.5)'
      
c     Define the physical properties
      type PhysProp
c       Averged molecular weight
        real*8 :: AMW
c       Density, viscosity, and specific heat
        real*8 :: rho, mu, cp
c       Thermal conductivity  
        real*8 :: KM
c       Membrane distillation coefficient
        real*8 :: CM
C       Membrane porosity
        real*8 :: porosity
C       Membrane tortuosity
        real*8 :: tortuosity
C       Membrane pore radius
        real*8 :: PoreRadius        
      end type PhysProp
      
c     Define the composition
      type Composition
        real*8 :: H2O, NaCl, Na, NaClS, Cl
      end type Composition  

c     Define the performance of MD module
      type ModPerfm
c       Permeative flux, and thermal efficiency
        real*8 :: JM, eta
      end type
      
      ! Define the properties of module
      type ModuleParam
c       numbers of hollow fiber membrane tubes
        real*8 :: NUM
c       housing size
        real*8 :: LEN, ID2, OD2
c       sizes of hollow fiber membrane tubes
        real*8 :: ID1, OD1, THK
c       packing density
        real*8 :: PHI
c       membrane area
        real*8 :: AM 
c       tube- and shell-side cross section areas
        real*8 :: CSA1, CSA2
c       physical properties of membrane
        type(PhysProp) :: Membrane
c       performance of MD module
        type(ModPerfm) :: Performance
      end type ModuleParam

      ! Define the properties of stream
      type StreamProp
c       stream id (eg, SW103, PW103...)
        character*5 :: ID
c       mass flow, temperature, and pressure
        real*8 :: W, T, P
c       mass flowing flux, velocity
        real*8 :: G, u
c       Reynolds number, and Nusselt number
        real*8 :: Re, Nu
c       physical properties
        type(PhysProp) :: PhysProp
c       mass fraction of composition 
        type(Composition) :: MassFrac
c       component's molar flowrate, [kmol/s]        
        type(Composition) :: MolarFlow 
      end type StreamProp

c     Define the cost
      type Cost
        real*8 :: Depreciation, Operation
      end type Cost
      
c     Define the process variables, which are calculated
      type ProcUnit
c       unit id (eg, X101, E101...)
        character*4 :: ID
c       unit catalog (eg, module, heater, cooler, cfpump)
        character*6 :: UC 
c       heat, and energy (eg, electricity or work)
        real*8 :: Q, E
c       tube- and shell-side inlet temperatures, [K]
        real*8 :: T1IN, T2IN 
c       tube- and shell-side outlet temperatures, [K]
        real*8 :: T1OUT, T2OUT
c       variable for cost estimation 
        real*8 :: VAR1
c       product cost 
        type(Cost) :: Cost
      end type ProcUnit
      
c     Define physical properties of thermoelectric element
      type TEHP
        character*4 :: ID ! dientifier in process flow diagram
        real*8 :: alpha ! Seebect coefficient (V/K)
        real*8 :: K ! thermal conductance (W/K)
        real*8 :: R ! electric resistance (ohm)
        integer :: NumNode ! number of p-n pair of nodes
        real*8 :: A(2), L(2)
        real*8 :: Km, Rc, Rm
        real*8 :: ThermConductivity(2) ! thermal conductivity of p- and n-node (W/m-K)
        real*8 :: ElecResistivity(2) ! electrical resistivity of p- and n-node (ohm/m)
      end type TEHP

c     Define common variables
      integer :: COM_OPT(16)
      type(ModuleParam) :: COM_MOD
      type(StreamProp), dimension(2) :: COM_SIN, COM_SOUT
      
      contains

c     Set the stream index as input argument
      subroutine SetStreamIndex(idx, stm_in, stm_out)
      integer, intent(in) :: idx
      real*8, dimension(:,:), intent(in) :: stm_in
      real*8, dimension(:,:), intent(out) :: stm_out
      integer, dimension(2) :: shp_stm
      integer :: i
      shp_stm = shape(stm_in)
c     Search the stream with non-zero Na
      do i = 1, shp_stm(2)
        if (stm_in(3,i) .NE. zero) then
          if (i .NE. idx) then
            stm_out(1:shp_stm(1),idx) = stm_in(1:shp_stm(1),i)
            stm_out(1:shp_stm(1),i)   = stm_in(1:shp_stm(1),idx)
          else
            stm_out = stm_in
          end if
        end if
      end do
      end subroutine

c     Initiate module parameters
      subroutine InitMOD(mod)
        type(ModuleParam), intent(inout) :: mod
!	  open(11, file="tmp_InitMOD.log")
!	  write(11, *) "Invoking InitMOD()"
!	  close(11)
        mod%NUM = zero
        mod%LEN = zero
        mod%ID1 = zero
        mod%OD1 = zero
        mod%ID2 = zero
        mod%OD2 = zero
        mod%THK = zero
        mod%PHI = zero
        mod%AM  = zero
        mod%CSA1 = zero
        mod%CSA2 = zero
        mod%membrane%rho = zero
        mod%membrane%mu = zero
        mod%membrane%cp = zero
        mod%membrane%KM = zero
        mod%membrane%CM = zero
        mod%membrane%porosity = zero
        mod%membrane%tortuosity = zero
        mod%membrane%PoreRadius = zero
        mod%performance%JM = zero
        mod%performance%eta = zero           
      end subroutine
      
c     Initiate stream
      subroutine InitStream(stm)
        type(StreamProp), intent(inout):: stm
!	  open(11, file="tmp_InitStream.log")
!	  write(11, *) "Invoking InitStream()"
!	  close(11)
        stm%ID = 'xxxxx'
        stm%W = zero
        stm%T = zero
        stm%P = zero
        stm%G = zero
        stm%u = zero
        stm%Re = zero
        stm%Nu = zero
        stm%PhysProp%rho = zero
        stm%PhysProp%mu = zero
        stm%PhysProp%cp = zero
        stm%PhysProp%KM = zero
        stm%PhysProp%CM = zero
        stm%PhysProp%porosity = zero
        stm%PhysProp%tortuosity = zero
        stm%PhysProp%PoreRadius = zero
        stm%MassFrac%H2O = zero
        stm%MassFrac%NaCl = zero
        stm%MolarFlow%H2O   = zero
        stm%MolarFlow%NaCl  = zero    
        stm%MolarFlow%Na    = zero
        stm%MolarFlow%NaClS = zero
        stm%MolarFlow%Cl    = zero
      end subroutine

c     Calculate the following values in the variable of module:
c     Thickness, membrane area, tube-side and shell-side cross section area,
c     shell-side inner diameter
      subroutine CalcModule(MOD)
        type(ModuleParam), intent(inout) :: MOD
c       Define the local variables
c       option for inner diameter
        integer*4 :: IDopt
!	  open(11, file="tmp_CalcModule.log")
!	  write(11, *) "Invoking CalcModule()"
!	  close(11)
c       Calculate the thickness    
        MOD%THK = (MOD%OD1 - MOD%ID1) / TWO
c       Calculate the membrane area based on inner diameter
        IDopt = COM_OPT(1)
        MOD%AM = MembrArea(MOD, IDopt)
        if ((MOD%PHI .NE. zero).and.(MOD%ID2 .EQ. zero)) then
c         Calculate the ID2 according to PHI
          MOD%ID2 = MOD%OD1 * DSQRT(MOD%NUM / MOD%PHI)
        else
c         Calculate the packing density according to ID2
          MOD%PHI = MOD%NUM*MOD%OD1**two/MOD%ID2**two
        end if
c       Calculate the cross section area
        MOD%CSA1 = MOD%NUM * PI * MOD%ID1**TWO / FOUR
        MOD%CSA2 = PI / FOUR * (MOD%ID2**TWO - MOD%NUM * MOD%OD1**TWO)
c       Calculate the tortuosity
        MOD%Membrane%tortuosity = MOD%Membrane%porosity/
     +                    (one-(one-MOD%Membrane%porosity)**(one/three))
      end subroutine
      
c     Calculate the membrane area        
      real(8) function MembrArea(MOD, opt) 
c       Input the geometric size of MD module
        type(ModuleParam), intent(in) :: MOD 
c       opt = 0 calculate the area based on ID
c           = 1                             OD
c           = 2                             ID+THK
        integer*4, intent(in) :: opt 
        select case(opt)
        case(0)
          MembrArea = PI * MOD%ID1 * MOD%LEN * MOD%NUM
        case(1)
          MembrArea = PI * MOD%OD1 * MOD%LEN * MOD%NUM
        case(2)
          MembrArea = PI * (MOD%ID1 + MOD%THK) * MOD%LEN * MOD%NUM
        end select
      end function
      
      subroutine CalcStream(mod, stm)
c       Calculate the mass flowing flux and Reynolds number
        type(ModuleParam), intent(in) :: mod
        type(StreamProp), dimension(:), intent(inout) :: stm
c       Define local variables
        real*8 :: VolFlow, Length, csa(2), wc(2)
        integer :: NumStream, is
!        open(11, file="tmp_CalcStream.log")
!	  write(11, *) "Invoking CalcModule()"
!	  close(11)
        NumStream = size(stm)

        csa(1) = mod%CSA1
        csa(2) = mod%CSA2
        wc(1) = mod%NUM * mod%ID1 * PI
        wc(2) = mod%NUM * mod%OD1 * PI + PI * mod%ID2

        do is = 1, NumStream
c         Volume flow, [m3/s]
          VolFlow = stm(is)%W / stm(is)%PhysProp%rho
          stm(is)%u = VolFlow/csa(is)
          stm(is)%G = stm(is)%u * stm(is)%PhysProp%rho
          Length = HydraDiam(csa(is), wc(is))
          stm(is)%Re = Re(stm(is)%G, Length, stm(is)%PhysProp%mu)
        end do

      end subroutine  
      
      real*8 function EquivDiam(CSA)
c       Calculate the equivalent diameter
c       (ie, diameter of cycle with the same area)
        real*8, intent(in) :: CSA
        EquivDiam = DSQRT(CSA * FOUR / PI)
      end function
      
      real*8 function HydraDiam(CSA, WC)
c       Calculate the hydraulic diameter
c       (ie, 4 times of the quotient of cross-section area, CSA,
c        and wetting circumference, WC)
        real*8, intent(in) :: CSA, WC
        HydraDiam = four*CSA/WC
      end function

c     Calculate the Reynolds number based on the given mass flux, length, 
c     and viscosity
      real(8) function Re(G, D, mu)
        real*8, intent(in) :: G, D, mu
        Re = G * D / mu
      end function
      
c     Calculate the Peclet Number with the specified fluid velocity, length
c     and thermal diffusivity
      real(8) function Pe(v, L, alpha)
        real*8, intent(in) :: v, L, alpha
        Pe = two*v*L/alpha
      end function
      
c     Calculate the thermal diffusivity with the specified thermal conductivity,
c     density and specific heat
      real(8) function alpha(kappa, rho, cp)
        real*8, intent(in) :: kappa, rho, cp
        alpha = kappa/rho/cp
      end function      
      
c     Convert the unit
      real(8) function UnitConvert(Val, FromUnit, ToUnit)
        real*8, intent(in) :: Val
        character(*), intent(in) :: FromUnit, ToUnit
        select case(FromUnit)
          case("h")
            if (ToUnit == "s") UnitConvert = Val * 3.6d3
          case("s")
            if (ToUnit == "h") UnitConvert = Val / 3.6d3
          case("h-1")
            if (ToUnit == "s-1") UnitConvert = Val / 3.6d3
          case("s-1")
            if (ToUnit == "h-1") UnitConvert = Val * 3.6d3
          case("K")
            if (ToUnit == "C") UnitConvert = Val - 2.7315d2
          case("C")
            if (ToUnit == "K") UnitConvert = Val + 2.7315d2
          case("kg")
            if (ToUnit == "t") UnitConvert = Val*1.0d-3
          case("t")
            if (ToUnit == "kg") UnitConvert = Val*1.0d3
          case("Pa")
            if (ToUnit == "mmHg") UnitConvert = Val/1.01325d5*7.6d2  
        end select
      end function
      
c     Calculate the LMTD
      real(8) function LMTD(DT1, DT2)
        real*8, intent(in) :: DT1, DT2
        LMTD = (DT1-DT2)/DLOG(DT1/DT2)
      end function

c     The averaged molecular weight [kg/kmol]
      real(8) function AvgMolWeight(n, comp_mw, comp_x)
c       Molecular weight vector of each component
        real*8, dimension(:), intent(in) :: comp_mw 
c       Molar fraction vector
        real*8, dimension(:), intent(in) :: comp_x 
        integer :: i
        AvgMolWeight = 0.
        do i = 1, n
          AvgMolWeight = AvgMolWeight+comp_mw(i)*comp_x(i)
        end do
      end function

c     A simple numerical integration for 1-d dispersed points
      real(8) function Integration(x, y)
        real*8, dimension(:), intent(in) :: x, y
        real*8 :: dx
        integer :: n, i
        Integration = zero
        n = size(x)
        do i = 1, n-1
          dx = x(i+1)-x(i)
          Integration = Integration+dx*(y(i+1)+y(i))/two
        end do
      end function
      
c     Grid the domain
      subroutine Grid1D(x0, x1, A, opt)
        real*8, intent(in) :: x0, x1
        real*8, intent(inout), dimension(:) :: A
        integer, intent(out), optional :: opt
        integer :: xPointNum, xSectionNum
        integer :: i
        real*8, allocatable, dimension(:) :: xSpan
!	  open(11, file="tmp_Grid1D.log")
!	  write(11, *) "Invoking Grid1D()"
!	  close(11)
c       Get the x-dimensional size
        xPointNum = size(A, dim=1)
        xSectionNum = xPointNum-1
        A(1) = x0
        allocate (xSpan(xSectionNum))
        do i = 1, xSectionNum
c         Calculate the x-dimensional spans
          xSpan(i) = (x1-x0)/xSectionNum
c         Generate the grid data
          A(i+1) = A(i)+xSpan(i)
        end do
        if (IsRelDiff(A(xPointNum), x1, 1.d-5)) then
          opt = 0
        else
          opt = 1
        end if
        deallocate(xSpan) 
      end subroutine
      
      logical(2) function IsRelDiff(x1, x2, RelTol)
        real*8, intent(in) :: x1, x2, RelTol
        if (dabs((x1-x2)/x1) .le. RelTol) then
          IsRelDiff = .true.
        else
          IsRelDiff = .false.
        end if
      end function
      
      ! Calculate the average of input vector
      real(8) function DAVG(fx)
        real*8, intent(in) :: fx(:) 
        integer :: size_fx
        size_fx = size(fx, dim=1)
        DAVG = DSUM(fx)/size_fx
      end function
      
      ! Calculate the sum of input vector
      real(8) function DSUM(fx)
        real*8, intent(in) :: fx(:)
        integer :: size_fx, i
        size_fx = size(fx, dim=1)
        tmp = zero
        do i = 1, size_fx
          tmp = tmp+fx(i)
        end do
        DSUM = tmp
      end function
      
      end module
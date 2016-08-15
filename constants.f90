module constants

	real*8, parameter :: mp_g = 1.*1836.*9.11e-28; !proton mass in grams
	real*8, parameter :: me = 9.11e-31; !kg
    real*8, parameter :: me_g = 9.11e-28
	real*8, parameter :: qe_C = 1.6e-19;
	real*8, parameter :: qe = 4.8032e-10; !statcoulomb
	real*8, parameter :: g = 5./3. !gamma
	real*8, parameter :: eps0 = 8.854187e-12
	integer, parameter :: unt = 50 !unit for opening files
	
	integer :: nspec, neqi, nregions, nsubcycle
	real*8 :: v_piston, v_shell, n_left, T_left, n2n1, T2T1, n_shocked, T_shocked, T_ablation, Tmin
	real*8 :: n_right, T_right
	real*8 :: L, Lab_min
	integer :: maxind, nz, dnz, dt_multi, nquiet, nsmooth, vsmooth, tsmooth, nz0, nz00
	integer :: viscous_species
	real*8 :: dr, eps_visc_max, eps_compress, eps_n, phi, tm_quiet
	logical :: heattransfer_switch, Efield_switch, electron_switch, friction_switch
	logical :: smoothing, bremsstrahlung, ion_viscosity, restart, electron_heat_conduction
	logical :: check_frequency, ifprint
	real*8 :: flimit, ni_frac, nmin, rmin, dt_print
	real*8 :: lm, dtm, dtm_, CFL, cs, coeff_CFL, Tmax, smax, maxTime
	character*9 :: geom
	
	real*8, dimension(:), allocatable :: mi, mi_g, qi, Ai, nleft, nmiddle, nright
	real*8, dimension(:), allocatable :: Zi, r_shell, dr_shell
    real*8, dimension(:,:), allocatable :: p, u, rho, T, ni, r0, den0, temp0, vel0
    real*8, dimension(:,:), allocatable :: U1D, F1D, R1D, U1D_p, U1D_c, G1D, Q_DT, eps_visc, C1D
    real*8, dimension(:,:), allocatable :: T0, N0, V0, mu_ion, AA, AA2
    real*8, dimension(:), allocatable   :: LDIAG, DIAG, UDIAG, BB, IPIV
    real*8, dimension(:), allocatable   :: LDIAG2, DIAG2, UDIAG2, BB2, IPIV2
    real*8, dimension(:), allocatable :: r, xi, xi12, nu0, du, ke, vth_e, q_FS, q_diff, Te, Efield
    real*8, dimension(:), allocatable :: ne, Qextra 	
    real*8, dimension(:,:,:), allocatable :: k_DT, xiab, L_ab	
    real*8 :: tm = 0.0, xmax_visc
    
    save
	
	contains
	
	subroutine initialize()
		implicit none
		integer :: i
		real :: dr_smooth

		geom = "spherical" 

		flimit = 0.06!1.!0.06 !flux limit


		heattransfer_switch = .true. !if .false., no thermal equilibration
		Efield_switch = .true.
		electron_switch = .true.   ! if .false., electrons are uncoupled to ions (no E-field, no thermalization with electrons)	
		electron_heat_conduction = .true.
		friction_switch = .true.
		bremsstrahlung = .true.
		check_frequency = .true. !check if put a max limit on frequency for friction
		smax = 5. ! max time step for friction (only used when check_frequency==true)
		ion_viscosity = .false.
		restart = .false.  !if true, will look at .dat files
		viscous_species = 4
		nsubcycle = 1
		smoothing = .true.
		dr_smooth = 10.e-6 !smoothing range

		coeff_CFL = 5.e-4 !it was 5.-5 yesterday night, with visc=1.e-3
		dt_print = 2.e-12

		nspec = 4
		nz = 2000 !number of zones
		neqi = 3*nspec !total # of equations for the ions
		dnz = 40
		Tmax = 500.e3 * qe_C !in Joules
		Lab_min = 2.

		nregions = 4
		
		L = 5.8e-04
		maxTime = 1.5e-9
		
		call allocate_arrays()
		
		Tmin = 0.5 * qe_C
				
!M =  1.0
!  <A> =  6.5
!  Z   =  2.5
!  T1  =  10.0
!  v_piston =  0.
!  n2n1 =  1.
!  T2T1 =  1.
 
  		!C ions (+4)
		Ai(1) = 12.0
		Zi(1) = 4.0
  		!H ions 
		Ai(2) = 1.0
		Zi(2) = 1.0
		!Ti ions
		Ai(3) = 47.87
		Zi(3) = 4.
		!D ions
		Ai(4) = 1.
		Zi(4) = 1.

		T_ablation = 1.0e3 * qe_C !temperature of expanding gas

		v_piston = 0.
		v_shell  = -225.e3
	
		n_left = 3.0e23 * 1.e6!m-3 this is the density of CH
		T_left = 5. * qe_C!J
	
		n2n1 = 1.
		T2T1 = 1.
		n_shocked = n_left * n2n1
		T_shocked = T_left * T2T1
		
		n_right = 2.4e20 * 1.e6!!m-3 this is the density of D ions 
		T_right = T_left !J
		
		eps_n = 1.e-5  !should not be too much lower than 1.e-4	
		nmin = n_right * eps_n

		!-------------------------------
		!-- parameters for collisions --
		!-------------------------------
		do i = 1, nspec
			mi_g(i) = Ai(i) * mp_g; !mass of <DT> in grams
			mi(i) = 1.e-3 * mi_g(i)
			qi(i) = Zi(i) * qe; !charge of species i in statcoulomb
		enddo
		
		
		!define initial density and temperature
		!--------------------------------------

		!C
		r0(1,:)  = (/ double precision :: 0.0, 310.e-6, 311.e-6, 315.e-6, L /)
		den0(1,:)   = (/ double precision :: nmin, nmin, 0.5*n_left , 2.5e-4*0.5*n_left /)
		temp0(1,:)  = (/ double precision :: Tmin, T_left, T_shocked, T_ablation /)
		vel0(1,:)     = (/ double precision :: 0., v_shell, v_shell, v_shell /)

		!H
		r0(2,:)  = (/ double precision :: 0.0, 310.e-6, 311.e-6, 315.e-6, L /)
		den0(2,:)   = (/ double precision :: nmin, nmin, 0.5*n_left , 2.5e-4*0.5*n_left /)
		temp0(2,:)  = (/ double precision :: Tmin, T_left, T_shocked, T_ablation /)
		vel0(2,:)     = (/ double precision :: 0., v_shell, v_shell, v_shell /)
		
		!Ti		
		r0(3,:)  = (/ double precision :: 0.0, 310.e-6, 310.25e-6, 315.e-6, L /)
		den0(3,:)   = (/ double precision :: nmin, 0.018*n_left, 100.*nmin, 100.*nmin /)
		temp0(3,:)  = (/ double precision :: Tmin, T_left, T_shocked, T_ablation /)
		vel0(3,:)     = (/ double precision :: 0., v_shell, v_shell, v_shell /)
		
		!D		
		r0(4,:)  = (/ double precision :: 0.0, 310.e-6, 311.e-6, 315.e-6, L /)
		den0(4,:)   = (/ double precision :: n_right, 100.*nmin, 100.*nmin, 100.*nmin/)
		temp0(4,:)  = (/ double precision :: Tmin, T_left, T_shocked, T_ablation /)
		vel0(4,:)     = (/ double precision :: 0., v_shell, v_shell, v_shell /)
				
		xmax_visc = 180.e-6 !for r > xmax_visc, ion_visc = ion_visc * 1.e-3

		maxind  = 21200000
		dt_multi = 1
		dr = - L / nz !kg/m^2 !note that dr is negative
		

eps_visc_max = 8.e-3!artificial viscosity
eps_compress = 1.0
		if (geom=="slab") then
			rmin = 0.
		else !spherical
			rmin = 2.5e-6
		endif

		phi = 0.5		
		tm_quiet = 0.
		nquiet = floor(  50.e-6 / abs(dr) )

		cs = max( sqrt( g * T_left / minval(mi) ), sqrt( g * T_right / minval(mi) ) )

		CFL = abs(dr) / cs

		dtm = coeff_CFL * CFL !may need to be varied to avoid NaN when tm > tm_coll
		
		lm = dtm / dr

		if (ion_viscosity) then
			dtm_ = 0.5 * dtm
			eps_visc_max = 1./3 * eps_visc_max!artificial viscosity
		else
			dtm_ = dtm
		endif		

		nsmooth = abs(int(dr_smooth/dr)) !only used if smoothing==.true.
		vsmooth = nsmooth  !smoothing for velocity
		tsmooth = nsmooth  !smoothing for temperature
		
		write(*,*) "smoothing over", nsmooth, "grid points"

	end subroutine initialize

end module
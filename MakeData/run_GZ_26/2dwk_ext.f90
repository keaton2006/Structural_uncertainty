!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Extended 2D-WK equation solver 
!
!                 M. Nakata, Jan. 2022
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE header 

  implicit none

    integer, parameter :: DP = selected_real_kind(14)   ! double precision, e.g., 2.25_DP

!------------------------------------------------------ constants
    real(kind=DP),    parameter :: pi  = 3.141592653589793_DP,  &
                                   twopi  = pi * 2._DP,         &
                                   fourpi = pi * 4._DP,         &
                                   eps = 0.0000000001_DP
    complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )
 
!------------------------------------------------------ numerical parameters
    integer, parameter :: nxx  = 150    ! grid number in x 
    integer, parameter :: nkx  = 200    ! grid number in kx

    ! data-output interval, i.e., 1output per iout steps
    integer, parameter :: iout_0d = 200     ! for eneg.dat
    integer, parameter :: iout_2d = 2000    ! for x-time.dat, x-kx.dat

!------------------------------------------------------ physical parameters
    real(kind=DP),    parameter :: lx   = pi * 8._DP       ! domain of x  (periodic B.C.)
    real(kind=DP),    parameter :: lkx  = 20._DP       ! domain of kx (damping bufferfor |kx|>lkx*0.9)

    real(kind=DP),    parameter :: tend = 4000._DP     ! simulation end time
    real(kind=DP),    parameter :: dt   = 5.d-3       ! time step

    real(kind=DP),    parameter :: ky   = 1._DP       ! poloidal wavenumber (given)

    real(kind=DP),    parameter :: mu_x    = 0.01_DP  ! mu for I
    real(kind=DP),    parameter :: mu_n    = 0.01_DP  ! mu for N0prime
    real(kind=DP),    parameter :: mu_para = 0.01_DP  ! mu_parallel
    real(kind=DP),    parameter :: mu_perp = 0.01_DP  ! mu_perpendicular
    real(kind=DP),    parameter :: dk_l    = 0.01_DP  ! Dk^{L}
    real(kind=DP),    parameter :: dk_nl   = 0.0_DP  ! Dk^{NL}

    real(kind=DP),    parameter :: q0    = 2._DP      ! safety factor
    real(kind=DP),    parameter :: eps_n = 0.070711356_DP    ! Ln/R
    real(kind=DP),    parameter :: CC    = 3._DP     ! C=D_para*k_para
    real(kind=DP),    parameter :: S0    = 0.02_DP      ! mean-flow shearing rate ~0.05
    !real(kind=DP),    parameter :: S0    = 0.0_DP      ! mean-flow shearing rate ~0.05

!------------------------------------------------------ end of physical parameters


!------------------------------------------------------ Physical switches 
    integer, parameter :: idg  = 1    ! 1: coupled with N0prime,   0: decoupled  
    integer, parameter :: ivp  = 1    ! 1: coupled with Vpara,     0: decoupled 
    integer, parameter :: ivx2 = 1    ! 1: d2 Vy/dx2 effect is ON, 0: OFF


!------------------------------------------------------ others
    real(kind=DP),    parameter :: dx   = lx/nxx     ! grid size in x
    real(kind=DP),    parameter :: dkx  = lkx/nkx    ! grid size in kx

    real(kind=DP),dimension(-nxx:nxx) :: xx          !  x-grid 
    real(kind=DP),dimension(-nkx:nkx) :: kx          ! kx-grid 

!------------------------------------------------------ unit number for data output 
    integer, parameter :: oxkx = 500    ! x-kx.dat
    integer, parameter :: oxtm = 501    ! x-time.dat
    integer, parameter :: oeng = 502    ! eneg.dat

    ! --- parameters for the initial condition specified in the main routine
    real(kind=DP) :: gamma0, delomg, deltak, k0, q_x, omega0 

END MODULE 


PROGRAM WK2Dext

  use header 

    real(kind=DP), dimension(-nxx:nxx,-nkx:nkx) :: ff  ! I(x,kx)
    real(kind=DP), dimension(-nxx:nxx) :: gg           ! Vy(x)
    real(kind=DP), dimension(-nxx:nxx) :: hh           ! Vpara(x)                
    real(kind=DP), dimension(-nxx:nxx) :: pp           ! N0prime(x)      
    real(kind=DP), dimension(-nxx:nxx) :: uu           ! Ns(x)        
    real(kind=DP), dimension(-nxx:nxx,-nkx:nkx) :: df, qf
    real(kind=DP), dimension(-nxx:nxx) :: dg, qg
    real(kind=DP), dimension(-nxx:nxx) :: dh, qh
    real(kind=DP), dimension(-nxx:nxx) :: dpp, qp
    real(kind=DP), dimension(-nxx:nxx) :: du, qu
    real(kind=DP), dimension(-nxx:nxx) :: ee, phi2

    real(kind=DP) :: time
    real(kind=DP) :: ee_total, vy2_total, phi2_total
    real(kind=DP) :: vp2_total, N02_total, Ns2_total
    real(kind=DP) :: a0

    integer :: ix, mx, istep, loop
    integer :: icount

! --- grid setting
      do ix = -nxx, nxx
        xx(ix)   = dx * real( ix ) 
      end do 

      do mx = -nkx, nkx
        kx(mx)   = dkx * real( mx ) 
      end do 

! --- initial condition 
      gamma0  = 0.5_DP   ! amplitude
      delomg  = 0.4_DP   ! delta omega
      deltak  = 3._DP    ! delta k
      k0      = 0.0_DP   ! k_0
      q_x     = 0.25_DP   ! k for initial perturb.  
!      q_x     = 0.1_DP   ! k for initial perturb.  
      omega0  = 0.1_DP   ! omega for init. perturb. 

      ff(:,:) = (0._DP)  ! --- zero clear 
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          !ff(ix,mx) = gamma0*exp(-(kx(mx)-k0)**2/deltak**2)/delomg   ! ff = I(x,kx) 
	  ff(ix,mx) = ky*(kx(mx)**2+ky**2)/CC/( 1._DP + kx(mx)**2 + ky**2)**3*exp(-(kx(mx)-k0)**2/deltak**2)/delomg   ! ff = I(x,kx)
        end do
      end do

      a0 = 0.01_DP  ! infinitesimal perturbation amplitude
      do ix = -nxx, nxx
        gg(ix) = a0*sin(q_x*xx(ix))  ! standing wave Vy=A*cos(omega0*t)*sin(q_x*x)
        hh(ix) = (q0/eps_n)*(omega0**2/2._DP/eps_n-eps_n)*a0*sin(q_x*xx(ix))  ! Vpara
        pp(ix) = 0._DP  ! N0prime
        uu(ix) = 0._DP  ! Ns
      end do

      time = 0._DP
      loop = 0
      print*, "dt = ", dt
 
! --- file open 
      open(unit=oxkx, file="x-kx.dat")  
      open(unit=oxtm, file="x-time.dat")  
      open(unit=oeng, file="eneg.dat")  

! --- initial output 
      icount = 0 
      !write(unit=oxkx, fmt=*) "## time, index = ", time, icount
 
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          !write(unit=oxkx, fmt="(3f15.8,SP,256E24.14e3)") time, xx(ix), kx(mx), & 
          !                           ff(ix,mx), gg(ix), hh(ix), pp(ix), uu(ix) 
	 write(unit=oxkx, fmt="(3f15.8,SP,256E24.14e3)") time, xx(ix), kx(mx), ff(ix,mx)                                  
        end do                                             
        !write(unit=oxkx, fmt=*) 
      end do                                             
      !write(unit=oxkx, fmt=*) 
      !write(unit=oxkx, fmt=*) 


      write(unit=oxtm, fmt=*) "## time, index = ", time, icount

      ee(:)   = 0._DP
      phi2(:) = 0._DP
      do ix = -nxx, nxx
        do mx = -nkx, nkx
          ee(ix)   = ee(ix)   + ff(ix,mx) * dkx
          phi2(ix) = phi2(ix) + ff(ix,mx)/(1._DP+kx(mx)**2+ky**2)**2 * dkx
        end do                                             
      end do                                             
       
      do ix = -nxx, nxx
        write(unit=oxtm, fmt="(2f15.8,SP,256E24.14e3)") time, xx(ix), ee(ix), phi2(ix), &
                                                      gg(ix), hh(ix), pp(ix), uu(ix)
      end do  
      write(unit=oxtm,fmt=*)         

      ee_total   = 0._DP
      phi2_total = 0._DP
      vy2_total  = 0._DP
      vp2_total  = 0._DP
      N02_total  = 0._DP
      Ns2_total  = 0._DP

      do ix = -nxx, nxx
        ee_total   = ee_total   + ee(ix)    * dx 
        phi2_total = phi2_total + phi2(ix)  * dx 
        vy2_total  = vy2_total  + gg(ix)**2 * dx 
        vp2_total  = vp2_total  + hh(ix)**2 * dx 
        N02_total  = N02_total  + pp(ix)**2 * dx 
        Ns2_total  = Ns2_total  + uu(ix)**2 * dx 
      end do                                             
                                            
      write(unit=oeng, fmt="(115a)") "## time, ee_tot, ee(0), phi2_tot, phi2(0), "// & 
                                     "vy2_tot, vy(0), vp2_tot, vp(0), N0p2_tot, N0p(0), Ns2_tot, Ns(0)" 
      write(unit=oeng, fmt="(f15.8,SP,256E24.14e3)") time, ee_total, ee(0),   & 
                                                         phi2_total, phi2(0), & 
                                                          vy2_total, gg(0),   &
                                                          vp2_total, hh(0),   &
                                                          N02_total, pp(0),   &
                                                          Ns2_total, uu(0)

! --- loop for time integaration 
   do 
      if ( time > tend - eps ) exit

        time   = time + dt
        loop   = loop + 1

! --- zero clear for work array
        df(:,:) = ( 0._DP )
        dg(:)   = ( 0._DP )
        dh(:)   = ( 0._DP )
        qf(:,:) = ( 0._DP )
        qg(:)   = ( 0._DP )
        qh(:)   = ( 0._DP )

! --- zero clear for high-kx damping 
        do ix = -nxx, nxx
          do mx = -nkx, nkx
            if ( abs(kx(mx)) > lkx*0.9_DP ) then 
              ff(ix,mx) = 0._DP
              df(ix,mx) = 0._DP
            end if
          end do                                             
        end do                                             

! --- for RKG4
        do istep = 1, 4
    
          call calc_delta( ff, gg, hh, pp, uu, df, dg, dh, dpp, du )
 
! --- zero clear for high-kx damping 
          do ix = -nxx, nxx
            do mx = -nkx, nkx
              if ( abs(kx(mx)) > lkx*0.9_DP ) then 
                ff(ix,mx) = 0._DP
                df(ix,mx) = 0._DP
              end if
            end do                                             
          end do                                             
   
          call rkg4( ff, gg, hh, pp, uu, & 
                     df, dg, dh, dpp, du, & 
                     qf, qg, qh, qp, qu, istep )
    
        end do 

! --- data output for field quantities
      if (mod(loop+iout_2d,iout_2d) == 0) then
 
        icount = icount + 1 

        ee(:)   = 0._DP
        phi2(:) = 0._DP
        do ix = -nxx, nxx
          do mx = -nkx, nkx
            ee(ix)   = ee(ix) + ff(ix,mx) * dkx
            phi2(ix) = phi2(ix) + ff(ix,mx)/(1._DP+kx(mx)**2+ky**2)**2 * dkx
          end do                                             
        end do                                             
        
        !write(unit=oxkx, fmt=*) "## time = ", time, icount
        do mx = -nkx, nkx
          do ix = -nxx, nxx
            !write(unit=oxkx, fmt="(3f15.8,SP,256E24.14e3)") time, xx(ix), kx(mx), & 
            !                           ff(ix,mx), gg(ix), hh(ix), pp(ix), uu(ix), &
            !                              ee(ix), phi2(ix)
	   write(unit=oxkx, fmt="(3f15.8,SP,256E24.14e3)") time, xx(ix), kx(mx), ff(ix,mx) 
          end do                                             
          !write(unit=oxkx, fmt=*) 
        end do                                             
        !write(unit=oxkx, fmt=*) 
        !write(unit=oxkx, fmt=*) 
 

        write(unit=oxtm, fmt=*) "## time, index = ", time, icount
        do ix = -nxx, nxx
          write(unit=oxtm, fmt="(2f15.8,SP,256E24.14e3)") time, xx(ix), ee(ix), phi2(ix), & 
                                                        gg(ix), hh(ix), pp(ix), uu(ix),   &
                                                        ee(ix), phi2(ix)         
        end do  
        write(unit=oxtm,fmt=*)         
  
        call flush(oxkx)
        call flush(oxtm)
 
      end if

! --- data output for integrated qauntities
      if (mod(loop+iout_0d,iout_0d) == 0) then

        ee(:)   = 0._DP
        phi2(:) = 0._DP
        do ix = -nxx, nxx
          do mx = -nkx, nkx
            ee(ix)   = ee(ix) + ff(ix,mx) * dkx
            phi2(ix) = phi2(ix) + ff(ix,mx)/(1._DP+kx(mx)**2+ky**2)**2 * dkx
          end do                                             
        end do                                             

        ee_total   = 0._DP
        phi2_total = 0._DP
        vy2_total  = 0._DP
        vp2_total  = 0._DP
        N02_total  = 0._DP
        Ns2_total  = 0._DP

        do ix = -nxx, nxx
          ee_total   = ee_total   + ee(ix)    * dx 
          phi2_total = phi2_total + phi2(ix)  * dx 
          vy2_total  = vy2_total  + gg(ix)**2 * dx 
          vp2_total  = vp2_total  + hh(ix)**2 * dx 
          N02_total  = N02_total  + pp(ix)**2 * dx 
          Ns2_total  = Ns2_total  + uu(ix)**2 * dx 
        end do                                             
                                              
        write(unit=oeng, fmt="(f15.8,SP,256E24.14e3)") time, ee_total, ee(0),   & 
                                                           phi2_total, phi2(0), & 
                                                            vy2_total, gg(0),   &
                                                            vp2_total, hh(0),   &
                                                            N02_total, pp(0),   &
                                                            Ns2_total, uu(0)
        call flush(oeng)

        print*, "OUTPUT! at loop, time = ", loop, time, ee_total, vy2_total

! --- check diverging behavior
        if ( ee(0) /= ee(0) ) then 
          print*, "NaN DIVERGE! at time =", time, "ee(0)=", ee(0)
          stop 
        end if

      end if

  end do

      close(oxkx)
      close(oxtm)
      close(oeng)
      print*, "DONE!!"
 
   stop 

END PROGRAM WK2Dext


!--------------------------------------
  SUBROUTINE calc_delta( ff, gg, hh, pp, uu, df, dg, dh, dpp, du )
!--------------------------------------
!     increment of delta-f within a time step

 use header

    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx,-nkx:nkx) :: ff

    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx) :: gg, hh, pp, uu

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx,-nkx:nkx) :: df

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx) :: dg, dh

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx) :: dpp, du

! --- local variables
    real(kind=DP), &
      dimension(-nxx:nxx,-nkx:nkx) :: dfdx, dfdk, ddfdx
    real(kind=DP), &
      dimension(-nxx:nxx) :: dpdx, dgdx
    real(kind=DP), &
      dimension(-nxx:nxx) :: Gammax, Pixy
    real(kind=DP), &
      dimension(-nxx:nxx) :: dgamdx2, dPidx
    real(kind=DP), &
      dimension(-nxx:nxx) :: dhdx2, dgdx2, dpdx2, dgdx3
    real(kind=DP), &
      dimension(-nxx:nxx,-nkx:nkx) :: dfdx2, dfdk2
    real(kind=DP), &
      dimension(-nxx:nxx,-nkx:nkx) :: domegadk, domegadx, gamma_L

    integer  ::  ix, mx

! ---  calculation of df = RHS * dt


      ! --- evaluate the derivatives 
      call dfdx_CFD6  ( ff, dfdx  )    ! dI/dx
      call dfdk_CFD6  ( ff, dfdk  )    ! dI/dk 
      call dfdx2_CFD4 ( ff, dfdx2 )    ! d2 I/dx2       
      call dfdk2_CFD4 ( ff, dfdk2 )    ! d2 I/dk2       

      call dfdx_CFD6_x  ( pp, dpdx  )  ! dN0p/dx 
      call dfdx_CFD6_x  ( gg, dgdx  )  ! dVy/dx 
      call dfdx2_CFD4_x ( gg, dgdx2 )  ! d2 Vy/dx2     
      call dfdx2_CFD4_x ( pp, dpdx2 )  ! d2 N0p/dx2 
      call dfdx2_CFD4_x ( hh, dhdx2 )  ! d2 V_para/dx2     
      call dfdx3_CFD4_x ( gg, dgdx3 )  ! d3 Vy/dx3 

             
      ! --- domega/dk, domega/dx, gamma_L
      do mx = -nkx, nkx               
        do ix = -nxx, nxx
          domegadk(ix,mx)  = -2._DP*kx(mx)*ky / ( 1._DP + kx(mx)**2 + ky**2 )**2            &
                                              * ( 1._DP - real(idg)*pp(ix) + real(ivx2)*dgdx2(ix) )

          domegadx(ix,mx)  = -ky/( 1._DP + kx(mx)**2 + ky**2 )*( real(idg)*dpdx(ix) - real(ivx2)*dgdx3(ix) )    &
                                              +  ky*( dgdx(ix) + S0 )

          !gamma_L(ix,mx)   =  ky*(kx(mx)**2+ky**2)/CC/( 1._DP + kx(mx)**2 + ky**2)**3      & 
          !                     * ( 1._DP - real(idg)*2._DP*pp(ix)      &
	  !			- real(ivx2)*(1._DP - kx(mx)**2 - ky**2)/(kx(mx)**2+ky**2)*dgdx2(ix) )
	  gamma_L(ix,mx)   =  ky*(kx(mx)**2+ky**2)/CC/( 1._DP + kx(mx)**2 + ky**2)**3*exp(-(kx(mx)-k0)**2/deltak**2)      & 
                                * ( 1._DP - real(idg)*2._DP*pp(ix)      &
				- real(ivx2)*(1._DP - kx(mx)**2 - ky**2)/(kx(mx)**2+ky**2)*dgdx2(ix) )
         ! gamma_L(ix,mx) = gamma0*exp(-(kx(mx)-k0)**2/deltak**2)
        end do  
      end do  

      ! --- calclation of Gamma_x and Pi_xy
      Gammax(:) = 0._DP
      Pixy(:)   = 0._DP
      do ix = -nxx, nxx
        do mx = -nkx, nkx
          Gammax(ix) = Gammax(ix) + dkx * ff(ix,mx)                                     &
                                  * ky**2/CC/( 1._DP + kx(mx)**2 + ky**2)**3            &
                                  * ( (kx(mx)**2+ky**2)*(1._DP - real(idg)*pp(ix)) - real(ivx2)*dgdx2(ix) ) 

          Pixy(ix) = Pixy(ix) + dkx * ff(ix,mx) * ( -kx(mx)*ky/( 1._DP + kx(mx)**2 + ky**2)**2 )           
        end do 
      end do 

      call dfdx2_CFD4_x ( Gammax, dgamdx2 )  ! d2 Gamma_x/dx2     
      call dfdx_CFD6_x  (   Pixy, dPidx   )  ! dPi_xy/dx 


      ! --- equation for I(x,kx)
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          df(ix,mx) = dt *  ( - domegadk(ix,mx)*dfdx(ix,mx) + domegadx(ix,mx)*dfdk(ix,mx)      & 
                              + gamma_L(ix,mx)*ff(ix,mx)    - delomg*ff(ix,mx)**2              &
                              + dk_nl*dfdk(ix,mx)**2        + dk_nl*ff(ix,mx)*dfdk2(ix,mx)     &
                              + mu_x*dfdx2(ix,mx)                                              &
                              + dk_l*dfdk2(ix,mx)                                              &
                            ) 

        end do
      end do

      ! --- equation for N0prime, Vy, Ns, V_para
      do ix = -nxx, nxx
        dpp(ix) = dt *  ( -dgamdx2(ix) + mu_n*dpdx2(ix) ) * real(idg)            ! for N0prime
        dg(ix)  = dt *  ( -2._DP*eps_n*uu(ix) - dPidx(ix) + mu_perp*dgdx2(ix) )  ! for Vy 
        du(ix)  = dt *  (  eps_n*gg(ix) + real(ivp)*eps_n*hh(ix)/q0 )            ! for Ns 
        dh(ix)  = dt *  ( -eps_n*uu(ix)/q0 + mu_para*dhdx2(ix) ) * real(ivp)     ! for V_para
      end do

    return
 
  END SUBROUTINE calc_delta


!--------------------------------------
  SUBROUTINE rkg4( ff, gg, hh, pp, uu, df, dg, dh, dpp, du, qf, qg, qh, qp, qu, istep )
!--------------------------------------
!     4th order Runge-Kutta-Gill

 use header

    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx,-nkx:nkx) :: ff
    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx) :: gg, hh, pp, uu
    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx, -nkx:nkx) :: df, qf
    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx) :: dg, qg
    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx) :: dh, qh
    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx) :: dpp, qp
    real(kind=DP), intent(inout), &
      dimension(-nxx:nxx) :: du, qu

    integer, intent(in) :: istep

    real(kind=DP) :: c1, c2, cq, c0


      if      ( istep == 1 ) then
        c1  =  0.5_DP
        c2  = -1._DP
        cq  = -2._DP
        c0  =  1._DP
      else if ( istep == 2 ) then
        c1  =  1._DP - sqrt( 0.5_DP )
        c2  = -c1
        cq  =  1._DP - 3._DP * c1
        c0  =  2._DP * c1
      else if ( istep == 3 ) then
        c1  =  1._DP + sqrt( 0.5_DP )
        c2  = -c1
        cq  =  1._DP - 3._DP * c1
        c0  =  2._DP * c1
      else if ( istep == 4 ) then
        c1  =  1._DP / 6._DP
        c2  = -1._DP / 3._DP
        cq  =  0._DP
        c0  =  0._DP
      end if

      do mx = -nkx, nkx
        do ix = -nxx, nxx
          ff(ix,mx) = ff(ix,mx)         &
                       + c1 * df(ix,mx) &
                       + c2 * qf(ix,mx)
          qf(ix,mx) = cq * qf(ix,mx)    &
                       + c0 * df(ix,mx)
        end do
      end do

      do ix = -nxx, nxx
        gg(ix) = gg(ix)                 &
                  + c1 * dg(ix)         &
                  + c2 * qg(ix)
        qg(ix) = cq * qg(ix)            &
                  + c0 * dg(ix)
        hh(ix) = hh(ix)                 &
                  + c1 * dh(ix)         &
                  + c2 * qh(ix)
        qh(ix) = cq * qh(ix)            &
                  + c0 * dh(ix)
        pp(ix) = pp(ix)                 &
                  + c1 * dpp(ix)        &
                  + c2 * qp(ix)
        qp(ix) = cq * qp(ix)            &
                  + c0 * dpp(ix)
        uu(ix) = uu(ix)                 &
                  + c1 * du(ix)         &
                  + c2 * qu(ix)
        qu(ix) = cq * qu(ix)            &
                  + c0 * du(ix)
      end do

  END SUBROUTINE rkg4


!--------------------------------------
  SUBROUTINE dfdx_CFD6( ff0, dfdx )
!--------------------------------------
!     6th order finite difference for df/dx

  use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx,-nkx:nkx) :: ff0

    real(kind=DP), &
      dimension(-nxx-3:nxx+3,-nkx:nkx) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx,-nkx:nkx) :: dfdx

    real(kind=DP) :: cef6
    integer  ::  ix, mx

      cef6  = 1._DP / ( 60._DP * dx )

! --- data copy & periodic B.C. 
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          ff(ix,mx) = ff0(ix,mx)
        end do
      end do
       
      do mx = -nkx, nkx
       ff(-nxx-1,mx) = ff0(nxx-1 ,mx)    
       ff(-nxx-2,mx) = ff0(nxx-2 ,mx)    
       ff(-nxx-3,mx) = ff0(nxx-3 ,mx)    
       ff(nxx+1 ,mx) = ff0(-nxx+1,mx)    
       ff(nxx+2 ,mx) = ff0(-nxx+2,mx)    
       ff(nxx+3 ,mx) = ff0(-nxx+3,mx)    
      end do


! --- calc. the derivative
      do mx = -nkx, nkx
        do ix = -nxx, nxx
                dfdx(ix,mx) =                                       &
                           ( +           ff(ix+3,mx)                &
                             -  9._DP  * ff(ix+2,mx)                &
                             +  45._DP * ff(ix+1,mx)                &
                             -  45._DP * ff(ix-1,mx)                &
                             +  9._DP  * ff(ix-2,mx)                &
                             -           ff(ix-3,mx)                &
                           ) * cef6
        end do
      end do

  return

  END SUBROUTINE dfdx_CFD6


!--------------------------------------
  SUBROUTINE dfdx_CFD6_x( ff0, dfdx )
!--------------------------------------
!     6th order finite difference for df/dx

  use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx) :: ff0

    real(kind=DP), &
      dimension(-nxx-3:nxx+3) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx) :: dfdx

    real(kind=DP) :: cef6
    integer  ::  ix, mx


      cef6  = 1._DP / ( 60._DP * dx )

! --- data copy & periodic B.C. 
      do ix = -nxx, nxx
        ff(ix) = ff0(ix)
      end do
         
      ff(-nxx-1) = ff0(nxx-1 )    
      ff(-nxx-2) = ff0(nxx-2 )    
      ff(-nxx-3) = ff0(nxx-3 )    
      ff(nxx+1 ) = ff0(-nxx+1)    
      ff(nxx+2 ) = ff0(-nxx+2)    
      ff(nxx+3 ) = ff0(-nxx+3)    

! --- calc. the derivative
      do ix = -nxx, nxx
              dfdx(ix) =                                       &
                         ( +           ff(ix+3)                &
                           -  9._DP  * ff(ix+2)                &
                           +  45._DP * ff(ix+1)                &
                           -  45._DP * ff(ix-1)                &
                           +  9._DP  * ff(ix-2)                &
                           -           ff(ix-3)                &
                         ) * cef6
      end do

  return

  END SUBROUTINE dfdx_CFD6_x



!--------------------------------------
  SUBROUTINE dfdk_CFD6( ff0, dfdk )
!--------------------------------------
!     6th order finite difference for df/dkx

  use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx,-nkx:nkx) :: ff0

    real(kind=DP), &
      dimension(-nxx:nxx,-nkx-3:nkx+3) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx,-nkx:nkx) :: dfdk

    real(kind=DP) :: cef6
    integer  ::  ix, mx 


      cef6  = 1._DP / ( 60._DP * dkx )

! --- data copy & periodic B.C. 
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          ff(ix,mx) = ff0(ix,mx)
        end do
      end do
       
      do ix = -nxx, nxx
       ff(ix,-nkx-1) = ff0(ix,nkx-1 )    
       ff(ix,-nkx-2) = ff0(ix,nkx-2 )    
       ff(ix,-nkx-3) = ff0(ix,nkx-3 )    
       ff(ix,nkx+1 ) = ff0(ix,-nkx+1)    
       ff(ix,nkx+2 ) = ff0(ix,-nkx+2)    
       ff(ix,nkx+3 ) = ff0(ix,-nkx+3)    
      end do

! --- calc. the derivative
      do mx = -nkx, nkx
        do ix = -nxx, nxx
                dfdk(ix,mx) =                                       &
                           ( +           ff(ix,mx+3)                &
                             -  9._DP  * ff(ix,mx+2)                &
                             +  45._DP * ff(ix,mx+1)                &
                             -  45._DP * ff(ix,mx-1)                &
                             +  9._DP  * ff(ix,mx-2)                &
                             -           ff(ix,mx-3)                &
                           ) * cef6
        end do
      end do

  return

  END SUBROUTINE dfdk_CFD6


!--------------------------------------
  SUBROUTINE dfdx2_CFD4( ff0, dfdx2 )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    calculation of df/dv_perp term
!
!    by M. Nunami and M. Nakata, April 2014
!
!-------------------------------------------------------------------------------
    use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx,-nkx:nkx) :: ff0

    real(kind=DP), &
      dimension(-nxx-2:nxx+2,-nkx:nkx) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx,-nkx:nkx) :: dfdx2

    real(kind=DP) :: cef4
    integer  ::  ix, mx


      cef4  = 1._DP / ( 12._DP * dx * dx )

! data copy & periodic B.C. 
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          ff(ix,mx) = ff0(ix,mx)
        end do
      end do
       
      do mx = -nkx, nkx
       ff(-nxx-1,mx) = ff0(nxx-1 ,mx)    
       ff(-nxx-2,mx) = ff0(nxx-2 ,mx)    
       ff(nxx+1 ,mx) = ff0(-nxx+1,mx)    
       ff(nxx+2 ,mx) = ff0(-nxx+2,mx)    
      end do

! --- calc. the derivative
      do mx = -nkx, nkx
        do ix = -nxx, nxx
                dfdx2(ix,mx) =                                      &
                           ( -           ff(ix+2,mx)                &
                             +  16._DP * ff(ix+1,mx)                &
                             -  30._DP * ff(ix  ,mx)                &
                             +  16._DP * ff(ix-1,mx)                &
                             -           ff(ix-2,mx)                &
                           ) * cef4
        end do
      end do


  return

  END SUBROUTINE dfdx2_CFD4


!--------------------------------------
  SUBROUTINE dfdx2_CFD4_x( ff0, dfdx2 )
!--------------------------------------
!     4th order finite difference for d^2f/dx^2

  use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx) :: ff0

    real(kind=DP), &
      dimension(-nxx-2:nxx+2) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx) :: dfdx2

    real(kind=DP) :: cef4
    integer  ::  ix, mx


      cef4  = 1._DP / ( 12._DP * dx * dx )

! --- data copy & periodic B.C. 
      do ix = -nxx, nxx
        ff(ix) = ff0(ix)
      end do
       
      ff(-nxx-1) = ff0(nxx-1 )    
      ff(-nxx-2) = ff0(nxx-2 )    
      ff(nxx+1 ) = ff0(-nxx+1)    
      ff(nxx+2 ) = ff0(-nxx+2)    

! --- calc. the derivative
      do ix = -nxx, nxx 
              dfdx2(ix) =                                      &
                         ( -           ff(ix+2)                &
                           +  16._DP * ff(ix+1)                &
                           -  30._DP * ff(ix  )                &
                           +  16._DP * ff(ix-1)                &
                           -           ff(ix-2)                &
                         ) * cef4
      end do

  return

  END SUBROUTINE dfdx2_CFD4_x


!--------------------------------------
  SUBROUTINE dfdx3_CFD4_x( ff0, dfdx3 )
!--------------------------------------
!     4th order finite difference for d^3f/dx^3

  use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx) :: ff0

    real(kind=DP), &
      dimension(-nxx-3:nxx+3) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx) :: dfdx3

    real(kind=DP) :: cef4
    integer  ::  ix, mx


      cef4  = 1._DP / ( 8._DP * dx * dx * dx )

! --- data copy & periodic B.C. 
      do ix = -nxx, nxx
        ff(ix) = ff0(ix)
      end do
       
      ff(-nxx-1) = ff0(nxx-1 )    
      ff(-nxx-2) = ff0(nxx-2 )    
      ff(-nxx-3) = ff0(nxx-3 )    
      ff(nxx+1 ) = ff0(-nxx+1)    
      ff(nxx+2 ) = ff0(-nxx+2)    
      ff(nxx+3 ) = ff0(-nxx+3)    

! --- calc. the derivative
      do ix = -nxx, nxx
              dfdx3(ix) =                                      &
                         ( -           ff(ix+3)                &
                           +  8._DP  * ff(ix+2)                &
                           -  13._DP * ff(ix+1)                &
                           +  13._DP * ff(ix-1)                &
                           -  8._DP  * ff(ix-2)                &
                           +           ff(ix-3)                &
                         ) * cef4
      end do

  return

  END SUBROUTINE dfdx3_CFD4_x


!--------------------------------------
  SUBROUTINE dfdk2_CFD4( ff0, dfdk2 )
!--------------------------------------
!     4th order finite difference for d^2f/dkx^2

  use header

    real(kind=DP), intent(in), &
      dimension(-nxx:nxx,-nkx:nkx) :: ff0

    real(kind=DP), &
      dimension(-nxx:nxx,-nkx-2:nkx+2) :: ff

    real(kind=DP), intent(out), &
      dimension(-nxx:nxx,-nkx:nkx) :: dfdk2

    real(kind=DP) :: cef4
    integer  ::  ix, mx

      cef4  = 1._DP / ( 12._DP * dkx * dkx )

! --- data copy & periodic B.C. 
      do mx = -nkx, nkx
        do ix = -nxx, nxx
          ff(ix,mx) = ff0(ix,mx)
        end do
      end do
       
      do ix = -nxx, nxx
       ff(ix,-nkx-1) = ff0(ix,nkx-1 )    
       ff(ix,-nkx-2) = ff0(ix,nkx-2 )    
       ff(ix,nkx+1 ) = ff0(ix,-nkx+1)    
       ff(ix,nkx+2 ) = ff0(ix,-nkx+2)    
      end do

! --- calc. the derivative
      do mx = -nkx, nkx
        do ix = -nxx, nxx
                dfdk2(ix,mx) =                                      &
                           ( -           ff(ix,mx+2)                &
                             +  16._DP * ff(ix,mx+1)                &
                             -  30._DP * ff(ix,mx  )                &
                             +  16._DP * ff(ix,mx-1)                &
                             -           ff(ix,mx-2)                &
                           ) * cef4
        end do
      end do

  return

  END SUBROUTINE dfdk2_CFD4

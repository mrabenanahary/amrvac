module mod_usr
  use mod_rho
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_types
  

  implicit none

  integer :: i_sol,i_err

  type(usr_params) :: usr_config

contains

  subroutine usr_init()
    use mod_variables

   !> to read and write user-defined parameters
   usr_set_parameters  => set_parameters

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => store_sol_err
    usr_print_log => print_min_max

    !> First, set user-defined parameters by default
    call usr_set_default_parameters

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    i_sol = var_set_extravar("solution", "solution")
    i_err = var_set_extravar("error", "error")
  end subroutine usr_init

  !> default user-defined parameters from files
  subroutine usr_set_default_parameters
    !-------------------------------------
    usr_config%parcel_method                   = 'fixed_numbers'
    usr_config%number_of_subregions            = 2 !< jet inlet + ambient medium
    
  end subroutine usr_set_default_parameters

  subroutine set_parameters
   use mod_variables
   use mod_obj_input_output
   implicit none
    
   !> Write parameters at the begining of each (re)start
   if(mype==0 .and. it==it_init) call usr_write_setting(usr_config)

  end subroutine set_parameters

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    ! initialize one grid 
    use mod_physics
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: rhoprofile(ixImin1:ixImax1,ixImin2:ixImax2)
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          if(ndim==2)then
             print *,'advection of VAC logo in 2D periodic box'
          else
             call mpistop("VAC logo advection is 2D, setup.pl -d=2")
          endif
       end if
       first=.false.
    end if
    call set_density_profile(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,0.0d0,x,rhoprofile)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=rhoprofile(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine initonegrid_usr

  subroutine set_density_profile(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,x,rhoprofile)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out) :: rhoprofile(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    logical::          maskv(ixImin1:ixImax1,ixImin2:ixImax2),&
       maska(ixImin1:ixImax1,ixImin2:ixImax2),maskc(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision:: rhoflat,rhosquare, xshift(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    double precision:: xc1,yc1,xa1,xa2,ya1,ya2,xb1,xb2,yb1,yb2,xc2,yc2, rad,&
       rad2,alp,nsig
    !----------------------------------------------------------------------------

    rhoflat  = 0.5d0 
    rhosquare= 2.0d0 
    
    
    
     ! account for periodicity, at least during one cycle....
     xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=x(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)-rho_v(1)*qt
     xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=x(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,2)-rho_v(2)*qt
     ! when v_x,v_y positive
     maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)<xprobmin1)
     where(maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)=x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-rho_v(1)*qt+(xprobmax1-xprobmin1)
     maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,2)<xprobmin2)
     where(maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,2)=x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-rho_v(2)*qt+(xprobmax2-xprobmin2)
     ! when v_x,v_y negative
     maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)>xprobmax1)
     where(maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)=x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-rho_v(1)*qt-(xprobmax1-xprobmin1)
     maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,2)>xprobmax2)
     where(maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)) xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,2)=x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-rho_v(2)*qt-(xprobmax2-xprobmin2)
     xc1=0.25d0
     yc1=0.50d0
     rad=0.23d0
     rad2=0.13d0
     alp=dpi/3.0d0
     xa1=xc1
     ya1=yc1-rad
     xa2=xc1-rad*cos(alp)
     ya2=yc1+rad*sin(alp)
     xb1=xa1
     yb1=ya1
     xb2=xc1+rad*cos(alp)
     yb2=yc1+rad*sin(alp)
     xc2=xc1
     yc2=ya2+sqrt(rad2**2-(xa2-xc2)**2)
     maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= ((xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)-xc1)**2+(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-yc1)**2 <= rad**2) .and.(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)>= (ya2-ya1)*(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xa1)/(xa2-xa1)+ya1) .and.(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)>= (yb2-yb1)*(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xb1)/(xb2-xb1)+yb1) .and.((xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xc2)**2+(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-yc2)**2 > rad2**2) 
     xc1=0.45d0
     yc1=0.475d0
     xa1=xc1
     ya1=yc1+rad
     xa2=xc1-rad*cos(alp)
     ya2=yc1-rad*sin(alp)
     xb1=xa1
     yb1=ya1
     xb2=xc1+rad*cos(alp)
     yb2=yc1-rad*sin(alp)
     xc2=xc1
     yc2=ya2-sqrt(rad2**2-(xa2-xc2)**2)
     maska(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= ((xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)-xc1)**2+(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-yc1)**2 <= rad**2) .and.(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)<= (ya2-ya1)*(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xa1)/(xa2-xa1)+ya1) .and.(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)<= (yb2-yb1)*(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xb1)/(xb2-xb1)+yb1) .and.((xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xc2)**2+(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-yc2)**2 > rad2**2) 
     xc1=0.75d0
     yc1=0.50d0
     alp=half*dpi-alp
     xa1=xc1-rad
     ya1=yc1
     xa2=xc1+rad*cos(alp)
     ya2=yc1+rad*sin(alp)
     xb1=xa1
     yb1=ya1
     xb2=xc1+rad*cos(alp)
     yb2=yc1-rad*sin(alp)
     yc2=yc1
     xc2=xa2+sqrt(rad2**2-(ya2-yc2)**2)
     maskc(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= ((xshift(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2,1)-xc1)**2+(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-yc1)**2 <= rad**2) .and.(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)<= (ya2-ya1)*(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xa1)/(xa2-xa1)+ya1) .and.(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)>= (yb2-yb1)*(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xb1)/(xb2-xb1)+yb1) .and.((xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1)-xc2)**2+(xshift(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        2)-yc2)**2 > rad2**2) 
     where(maskv(ixOmin1:ixOmax1,ixOmin2:ixOmax2).or.maska(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2).or.maskc(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        rhoprofile(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = rhosquare
     elsewhere
        rhoprofile(ixOmin1:ixOmax1,ixOmin2:ixOmax2)     = rhoflat
     endwhere
    

  end subroutine set_density_profile

  subroutine store_sol_err(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,w,x)
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: rhoprofile(ixImin1:ixImax1,ixImin2:ixImax2)

    call set_density_profile(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,qt,x,rhoprofile)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_sol) = rhoprofile(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_err) = dabs(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_sol))
  end subroutine store_sol_err

  subroutine print_min_max
    use mod_input_output, only: get_global_minima, get_global_maxima,&
        get_volume_average
    double precision   :: minvals(nw),maxvals(nw)

    integer :: iw
    double precision :: modes(nw,2), volume

    call get_global_minima(minvals)
    call get_global_maxima(maxvals)
    call get_volume_average(1,modes(:,1),volume)
    call get_volume_average(2,modes(:,2),volume)

    if (mype == 0) then
       write(*, "(A,4E16.8,A,3E16.8)") " time + rho min-max-tot:", global_time,&
           minvals(rho_), maxvals(rho_), modes(rho_,1), " Error: Linf-L1-L2:",&
           maxvals(i_err),modes(i_err,1),dsqrt(modes(i_err,2))
    end if
  end subroutine print_min_max

end module mod_usr


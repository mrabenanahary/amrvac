module mod_obj_geometry
  use mod_global_parameters
  use mod_geometry
  use mod_obj_global_parameters
  use mod_obj_types 
  use mod_obj_error_handling 
  use mod_obj_read_parameters
 
  implicit none                                            

  type(error_handler) :: obj_geometry_error
  
  contains



  function geom_to_coord(geom_var) result(coord_var)
    implicit none
    character(len=*),intent(in) :: geom_var
    integer                     :: coord_var
    !--------------------------------------------
    select case (geom_var)
     case ("Cartesian","Cartesian_1D","Cartesian_2D","Cartesian_3D",&
           "Cartesian_1.5D","Cartesian_1.75D","Cartesian_2.5D")
       coord_var=cartesian
     case ("Cartesian_1D_expansion")
       coord_var=cartesian_expansion
     case ("cylindrical","cylindrical_2D","cylindrical_3D",&
           "cylindrical_2.5D")
       coord_var=cylindrical
     case ("polar","polar_2D","polar_3D",&
           "polar_1.5D","polar_2.5D")
       coord_var=polar
     case ("spherical","spherical_2D","spherical_3D",&
           "spherical_2.5D")
       coord_var=spherical
     case default
       call mpistop("Unknown geometry input specified")
    end select

  end function geom_to_coord

  !> Computes the idim component (in coord_output coordinates) 
  !> of the v vector expressed (in coord_input coordinates)
  !> from the spatial coordinates x (in any coordinates)
  !> in any of the n_coord =1,2,3 dimensions
  !>
  !> Optional : check_input,check_output check that the input and output
  !>            vector components coordinates may be converted in AMRVAC
  function vector_comp(ixI^L,ixO^L,v,idim,coord_input,coord_output,&
                       x,x_coord,n_coord,&
                       check_input,check_output) result(v_out)   
    implicit none
    integer, intent(in)              :: ixI^L,ixO^L,n_coord,idim
    character(len=*),intent(in)      :: x_coord,coord_input,coord_output
    real(dp), intent(in)             :: x(ixI^S,1:n_coord)
    real(dp), intent(in)             :: v(ixI^S,1:n_coord)
    real(dp)                         :: v_out(ixI^S) 
    logical, intent(inout), optional :: check_input
    logical, intent(inout), optional :: check_output
    !--local---------
    integer                          :: icoord,ocoord,xcoord, i1
    logical, dimension(2)            :: check_list
    integer, dimension(2)            :: coord_list
    real(dp)                         :: radius(ixI^S),theta(ixI^S),phi(ixI^S)
    !-------------------------------------------------------

    if(.not.present(check_input)) check_input = .false. !< default to .false.
    if(.not.present(check_output)) check_output = .false. !< default to .false.


    icoord = geom_to_coord(coord_input)
    ocoord = geom_to_coord(coord_output)
    xcoord = geom_to_coord(x_coord)

    !> Trivial case : return the same components
    !> TO DO: check for Cartesian expansion/stretched case
    if(icoord==ocoord)then
      v_out = 1.0_dp
      return
    end if
    
    select case(n_coord)
     !> 1D cases
     case(1)
      
      !> check for available input and output coordinates
      check_list = [check_input,check_output]
      coord_list = [icoord,ocoord]
      do i1=1,2
        if(check_list(i1)) call obj_geometry_error%check_integer_in_list(coord_list(i1),&
        [cartesian,cartesian_expansion,cartesian_stretched,polar],&
        4)
      end do

      !> Start from input coordinates icoord to output coordinates ocoord
      select case(icoord)
       case(cartesian,cartesian_expansion,cartesian_stretched)
        ! TO DO: check for Cartesian expansion/stretched case
        select case(ocoord)
         case(polar) !< Cartesian -> Polar
          v_out = 1.0_dp
        end select
       case(polar)        
        select case(ocoord)
         case(cartesian,cartesian_expansion,cartesian_stretched) !< Polar -> Cartesian
          v_out = 1.0_dp
        end select

      end select
      
     !> 2D cases
     case(2)   

      !> all input and output coordinates avalaible so : no check needed

      !> Start from input coordinates icoord to output coordinates ocoord
      select case(icoord)
       case(spherical)

        select case(ocoord)
         case(cartesian,cartesian_expansion,&
         cartesian_stretched,cylindrical) !< Spherical -> Cartesian or cylindrical

          !TO DO: test those subroutines
          call compute_angle_theta(ixI^L,ixO^L,x,theta,xcoord,n_coord)                 

          ! TO DO: check for Cartesian expansion/stretched case
          !> (vec(V))_x = V_r*sin(theta)+V_theta*cos(theta)
          !> (vec(V))_rho = V_r*sin(theta)+V_theta*cos(theta)
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1) * dsin(theta(ixO^S))&
                        + v(ixO^S,2) * dcos(theta(ixO^S))
          !> (vec(V))_y = V_r*cos(theta)-V_theta*sin(theta) 
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,1) * dcos(theta(ixO^S))&
                        - v(ixO^S,2) * dsin(theta(ixO^S))

         case(polar) !< Spherical -> polar
          !> (vec(V))_r = V_r
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1)
          !> (vec(V))_phi = -V_theta 
          if(idim==2)&
           v_out(ixO^S) = - v(ixO^S,2)                         
        end select

       case(cylindrical)

        select case(ocoord)
         case(cartesian,cartesian_expansion,&
         cartesian_stretched) !< cylindrical -> Cartesian

          ! TO DO: check for Cartesian expansion/stretched case
          !> (vec(V))_x = V_rho
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1)
          !> (vec(V))_y = V_y 
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,2)

         case(polar) !< cylindrical -> polar

          !TO DO: test those subroutines
          call compute_angle_phi(ixI^L,ixO^L,x,phi,xcoord,n_coord)  

          !> (vec(V))_r = V_rho*cos(phi)-V_phi*sin(phi)
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1)*dcos(phi(ixO^S))&
                        - v(ixO^S,2)*dsin(phi(ixO^S))
          !> (vec(V))_phi = V_rho*sin(phi)+V_phi*cos(phi)
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,1)*dsin(phi(ixO^S))&
                        + v(ixO^S,2)*dcos(phi(ixO^S)) 

         case(spherical) !< cylindrical -> spherical

          !TO DO: test those subroutines
          call compute_angle_theta(ixI^L,ixO^L,x,theta,xcoord,n_coord)                 

          !> (vec(V))_r = V_rho*sin(theta)+V_y*cos(theta)
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1)*dsin(theta(ixO^S))&
                        + v(ixO^S,2)*dcos(theta(ixO^S))
          !> (vec(V))_theta = V_rho*cos(theta)-V_y*sin(theta)
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,1)*dcos(theta(ixO^S))&
                        - v(ixO^S,2)*dsin(theta(ixO^S))        

        end select

       case(cartesian,cartesian_expansion,cartesian_stretched)
       ! TO DO: check for Cartesian expansion/stretched case

        select case(ocoord)
         case(spherical) !< Cartesian -> spherical

          !TO DO: test those subroutines
          call compute_angle_theta(ixI^L,ixO^L,x,theta,xcoord,n_coord)

          !> (vec(V))_r = V_x*sin(theta) + V_y*cos(theta)
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1) * dsin(theta(ixO^S))&
                        + v(ixO^S,2) * dcos(theta(ixO^S))
          !> (vec(V))_theta = V_x*cos(theta) - V_y*sin(theta)
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,1) * dcos(theta(ixO^S))&
                        - v(ixO^S,2) * dsin(theta(ixO^S))

         case(cylindrical) !< Cartesian -> cylindrical

          !> (vec(V))_rho = V_x
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1)
          !> (vec(V))_y = V_y
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,2)

         case(polar) !< Cartesian -> polar

          !TO DO: test those subroutines
          call compute_angle_phi(ixI^L,ixO^L,x,phi,xcoord,n_coord)   

          !> (vec(V))_r = V_x*cos(phi)-V_y*sin(phi)
          if(idim==1)&
           v_out(ixO^S) = v(ixO^S,1)*dcos(phi(ixO^S))&
                        - v(ixO^S,2)*dsin(phi(ixO^S))
          !> (vec(V))_phi = V_x*sin(phi)+V_y*cos(phi)
          if(idim==2)&
           v_out(ixO^S) = v(ixO^S,1)*dsin(phi(ixO^S))&
                        + v(ixO^S,2)*dcos(phi(ixO^S))        

        end select

      end select        


     case default
      call mpistop('Unknown number of coords')
    end select

  end function vector_comp

  subroutine compute_angle_theta(ixI^L,ixO^L,x,theta,coord,n_coord) 
    implicit none
    integer, intent(in)              :: ixI^L,ixO^L
    real(dp), intent(in)             :: x(ixI^S,1:n_coord)
    integer,intent(in)               :: coord, n_coord
    real(dp), intent(out)            :: theta(ixI^S)
    !--local------------
    logical                          :: dummy_error(ixI^S)
    !---------------------------------------------------

    dummy_error = .false.

    select case(n_coord)
     case(1) !< 1D cases
      select case(coord)   
       case(Cartesian,cartesian_expansion,&
       cartesian_stretched,spherical,polar,cylindrical) 
        !> arbitrarily treat spherical/polar/cylindrical
        !> as what is for cartesian case in 1D-axis
        where(x(ixO^S,1)>=0.0_dp)
          theta(ixO^S) = 0.0_dp
        elsewhere
          theta(ixO^S) = dpi
        end where
      end select
     case(2) !< 2D cases
      select case(coord)
        case(spherical)

          where(dabs(x(ixO^S,1))<smalldouble * usr_unit_length)
            theta(ixO^S) = 0.0_dp !< (conventionally arbitrary)  
          elsewhere
            theta(ixO^S) = x(ixO^S,2) !< theta in spherical
          end where

        case(Cartesian,cartesian_expansion,cartesian_stretched)

          where((dabs(x(ixO^S,1))<smalldouble * usr_unit_length)&
          .and. (dabs(x(ixO^S,2))<smalldouble * usr_unit_length))
            !> theta = 0 if x=y=0 (conventionally arbitrary)
            theta(ixO^S) = 0.0_dp
          elsewhere
            !> theta = arcos(y/sqrt(x^2+y^2)) in [0,pi] forall (x,y) | x*y /= 0
            theta(ixO^S) = dacos(x(ixO^S,2) / &
            dsqrt(x(ixO^S,1) ** 2.0_dp + x(ixO^S,2) ** 2.0_dp))
          end where

        case(polar)

          where(      x(ixO^S,1) <  smalldouble * usr_unit_length)
            !> theta = 0 if r=0 (conventionally arbitrary)
            theta(ixO^S) = 0.0_dp
          else where((x(ixO^S,2) >= 0.0_dp     )&
              .and.  (x(ixO^S,2) <= half * dpi ))
            !> (x>=0 , y>=0) = (phi in [0,pi/2]) : 
            !> theta = pi/2 - phi in [0,pi/2]
            theta(ixO^S) = half * dpi - x(ixO^S,2)
          else where((x(ixO^S,2) >= 3.0_dp * half * dpi )&
               .and. (x(ixO^S,2) <  2.0_dp * dpi        ))
            !> (x>=0 , y<0) = (phi in [3pi/2,2*pi[) : 
            !> theta = 5pi/2 - phi in ]pi/2,pi]
            theta(ixO^S) = 5.0_dp * half * dpi - x(ixO^S,2)
          else where((x(ixO^S,2) >  half * dpi )&
               .and. (x(ixO^S,2) <= dpi        ))
            !> (x<0 , y>=0) = (phi in  ]pi/2,pi]) : 
            !> theta = phi - pi/2 in ]0,pi/2]
            theta(ixO^S) = x(ixO^S,2) - half * dpi
          else where((x(ixO^S,2) >  dpi                 )&
               .and. (x(ixO^S,2) <  3.0_dp * half * dpi )) 
            !> (x<0 , y<0) = (phi in ]pi,3pi/2[) : 
            !> theta = phi - pi/2 in ]pi/2,pi[]
            theta(ixO^S) = x(ixO^S,2) - half * dpi
          elsewhere
            dummy_error(ixO^S) = .true.
          end where

          !> if x contains some fishy untreated/unreal value !
          if(any(dummy_error))&
            call mpistop('Something fishy about conversion theta(polar->spherical) !')

        case(cylindrical)

          where(x(ixO^S,1) < smalldouble * usr_unit_length)
            !> theta = 0 if rho=0 (conventionally arbitrary)
            theta(ixO^S) = 0.0_dp
          elsewhere
           !> theta = acos(y^2/sqrt(y^2+rho^2) for all (rho,y) in plane
           theta(ixO^S) = dacos(x(ixO^S,2) /&
            sqrt(x(ixO^S,1) ** 2.0_dp + x(ixO^S,2) ** 2.0_dp) )
          end where

      end select
     case(3) !< 3D cases : surprisingly simpler
      select case(coord)   
       case(spherical)
        theta(ixO^S) = x(ixO^S,2)
       case(cartesian,cartesian_expansion,cartesian_stretched)

        where(dabs(x(ixO^S,1)) < smalldouble * usr_unit_length)
         !> theta = 0 if r=0 (conventionally arbitrary)
         theta(ixO^S) = 0.0_dp
        elsewhere
         !> theta = acos(z/sqrt(x^2+y^2+z^2))
         theta(ixO^S) = dacos(x(ixO^S,3) /&
          sqrt(x(ixO^S,1) ** 2.0_dp + &
               x(ixO^S,2) ** 2.0_dp + &
               x(ixO^S,3) ** 2.0_dp   ))
        end where

       case(cylindrical)

        where(dabs(x(ixO^S,1)) < smalldouble * usr_unit_length)
         !> theta = 0 if rho=0 (conventionally arbitrary)
         theta(ixO^S) = 0.0_dp
        elsewhere       
         !> theta = acos(z/sqrt(rho^2+z^2))
         theta(ixO^S) = dacos(x(ixO^S,2) /&
          sqrt(x(ixO^S,1) ** 2.0_dp + &
               x(ixO^S,2) ** 2.0_dp   ))
        end where

       case(polar)

        where(dabs(x(ixO^S,1)) < smalldouble * usr_unit_length)
         !> theta = 0 if r=0 (conventionally arbitrary)
         theta(ixO^S) = 0.0_dp
        elsewhere       
         !> theta = acos(z/sqrt(r^2+z^2))
         theta(ixO^S) = dacos(x(ixO^S,3) /&
           sqrt(x(ixO^S,1) ** 2.0_dp + &
                x(ixO^S,3) ** 2.0_dp   ))               
        end where

      end select
    end select





  end subroutine compute_angle_theta

  subroutine compute_angle_phi(ixI^L,ixO^L,x,phi,coord,n_coord) 
    implicit none
    integer, intent(in)              :: ixI^L,ixO^L
    real(dp), intent(in)             :: x(ixI^S,1:n_coord)
    integer,intent(in)               :: coord, n_coord
    real(dp), intent(out)            :: phi(ixI^S)
    !--local------------
    logical                          :: dummy_error(ixI^S)
    !---------------------------------------------------

    dummy_error = .false.

    select case(n_coord)
     case(1) !< 1D cases
      !> arbitrarily set to 0
      phi(ixO^S) = 0.0_dp
     case(2) !< 2D cases
      select case(coord)
        case(spherical)

          where(      dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length)
            !> phi = 0 if r=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where( dabs(x(ixO^S,2)) <  smalldouble)
            phi(ixO^S) = half * dpi
          else where( dabs(x(ixO^S,2)-dpi) <  smalldouble)
            phi(ixO^S) = 3.0_dp * half * dpi            
          else where((x(ixO^S,2) >= 0.0_dp     )&
               .and. (x(ixO^S,2) <= half * dpi ))
            !> (x>=0 , y>=0) = (theta in [0,pi/2]) : 
            !> phi = pi/2 - theta in [0,pi/2]
            phi(ixO^S) = half * dpi - x(ixO^S,2)
          else where((x(ixO^S,2) >= - half * dpi )&
               .and. (x(ixO^S,2) <   0.0_dp      ))
            !> (x>=0 , y<0) = (theta in [-pi/2,0[) : 
            !> phi = pi/2 - theta in ]pi/2,pi]
            phi(ixO^S) = half * dpi - x(ixO^S,2)
          else where((x(ixO^S,2) >  half * dpi )&
               .and. (x(ixO^S,2) <= dpi        ))
            !> (x<0 , y>=0) = (theta in  ]-pi,-pi/2]) : 
            !> phi = pi/2 - theta in ]pi,3pi/2]
            phi(ixO^S) = half * dpi - x(ixO^S,2)
          else where((x(ixO^S,2) >  half * dpi          )&
               .and. (x(ixO^S,2) <  dpi                 )) 
            !> (x<0 , y<0) = (theta in ]pi/2,pi[) : 
            !> phi = 5pi/2 - theta in ]3pi/2,2pi[]
            phi(ixO^S) = 5.0_dp * half * dpi - x(ixO^S,2)
          elsewhere
            dummy_error(ixO^S) = .true.
          end where

          !> if x contains some fishy untreated/unreal value !
          if(any(dummy_error))&
            call mpistop('Something fishy about conversion theta(spherical)->phi(spherical) in 2D !')


        case(Cartesian,cartesian_expansion,cartesian_stretched)

          where((dabs(x(ixO^S,1))<smalldouble * usr_unit_length)&
          .and. (dabs(x(ixO^S,2))<smalldouble * usr_unit_length))
            !> phi = 0 if x=y=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where( (dabs(x(ixO^S,1)) <  smalldouble ) &
          .and.       (x(ixO^S,2))      >  0.0_dp      )
            phi(ixO^S) = half * dpi
          else where( (dabs(x(ixO^S,1)) <  smalldouble ) &
          .and.       (x(ixO^S,2))      <  0.0_dp      )
            phi(ixO^S) = 3.0_dp * half * dpi              
          else where((x(ixO^S,2) >= 0.0_dp     ))
            !> (y>=0) : 
            !> phi = arcos(x/sqrt(x^2+y^2)) in [0,pi]
            phi(ixO^S) = dacos(x(ixO^S,1)/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          else where((x(ixO^S,2) < 0.0_dp     ))
            !> (y<0) : 
            !> phi = 2pi - arcos(x/sqrt(x^2+y^2)) in ]pi,2pi[
            phi(ixO^S) = 2.0_dp * dpi - &
            dacos(x(ixO^S,1)/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          elsewhere
            dummy_error(ixO^S) = .true.
          end where

          !> if x contains some fishy untreated/unreal value !
          if(any(dummy_error))&
            call mpistop('Something fishy about conversion phi(cartesian->spherical) in 2D !')


        case(polar)

          where(      dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length)
            !> phi = 0 if r=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where(      dabs(x(ixO^S,1)) >=  smalldouble * usr_unit_length)
            !> (r/=0) 
            !> phi in [0,2*pi]
            phi(ixO^S) = x(ixO^S,2)
          elsewhere
            dummy_error(ixO^S) = .true.
          end where

          !> if x contains some fishy untreated/unreal value !
          if(any(dummy_error))&
            call mpistop('Something fishy about conversion phi(polar->spherical) in 2D !')

        case(cylindrical)

          where((dabs(x(ixO^S,1)) < smalldouble * usr_unit_length )&
          .and.(dabs(x(ixO^S,2))  < smalldouble * usr_unit_length ))
            !> theta = 0 if rho=y=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where( (dabs(x(ixO^S,1)) <  smalldouble ) &
               .and.  (x(ixO^S,2))      >  0.0_dp      )
            phi(ixO^S) = half * dpi
          else where( (dabs(x(ixO^S,1)) <  smalldouble ) &
               .and.  (x(ixO^S,2))      <  0.0_dp      )
            phi(ixO^S) = 3.0_dp * half * dpi                 
          else where((x(ixO^S,1) >= 0.0_dp     )&
               .and. (x(ixO^S,2) >  0.0_dp     ))
            !> (rho>=0 , y>0) : 
            !> phi = arcos(rho/sqrt(rho^2+y^2)) in [0,pi/2]
            phi(ixO^S) = dacos(x(ixO^S,1)/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          else where((x(ixO^S,1) <  0.0_dp     )&
               .and. (x(ixO^S,2) >= 0.0_dp     ))
            !> (rho<0 , y>=0) : 
            !> phi = pi - arcos(abs(rho)/sqrt(rho^2+y^2)) in ]pi/2,pi]
            phi(ixO^S) = dpi - dacos(dabs(x(ixO^S,1))/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          else where((x(ixO^S,1) <= 0.0_dp     )&
               .and. (x(ixO^S,2) <  0.0_dp     ))
            !> (rho<=0 , y<0) : 
            !> phi = pi + arcos(abs(rho)/sqrt(rho^2+y^2)) in ]pi,3pi/2]
            phi(ixO^S) = dpi + dacos(dabs(x(ixO^S,1))/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          else where((x(ixO^S,1) >  0.0_dp     )&
               .and. (x(ixO^S,2) <  0.0_dp     ))
            !> (rho>0 , y<0) : 
            !> phi = 2pi - arcos(abs(rho)/sqrt(rho^2+y^2)) in ]pi,3pi/2]
            phi(ixO^S) = 2.0_dp * dpi - dacos(x(ixO^S,1)/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          elsewhere
            dummy_error(ixO^S) = .true.
          end where

          !> if x contains some fishy untreated/unreal value !
          if(any(dummy_error))&
            call mpistop('Something fishy about conversion phi(cylindrical->spherical) !')
          

      end select
     case(3) !< 3D cases : surprisingly simpler
      select case(coord)   
       case(spherical)
          where(      dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length)
            !> phi = 0 if r=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where((dabs(x(ixO^S,1))     >=  smalldouble * usr_unit_length) &
          .and.       ( dabs(x(ixO^S,2))   <   smalldouble))
            phi(ixO^S) = half * dpi
          else where((dabs(x(ixO^S,1))       >=  smalldouble * usr_unit_length) &
          .and.      (dabs(x(ixO^S,2)-dpi)   <   smalldouble))
            phi(ixO^S) = 3.0_dp * half * dpi    
          elsewhere
            phi(ixO^S) = x(ixO^S,3) !< phi in spherical
          end where       
       case(cartesian,cartesian_expansion,cartesian_stretched)
          where(      (dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length) &
          .and.       (dabs(x(ixO^S,2)) <  smalldouble * usr_unit_length) &
          .and.       (dabs(x(ixO^S,3)) <  smalldouble * usr_unit_length))
            !> phi = 0 if x=y=z=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where((dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length) &
          .and.      (dabs(x(ixO^S,2)) <  smalldouble * usr_unit_length) &
          .and.      (     x(ixO^S,2)  >  0.0_dp))
            phi(ixO^S) = half * dpi
          else where((dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length) &
          .and.      (dabs(x(ixO^S,2)) <  smalldouble * usr_unit_length) &
          .and.      (     x(ixO^S,2)  <  0.0_dp))
            phi(ixO^S) = 3.0_dp * half * dpi            
          else where(      x(ixO^S,2)  >= 0.0_dp )
            !> (y>=0) : 
            !> phi = arcos(x/sqrt(x^2+y^2)) in [0,pi]
            phi(ixO^S) = dacos(x(ixO^S,1)/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          else where(      x(ixO^S,2)  <  0.0_dp )
            !> (y<0) : 
            !> phi = 2pi - arcos(x/sqrt(x^2+y^2)) in [pi,2pi[
            phi(ixO^S) = 2.0_dp * dpi - dacos(x(ixO^S,1)/sqrt(&
            x(ixO^S,1) ** 2.0_dp + &
            x(ixO^S,2) ** 2.0_dp))
          elsewhere
            dummy_error(ixO^S) = .true.
          end where

          !> if x contains some fishy untreated/unreal value !
          if(any(dummy_error))&
            call mpistop('Something fishy about conversion phi(cartesian)->phi(spherical) in 2D !')



       case(cylindrical)
          where(      (dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length) &
          .and.       (dabs(x(ixO^S,2)) <  smalldouble * usr_unit_length))
            !> phi = 0 if rho=z=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where((dabs(x(ixO^S,1))    <=  smalldouble * usr_unit_length) &
          .and.      (     x(ixO^S,2)     >=  0.0_dp))
            phi(ixO^S) = half * dpi
          else where((dabs(x(ixO^S,1))    <=  smalldouble * usr_unit_length) &
          .and.      (     x(ixO^S,2)     <=  0.0_dp))
            phi(ixO^S) = 3.0_dp * half * dpi    
          elsewhere
            phi(ixO^S) = x(ixO^S,3) !< phi in cylindrical
          end where       


       case(polar)
          where(      (dabs(x(ixO^S,1)) <  smalldouble * usr_unit_length) &
          .and.       (dabs(x(ixO^S,2)) <  smalldouble * usr_unit_length))
            !> phi = 0 if r=z=0 (conventionally arbitrary)
            phi(ixO^S) = 0.0_dp
          else where((dabs(x(ixO^S,1))     <   smalldouble * usr_unit_length) &
          .and.      (     x(ixO^S,2)      >   0.0_dp))
            phi(ixO^S) = half * dpi
          else where((dabs(x(ixO^S,1))     <   smalldouble * usr_unit_length) &
          .and.      (     x(ixO^S,2)      <   0.0_dp))
            phi(ixO^S) = 3.0_dp * half * dpi    
          elsewhere
            phi(ixO^S) = x(ixO^S,2) !< phi in polar
          end where       

      end select
    end select

  end subroutine compute_angle_phi  

end module mod_obj_geometry

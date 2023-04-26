module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => rm1d_init_one_grid

        usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 

    call set_coordinate_system("polar_1.5D")
    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    ! iprob==1 rarefaction wave & shock
    if (iprob==1) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.5d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 1.5d0
        elsewhere
           w(ix^S,rho_)   = 0.5d0
           w(ix^S,mom(1)) = 0.0d0
           w(ix^S,e_)     = 0.5d0
        end where
    ! iprob==2  shock & shock
    else if (iprob==2) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 2.0d0
           w(ix^S,e_)     = 1.5d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -2.0d0
           w(ix^S,e_)     = 1.5d0
        end where
    ! iprob==3  rarefaction wave
    else if (iprob==3) then
        where (abs(x(ix^S,1))<0.5d0)
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = -0.5d0
           w(ix^S,e_)     = 1.0d0
        elsewhere
           w(ix^S,rho_)   = 1.0d0
           w(ix^S,mom(1)) = 0.5d0
           w(ix^S,e_)     = 1.0d0
        end where
    else
        call mpistop("iprob not available!")
    end if

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm1d_init_one_grid


    subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    ! .. local ..
    integer                    :: iw
    integer, parameter         :: iradius             = 1
    integer, parameter         :: itheta      = 2
    integer, parameter         :: iphi            = 3
    !----------------------------------------------------

    varnames='ir iphi iz'



  end subroutine specialvarnames_output

    !> special output
  subroutine specialvar_output(ixI^L,ixO^L,win,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
  !/
    use mod_physics
    use mod_dust
    implicit none
    integer, intent(in)        :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: win(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    ! .. local ..
    double precision                   :: w(ixI^S,nw)
    integer, parameter         :: iradius             = 1
    integer, parameter         :: itheta      = 2
    integer, parameter         :: iphi            = 3
    double precision                 :: mp,kB,me
    integer :: iw
    !-------------------------------------------------------


    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      me = const_me*1.0d-3
    else
      mp=mp_cgs
      kB=kB_cgs
      me = const_me
    end if


    w(ixI^S,1:nw) = win(ixI^S,1:nw)

    Loop_iw :  do iw = 1,nwauxio
    select case(iw)        
      case(iradius)
        normconv(nw+iradius)     = 1.0d0
        win(ixO^S,nw+iradius)    = x(ixO^S,1)*unit_length
      case(itheta)
        normconv(nw+itheta)     = 1.0d0
        win(ixO^S,nw+itheta)    = x(ixO^S,2)*unit_length
      case(iphi)
        normconv(nw+iphi)     = 1.0d0
        win(ixO^S,nw+iphi)    = x(ixO^S,3)*unit_length        
      case default
       write(*,*)'is not implimented at specialvar_output in mod_user'
       call mpistop('the code stops!')
    end select
  end do Loop_iw
    !----------------------------------------------------

   !w(ixI^S,1:nw)=win(ixI^S,1:nw)




  end subroutine specialvar_output




end module mod_usr

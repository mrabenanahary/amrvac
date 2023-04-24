! blast wave
module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0        => specialset_B0
    usr_init_vector_potential=>initvecpot_usr

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rbs,xc^D,Bloc(ixI^S,1:ndir)
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          print *,'MHD blast wave in Cartesian coordinate'
       end if
       first=.false.
    end if
    w(ixO^S,rho_)=1.d0
    w(ixO^S,p_)=1.d0
    rbs=0.2d0
    {xc^D=(xprobmin^D+xprobmax^D)*0.5d0\}
    where(^D&(x(ixO^S,^D)-xc^D)**2+ <rbs**2)
      w(ixO^S,p_)=100.d0
    endwhere
    if(B0field) then
      w(ixO^S,mag(:))=0.d0
    else if(stagger_grid) then
      call b_from_vector_potential(block%ixGs^L,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
    else
      call get_B(ixI^L,ixO^L,Bloc,x)
      w(ixO^S,mag(:))=Bloc(ixO^S,:)
    end if
    w(ixO^S,mom(:))=0.d0

    if(mhd_glm) w(ixO^S,psi_)=0.d0

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    use mod_global_parameters
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)

    if (idir==3) then
      A(ixC^S) = Busr*xC(ixC^S,2)
    else if(idir==2) then 
      A(ixC^S) = Busr*xC(ixC^S,1)
    else
      A(ixC^S) = Busr*xC(ixC^S,3)
    end if

  end subroutine initvecpot_usr

  subroutine get_B(ixI^L,ixO^L,B,x)
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out) :: B(ixI^S,1:ndir)

    B(ixO^S,:)=Busr

  end subroutine get_B


    !> this subroutine is ONLY to be used for computing auxiliary variables
  !> which happen to be non-local (like div v), and are in no way used for
  !> flux computations. As auxiliaries, they are also not advanced

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    ! .. local ..
    integer                    :: iw
    integer, parameter         :: iradius             = 1
    integer, parameter         :: itheta      = 2
    integer, parameter         :: iphi            = 3
    !----------------------------------------------------

    varnames='ix iy iz'



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

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here one can add a steady (time-independent) potential background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    double precision :: Bloc(ixI^S,1:ndir)

    call get_B(ixI^L,ixO^L,Bloc,x)
    wB0(ixO^S,:)=wB0(ixO^S,:)+Bloc(ixO^S,:)

  end subroutine specialset_B0

end module mod_usr

module mod_obj_types
  use mod_global_parameters

  implicit none

   !> User-defined parameters 

   type usr_params

    character(len=name_len) :: is_it_null
    procedure(n_arg), pointer :: set_value => null()
    

    contains

    procedure, pass(self) :: set_default_proc => set_default_procdr

   end type usr_params

contains

  subroutine n_arg(self)
    implicit none
    class(usr_params) :: self
    !----------------------------------
  end subroutine n_arg

    function check_null_proc(self,proc_name) result(is_it_null)
    implicit none
    class(usr_params) :: self
    character(len=*), intent(in) :: proc_name
    logical :: is_it_null
    !----------------------------------

    is_it_null = .false.

    select case(trim(proc_name))
    case('set_value')
      if(.not.associated(self%set_value))then
        is_it_null = .true.
      end if
    end select

  end function check_null_proc

  subroutine set_default_procdr(self,proc_name)
    implicit none
    class(usr_params) :: self
    character(len=*), intent(in) :: proc_name
    !----------------------------------

    if(check_null_proc(self,proc_name))self%is_it_null = 'yesyesyes'

  end subroutine set_default_procdr

end module mod_obj_types


module mod_obj_read_parameters
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_types  

  implicit none



contains

  subroutine new_procedure(self)
    implicit none
    class(usr_params) :: self
    !----------------------------------
    self%is_it_null = 'nonono'
  end subroutine new_procedure

  subroutine change_procedure(usr_config_var)
   implicit none
   type(usr_params), intent(inout) :: usr_config_var
   !-----------------------------------------------------
   usr_config_var%set_value => new_procedure
  end subroutine change_procedure

end module mod_obj_read_parameters


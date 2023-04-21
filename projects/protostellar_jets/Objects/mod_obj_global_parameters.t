module mod_obj_global_parameters
  use mod_global_parameters
  use mod_constants
  use mod_physics
 
  implicit none

  !> Global user-defined parameters for mod_usr 
  integer, parameter                           :: avail_prcl_mthd_nb = 1
  character(len=std_len), dimension(avail_prcl_mthd_nb), parameter :: &
                                                      parcel_method_list &
                                                      = ['fixed_numbers']
 
contains

end module mod_obj_global_parameters

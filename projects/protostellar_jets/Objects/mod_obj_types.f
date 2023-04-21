module mod_obj_types
  use mod_global_parameters

  implicit none

   !> User-defined parameters
   type usr_params
    !> Domain s parcelling method name. Possible choices are :
    !>  - 'fixed_numbers' (divide domain into fixed number of sub-regions)
    character(len=std_len) :: parcel_method
   
    !> If 'fixed_numbers', number of sub-regions
    integer                :: number_of_subregions

   end type usr_params

contains

end module mod_obj_types


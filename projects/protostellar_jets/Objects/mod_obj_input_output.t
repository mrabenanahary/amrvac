module mod_obj_input_output
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_types  
  use mod_obj_read_parameters
  use mod_obj_geometry

  implicit none

contains

  subroutine check_avail_strings(var,var_name,var_list,numb_avail)
   implicit none
   character(len=*), intent(in) :: var,var_name
   integer, intent(in)                :: numb_avail
   character(len=*), dimension(numb_avail), intent(in) :: var_list
   !----local----
   integer :: i
   logical :: avail_par
   !-------------------------------------

   avail_par = .false.

   do i=1,numb_avail
      if(trim(var)==trim(var_list(i)))then
         avail_par = .true.
      end if
   end do

   if(.not.avail_par)then
      write(*,*) 'You input ', trim(var_name), ' = ', trim(var)
      write(*,*) '> Possible inputs are : ', var_list
      call mpistop('> Please choose an available input !')
   end if
    
  end subroutine check_avail_strings  

  subroutine write_var_string(unit_config,var,var_name)
   implicit none
   integer,intent(in)           :: unit_config
   character(len=*), intent(in) :: var,var_name
   !-------------------------------------

   write(unit_config,*) trim(var_name),' = ', trim(var)

  end subroutine write_var_string  

  subroutine write_var_integer(unit_config,var,var_name)
   implicit none
   integer,intent(in)           :: unit_config
   integer, intent(in)          :: var
   character(len=*), intent(in) :: var_name
   !-------------------------------------

   write(unit_config,*) trim(var_name),' = ', var

  end subroutine write_var_integer 

  subroutine usr_write_setting(usr_config_var)
   implicit none
   type(usr_params), intent(inout) :: usr_config_var
   integer,parameter   :: unit_config = 12
   character(len=std_len)   :: filename_config
   !-------------------------------------
   filename_config=trim(base_filename)//'.config'

   open(unit_config,file=trim(filename_config), status='replace')
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%% Simulation configuration %%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

   write(unit_config,*)'*************************************************'
   write(unit_config,*)'************ User-defined parameters ************'
   write(unit_config,*)'*************************************************'
   write(unit_config,*)''
     
   call write_var_string(unit_config,parcel_method,&
         'Parcelling method name parcel_method ')

   call check_avail_strings(parcel_method,&
                            'parcel_method',&
                            parcel_method_list,&
                            avail_prcl_mthd_nb)
   
   select case(trim(trim(parcel_method)))
    case(trim(parcel_method_list(1))) !< 'fixed_numbers'
      call write_var_integer(unit_config,number_of_subregions,&
            'number of sub-regions number_of_subregions ')
   end select

   !< FOR TESTS ONLY to split parameters module with procedures module
   call usr_config_var%set_default_proc('set_value') !< check whether type s object-bound procedure isn t associated
                                                     !< in which case the corresponding parameter is set to a default value
   call write_var_string(unit_config,usr_config_var%is_it_null,&
         'is_it_null ')
   call change_procedure(usr_config_var) !< change the bound procedure from procedures module
   call usr_config_var%set_value() !< use type new object-bound procedure from procedures module
   call write_var_string(unit_config,usr_config_var%is_it_null,&
         'is_it_null ')

   
   write(unit_config,*)''
   write(unit_config,*)'*************************************************'
   write(unit_config,*)'******** End of user-defined parameters *********'
   write(unit_config,*)'*************************************************'

   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   write(unit_config,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   close(unit_config)

  end subroutine usr_write_setting

end module mod_obj_input_output


module mod_obj_input_output
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_types  

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
   type(usr_params), intent(in) :: usr_config_var
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
     
   call write_var_string(unit_config,usr_config_var%parcel_method,&
      'Parcelling method name')

   call check_avail_strings(usr_config_var%parcel_method,&
      'usr_config%parcel_method',parcel_method_list,avail_prcl_mthd_nb)
   
   select case(trim(trim(usr_config_var%parcel_method)))
    case(trim(parcel_method_list(1))) !< 'fixed_numbers'
      call write_var_integer(unit_config,usr_config_var%number_of_subregions,&
         'number of sub-regions')
   end select

   
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


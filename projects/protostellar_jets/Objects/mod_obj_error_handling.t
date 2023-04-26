module mod_obj_error_handling
  use mod_global_parameters
  use mod_obj_global_parameters
  use mod_obj_types  
 
  implicit none       

  type error_handler

  contains
    procedure, pass(self) :: check_string_in_list => search_string_in_list
    procedure, pass(self) :: check_integer_in_list => search_integer_in_list
  end type error_handler                                     
  
  contains

  !> \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ STRINGS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  function is_string_in_list(string,string_list,&
    n_elements) result(answer)
    implicit none
    integer                                  :: n_elements
    character(len=*), intent(in)             :: string
    character(len=*), dimension(n_elements)  :: string_list
    logical                                  :: answer
    !--local---------
    integer                                  :: i1
    !-------------------

    answer = .false.

    do i1=1,n_elements
      if(trim(string)==trim(string_list(i1))) answer = .true.
    end do

  end function is_string_in_list

  subroutine search_string_in_list(self,string,string_list,&
    n_elements)
    implicit none
    class(error_handler)                     :: self
    integer                                  :: n_elements
    character(len=*), intent(in)             :: string
    character(len=*), dimension(n_elements)  :: string_list
    !--------------------------------------------

    if(.not.is_string_in_list(string,string_list,n_elements))then
      write(*,*) string, ' is not in [', string_list, '] of size ', n_elements
      call mpistop('Unknown string passed to arguments') 
    end if

  end subroutine search_string_in_list

  !> \\\\\\\\\\\\\\\\\\\\\\\\\\\\ INTEGER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\

  function is_integer_in_list(i_input,integer_list,&
    n_elements) result(answer)
    implicit none
    integer                                  :: n_elements
    integer, intent(in)                      :: i_input
    integer, dimension(n_elements)           :: integer_list
    logical                                  :: answer
    !--local---------
    integer                                  :: i1
    !-------------------

    answer = .false.

    do i1=1,n_elements
      if(i_input==integer_list(i1)) answer = .true.
    end do

  end function is_integer_in_list

  subroutine search_integer_in_list(self,i_input,integer_list,&
    n_elements)
    implicit none
    class(error_handler)                     :: self
    integer                                  :: n_elements
    integer, intent(in)                      :: i_input
    integer, dimension(n_elements)           :: integer_list
    !--------------------------------------------

    if(.not.is_integer_in_list(i_input,integer_list,n_elements))then
      write(*,*) i_input, ' is not in [', integer_list, '] of size ', n_elements
      call mpistop('Unknown integer passed to arguments') 
    end if

  end subroutine search_integer_in_list

end module mod_obj_error_handling

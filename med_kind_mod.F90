module med_kind_mod

  use ESMF            , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError, ESMF_LOGMSG_ERROR, ESMF_MAXSTR
  use ESMF            , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE

  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------
  public
  integer,parameter :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte real
  integer,parameter :: SHR_KIND_RN = kind(1.0)              ! native real
  integer,parameter :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
  integer,parameter :: SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
  integer,parameter :: SHR_KIND_IN = kind(1)                ! native integer
  integer,parameter :: SHR_KIND_CS = 80                     ! short char
  integer,parameter :: SHR_KIND_CM = 160                    ! mid-sized char
  integer,parameter :: SHR_KIND_CL = 256                    ! long char
  integer,parameter :: SHR_KIND_CX = 512                    ! extra-long char
  integer,parameter :: SHR_KIND_CXX= 4096                   ! extra-extra-long char

  integer,  parameter :: med_constants_ispval_mask     = -987987 ! spval for RH mask values

  ! Module data
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

  contains

  logical function chkerr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module med_kind_mod

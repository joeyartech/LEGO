!public entities naming convention
! i : integer
! r : single precision real number (kind=4)
! d : double precision real number (kind=8)
! c : single precision complex number (kind=4)
! z : double precision complex number (kind=8)
! s : character, string
! m : module
! t : type, class
! o : optional
! is, if : logical

module m_global

    !default string length
    integer,parameter :: i_str_len  = 80
    integer,parameter :: i_str_xlen = 256
    integer,parameter :: i_str_xxlen = 512
    integer,parameter :: i_str_xxxlen = 1024
    character(*), parameter :: s_return = achar(13)
        !https://software.intel.com/en-us/forums/intel-fortran-compiler/topic/494946
        !char function: CHAR(10) is linefeed LF. CHAR(13) is carriage return CR. If you are a little paranoid, ACHAR(10) is better - this is a little more robust to the default character kind not being ascii.
        !The NEW_LINE standard intrinsic is even more robust. There's also the C_NEW_LINE and C_CARRIAGE_RETURN constants from the ISO_C_BINDING module for kind C_CHAR characters.

    !math constants
    real,parameter :: r_pi = 3.1415927
    double precision,parameter :: d_pi = 3.1415926535897932d0
    !real,parameter :: r_epsilon = epsilon(r_pi)
    
    !compilation info
    character(*),parameter :: s_commit = commit
    character(*),parameter :: s_branch = branch

#ifdef GNU
    character(*),parameter :: s_compiler = 'gfortran v' // __VERSION__
    character(*),parameter :: s_version = __VERSION__
    integer,parameter      :: i_endian = __BYTE_ORDER__
#endif

#ifdef INTEL
    character(*),parameter :: s_compiler = 'ifort'
    integer,parameter      :: s_version(2) = [__INTEL_COMPILER, __INTEL_COMPILER_UPDATE]
    character(*),parameter :: i_endian = 'unavailable'
#endif

end

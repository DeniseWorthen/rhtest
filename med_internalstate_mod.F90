module med_internalstate_mod

  use med_kind_mod , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8

  implicit none
  private

  ! Components
  integer, public :: compmed = 1
  integer, public :: compatm = 2
  integer, public :: complnd = 3
  integer, public :: compocn = 4
  integer, public :: compice = 5
  integer, public :: comprof = 6
  integer, public :: compwav = 7
  integer, public :: ncomps =  7 ! this will be incremented if the size of compglc is > 0
  integer, public, allocatable :: compglc(:)

  ! Generic component name (e.g. atm, ocn...)
  character(len=CS), public, allocatable :: compname(:)

  ! Mapping
  integer , public, parameter :: mapunset          = 0
  integer , public, parameter :: mapbilnr          = 1
  integer , public, parameter :: mapconsf          = 2
  integer , public, parameter :: mapconsd          = 3
  integer , public, parameter :: mappatch          = 4
  integer , public, parameter :: mapfcopy          = 5
  integer , public, parameter :: mapnstod          = 6  ! nearest source to destination
  integer , public, parameter :: mapnstod_consd    = 7  ! nearest source to destination followed by conservative dst
  integer , public, parameter :: mapnstod_consf    = 8  ! nearest source to destination followed by conservative frac
  integer , public, parameter :: mappatch_uv3d     = 9  ! rotate u,v to 3d cartesian space, map from src->dest, then rotate back
  integer , public, parameter :: mapbilnr_uv3d     = 10 ! rotate u,v to 3d cartesian space, map from src->dest, then rotate back
  integer , public, parameter :: map_rof2ocn_ice   = 11 ! custom smoothing map to map ice from rof->ocn (cesm only)
  integer , public, parameter :: map_rof2ocn_liq   = 12 ! custom smoothing map to map liq from rof->ocn (cesm only)
  integer , public, parameter :: map_glc2ocn_liq   = 13 ! custom smoothing map to map liq from glc->ocn (cesm only)
  integer , public, parameter :: map_glc2ocn_ice   = 14 ! custom smoothing map to map ice from glc->ocn (cesm only)
  integer , public, parameter :: mapfillv_bilnr    = 15 ! fill value followed by bilinear
  integer , public, parameter :: mapbilnr_nstod    = 16 ! bilinear with nstod extrapolation
  integer , public, parameter :: mapconsf_aofrac   = 17 ! conservative with aofrac normalization (ufs only)
  integer , public, parameter :: nmappers          = 17

  character(len=*) , public, parameter :: mapnames(nmappers) = &
       (/'bilnr       ',&
         'consf       ',&
         'consd       ',&
         'patch       ',&
         'fcopy       ',&
         'nstod       ',&
         'nstod_consd ',&
         'nstod_consf ',&
         'patch_uv3d  ',&
         'bilnr_uv3d  ',&
         'rof2ocn_ice ',&
         'rof2ocn_liq ',&
         'glc2ocn_ice ',&
         'glc2ocn_liq ',&
         'fillv_bilnr ',&
         'bilnr_nstod ',&
         'consf_aofrac'/)

  ! Coupling mode
  character(len=CS), public :: coupling_mode ! valid values are [cesm,nems_orig,nems_frac,nems_orig_data,hafs]

  logical, public :: dststatus_print = .true.
end module med_internalstate_mod

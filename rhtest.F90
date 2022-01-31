  program rhtest

  use ESMF
  use aMesh
  use med_internalstate_mod
  use med_map_mod

  implicit none

  type(ESMF_VM)           :: vm
  type(ESMF_DistGrid)     :: distGrid
  type(ESMF_Mesh)         :: atmMesh
  type(ESMF_Mesh)         :: ocnMesh
  type(ESMF_Field)        :: aField, oField
  type(ESMF_End_Flag)     :: endflag
  type(ESMF_LogKind_Flag) :: logkindflag
  type(ESMF_RouteHandle), pointer  :: RH(:,:,:)

  integer :: petCount, localPet
  integer :: rc
  character(len=5) :: mapfile = 'unset'

  call ESMF_Initialize(vm=vm, logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  allocate(RH(ncomps,ncomps,nmappers))

  allocate(compname(ncomps))
  compname(compmed) = 'med'
  compname(compatm) = 'atm'
  compname(complnd) = 'lnd'
  compname(compocn) = 'ocn'
  compname(compice) = 'ice'
  compname(comprof) = 'rof'
  compname(compwav) = 'wav'

  coupling_mode = 'nems_frac'

  call setupAtmMesh(atmMesh, rc=rc)
 
  ocnMesh = ESMF_MeshCreate(filename='data/mesh.mx100.nc', fileformat=ESMF_FILEFORMAT_ESMFMESH,rc=rc)

  aField = ESMF_FieldCreate(atmMesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
  oField = ESMF_FieldCreate(ocnMesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)

  call med_map_routehandles_initfrom_field(compocn, compatm, oField, aField, mapconsf, RH(compocn,compatm,:), mapfile, rc)

  call med_map_routehandles_initfrom_field(compocn, compatm, oField, aField, mapconsd, RH(compocn,compatm,:), mapfile, rc)

  call med_map_routehandles_initfrom_field(compatm, compocn, aField, oField, mapbilnr_nstod, RH(compatm,compocn,:), mapfile, rc)

  end program rhtest

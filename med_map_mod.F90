module med_map_mod

  use med_kind_mod          , only : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_kind_mod          , only : I4=>SHR_KIND_I4
  use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
  use ESMF                  , only : ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF                  , only : ESMF_Field
  use med_kind_mod          , only : chkerr

  implicit none
  private

  ! public routines
  public :: med_map_routehandles_initfrom_field
 
  character(*),parameter :: u_FILE_u = &
       __FILE__
  contains

  !================================================================================
  subroutine med_map_routehandles_initfrom_field(n1, n2, fldsrc, flddst, mapindex, routehandles, mapfile, rc)

    use ESMF                  , only : ESMF_RouteHandle, ESMF_RouteHandlePrint, ESMF_Field, ESMF_MAXSTR
    use ESMF                  , only : ESMF_PoleMethod_Flag, ESMF_POLEMETHOD_ALLAVG, ESMF_POLEMETHOD_NONE
    use ESMF                  , only : ESMF_FieldSMMStore, ESMF_FieldRedistStore, ESMF_FieldRegridStore
    use ESMF                  , only : ESMF_RouteHandleIsCreated, ESMF_RouteHandleCreate
    use ESMF                  , only : ESMF_REGRIDMETHOD_BILINEAR, ESMF_REGRIDMETHOD_PATCH
    use ESMF                  , only : ESMF_REGRIDMETHOD_CONSERVE, ESMF_NORMTYPE_DSTAREA, ESMF_NORMTYPE_FRACAREA
    use ESMF                  , only : ESMF_UNMAPPEDACTION_IGNORE, ESMF_REGRIDMETHOD_NEAREST_STOD
    use ESMF                  , only : ESMF_EXTRAPMETHOD_NEAREST_STOD
    use ESMF                  , only : ESMF_Mesh, ESMF_MeshLoc, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_I4
    use ESMF                  , only : ESMF_MeshGet, ESMF_DistGridGet, ESMF_DistGrid, ESMF_TYPEKIND_R8
    use ESMF                  , only : ESMF_FieldGet, ESMF_FieldCreate, ESMF_FieldWrite, ESMF_FieldDestroy
    use ESMF                  , only : ESMF_Array, ESMF_ArrayCreate, ESMF_ArrayWrite, ESMF_FieldFill
    use med_internalstate_mod , only : mapbilnr, mapconsf, mapconsd, mappatch, mappatch_uv3d, mapbilnr_uv3d, mapfcopy
    use med_internalstate_mod , only : mapunset, mapnames, nmappers
    use med_internalstate_mod , only : mapnstod, mapnstod_consd, mapnstod_consf, mapnstod_consd
    use med_internalstate_mod , only : mapfillv_bilnr, mapbilnr_nstod, mapconsf_aofrac
    use med_internalstate_mod , only : ncomps, compatm, compice, compocn, compwav, compname
    use med_internalstate_mod , only : coupling_mode, dststatus_print
    use med_kind_mod          , only : ispval_mask => med_constants_ispval_mask

    ! input/output variables
    integer                    , intent(in)    :: n1
    integer                    , intent(in)    :: n2
    type(ESMF_Field)           , intent(inout) :: fldsrc
    type(ESMF_Field)           , intent(inout) :: flddst
    integer                    , intent(in)    :: mapindex
    type(ESMF_RouteHandle)     , intent(inout) :: routehandles(:)
    character(len=*), optional , intent(in)    :: mapfile
    integer                    , intent(out)   :: rc

    ! local variables
    type(ESMF_Mesh)            :: dstmesh
    type(ESMF_Field)           :: dststatusfield, doffield
    type(ESMF_DistGrid)        :: distgrid
    type(ESMF_Array)           :: elemMaskArray
    character(len=CS)          :: string
    character(len=CS)          :: mapname
    character(len=CL)          :: fname
    integer                    :: srcMaskValue
    integer                    :: dstMaskValue
    character(len=ESMF_MAXSTR) :: lmapfile
    logical                    :: rhprint = .false., ldstprint = .false.
    integer                    :: ns
    integer(I4), pointer       :: dof(:)
    integer(I4), pointer       :: meshmask(:)
    integer                    :: srcTermProcessing_Value = 0
    type(ESMF_PoleMethod_Flag) :: polemethod
    integer             :: logunit = 6
    integer             :: mastertask = 0
    character(len=*), parameter :: subname=' (module_med_map: med_map_routehandles_initfrom_field) '
    !---------------------------------------------

    lmapfile = 'unset'
    if (present(mapfile)) then
       lmapfile = trim(mapfile)
    end if

    mapname = trim(mapnames(mapindex))
    call ESMF_LogWrite(trim(subname)//": mapname "//trim(mapname), ESMF_LOGMSG_INFO)

    ! create a field to retrieve the dststatus field
    call ESMF_FieldGet(flddst, mesh=dstmesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    dststatusfield = ESMF_FieldCreate(dstmesh, ESMF_TYPEKIND_I4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldFill(dststatusfield, dataFillScheme="const", const1=-987987.0_R8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! set local flag to false
    ldstprint = .false.

    polemethod=ESMF_POLEMETHOD_ALLAVG
    if (trim(coupling_mode) == 'cesm') then
       dstMaskValue = ispval_mask
       srcMaskValue = ispval_mask
       if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
       if (n2 == compocn .or. n2 == compice) dstMaskValue = 0
       if (n1 == compwav .and. n2 == compocn) then
         srcMaskValue = 0
         dstMaskValue = ispval_mask
      endif
      if (n1 == compwav .or. n2 == compwav) then
        polemethod = ESMF_POLEMETHOD_NONE ! todo: remove this when ESMF tripolar mapping fix is in place.
      endif
    else if (coupling_mode(1:4) == 'nems') then
       if ( (n1 == compocn .or. n1 == compice .or. n1 == compwav) .and. &
            (n2 == compocn .or. n2 == compice .or. n2 == compwav) ) then
          srcMaskValue = 0
          dstMaskValue = 0
       else if (n1 == compatm .and. (n2 == compocn .or. n2 == compice .or. n2 == compwav)) then
          srcMaskValue = 1
          dstMaskValue = 0
       else if (n2 == compatm .and. (n1 == compocn .or. n1 == compice .or. n1 == compwav)) then
          srcMaskValue = 0
          dstMaskValue = 1
          !dstMaskValue = ispval_mask
       else
          ! TODO: what should the condition be here?
          dstMaskValue = ispval_mask
          srcMaskValue = ispval_mask
       end if
    else if (trim(coupling_mode) == 'hafs') then
       dstMaskValue = ispval_mask
       srcMaskValue = ispval_mask
       if (n1 == compocn .or. n1 == compice) srcMaskValue = 0
       if (n2 == compocn .or. n2 == compice) dstMaskValue = 0
       if (n1 == compatm .and. n2 == compocn) then
          dstMaskValue = 0
       elseif (n1 == compocn .and. n2 == compatm) then
          srcMaskValue = 0
          dstMaskValue = ispval_mask
       elseif (n1 == compatm .and. n2 == compwav) then
          dstMaskValue = 0
       elseif (n1 == compwav .and. n2 == compatm) then
          srcMaskValue = 0
          dstMaskValue = ispval_mask
       endif
    end if

    write(string,'(a,i8,a,i8)') trim(compname(n1))//' to '//trim(compname(n2))//' srcMask = ', &
               srcMaskValue,' dstMask = ',dstMaskValue
    call ESMF_LogWrite(trim(string), ESMF_LOGMSG_INFO)

    ! Create route handle
    if (mapindex == mapfcopy) then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH redist for '//trim(string)
       end if
       call ESMF_FieldRedistStore(fldsrc, flddst, routehandle=routehandles(mapfcopy), &
            ignoreUnmatchedIndices = .true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (lmapfile /= 'unset') then
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//&
               ' via input file '//trim(mapfile)//' for '//trim(string)
       end if
       call ESMF_FieldSMMStore(fldsrc, flddst, mapfile, routehandle=routehandles(mapindex), &
            ignoreUnmatchedIndices=.true., &
            srcTermProcessing=srcTermProcessing_Value, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else if (mapindex == mapbilnr .or. mapindex == mapbilnr_uv3d) then
             write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
          call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapbilnr), &
               srcMaskValues=(/srcMaskValue/), &
               dstMaskValues=(/dstMaskValue/), &
               regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
               polemethod=polemethod, &
               srcTermProcessing=srcTermProcessing_Value, &
               ignoreDegenerate=.true., &
               dstStatusField=dststatusfield, &
               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ldstprint = .true.
    else if (mapindex == mapfillv_bilnr) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapfillv_bilnr), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
            polemethod=polemethod, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mapbilnr_nstod) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapbilnr_nstod), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
            polemethod=polemethod, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mapconsf .or. mapindex == mapnstod_consf) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapconsf), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_FRACAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mapconsf_aofrac) then
             write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
          call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapconsf_aofrac), &
               srcMaskValues=(/srcMaskValue/), &
               dstMaskValues=(/dstMaskValue/), &
               regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
               normType=ESMF_NORMTYPE_FRACAREA, &
               srcTermProcessing=srcTermProcessing_Value, &
               ignoreDegenerate=.true., &
               dstStatusField=dststatusfield, &
               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
               rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          ldstprint = .true.
    else if (mapindex == mapconsd .or. mapindex == mapnstod_consd) then
          write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapconsd), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            normType=ESMF_NORMTYPE_DSTAREA, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.
    else if (mapindex == mappatch .or. mapindex == mappatch_uv3d) then
             write(logunit,'(A)') trim(subname)//' creating RH '//trim(mapname)//' for '//trim(string)
          call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mappatch), &
               srcMaskValues=(/srcMaskValue/), &
               dstMaskValues=(/dstMaskValue/), &
               regridmethod=ESMF_REGRIDMETHOD_PATCH, &
               polemethod=polemethod, &
               srcTermProcessing=srcTermProcessing_Value, &
               ignoreDegenerate=.true., &
               dstStatusField=dststatusfield, &
               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       if (mastertask) then
          write(logunit,'(A)') trim(subname)//' mapindex '//trim(mapname)//' not supported for '//trim(string)
       end if
       call ESMF_LogWrite(trim(subname)//' mapindex '//trim(mapname)//' not supported ', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
       rc = ESMF_FAILURE
       return
    end if

    ! Output destination status field to file if requested
    if (dststatus_print .and. ldstprint) then
      fname = 'dststatus.'//trim(compname(n1))//'.'//trim(compname(n2))//'.'//trim(mapname)//'.nc'
      call ESMF_LogWrite(trim(subname)//": writing dstStatusField to "//trim(fname), ESMF_LOGMSG_INFO)

      call ESMF_FieldWrite(dststatusfield, filename=trim(fname), variableName='dststatus', &
           overwrite=.true., rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! the sequence index in order to sort the dststatus field
      call ESMF_MeshGet(dstmesh, elementDistgrid=distgrid, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_DistGridGet(distgrid, localDE=0, elementCount=ns, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(dof(ns))
      call ESMF_DistGridGet(distgrid, localDE=0, seqIndexList=dof, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      doffield = ESMF_FieldCreate(dstmesh, dof, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldWrite(doffield, fileName='dof.'//trim(compname(n2))//'.nc', variableName='dof', &
           overwrite=.true., rc=rc)

      allocate(meshmask(ns))
      meshmask(:) = 0
      ! get the meshmask
      elemMaskArray = ESMF_ArrayCreate(Distgrid, farrayPtr=meshmask, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_MeshGet(dstmesh, elemMaskArray=elemMaskArray, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_ArrayWrite(elemMaskArray, 'dststatus.meshmask.'//trim(compname(n2))//'.nc', variableName = 'mask', &
           overwrite=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      deallocate(dof)
      deallocate(meshmask)
      call ESMF_FieldDestroy(doffield, rc=rc, noGarbage=.true.)
    end if

    ! consd_nstod method requires a second routehandle
    if (mapindex == mapnstod .or. mapindex == mapnstod_consd .or. mapindex == mapnstod_consf) then
       call ESMF_FieldRegridStore(fldsrc, flddst, routehandle=routehandles(mapnstod), &
            srcMaskValues=(/srcMaskValue/), &
            dstMaskValues=(/dstMaskValue/), &
            regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., &
            dstStatusField=dststatusfield, &
            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
            rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ldstprint = .true.

       ! Output destination status field to file if requested
       if (dststatus_print .and. ldstprint) then
          fname = 'dststatus.'//trim(compname(n1))//'.'//trim(compname(n2))//'.'//trim(mapname)//'_2.nc'
          call ESMF_LogWrite(trim(subname)//": writing dstStatusField to "//trim(fname), ESMF_LOGMSG_INFO)

          call ESMF_FieldWrite(dststatusfield, filename=trim(fname), variableName='dststatus', overwrite=.true., rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Output route handle to file if requested
    if (rhprint) then
       if (mastertask) then
          write(logunit,'(a)') trim(subname)//trim(string)//": printing  RH for "//trim(mapname)
       end if
       call ESMF_RouteHandlePrint(routehandles(mapindex), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_FieldDestroy(dststatusfield, rc=rc, noGarbage=.true.)

  end subroutine med_map_routehandles_initfrom_field
end module med_map_mod

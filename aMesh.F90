module aMesh

   use ESMF
   use med_kind_mod          , only : chkerr

   implicit none


   character(*),parameter :: u_FILE_u = &
        __FILE__

   contains

   subroutine setupAtmMesh(atmMesh, rc)

     type(ESMF_Mesh) :: atmMesh
     integer, intent(out) :: rc

     ! local variables
     type(ESMF_Grid)         :: grid2D
     type(ESMF_DistGrid)     :: distGrid
     type(ESMF_Mesh)         :: mesh0
     type(ESMF_Array)        :: elemMaskArray
     type(ESMF_Field)        :: doffield, maskfield
     type(ESMF_Decomp_Flag)  :: decompflagPTile(2,6)

     integer :: ncnt, tl
     integer, pointer :: dof(:)
     integer, pointer :: meshmask(:)
     integer, pointer :: lmask(:)

     integer, dimension(2,6) :: decomptile

     rc = ESMF_SUCCESS

     ! Set up decomposition for each tile
     !decomptile(:,1)=(/2,2/) ! Tile 1
     !decomptile(:,2)=(/2,2/) ! Tile 2
     !decomptile(:,3)=(/1,2/) ! Tile 3
     !decomptile(:,4)=(/1,2/) ! Tile 4
     !decomptile(:,5)=(/1,2/) ! Tile 5
     !decomptile(:,6)=(/1,2/) ! Tile 6
     do tl = 1,6
              decomptile(1,tl) = 3
              decomptile(2,tl) = 8
              decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
     end do
     ! Create cubed sphere grid and read in the center and corner stagger coordinates
     ! from the tile files
     grid2D = ESMF_GridCreateMosaic(filename='data/C96_mosaic.nc', &
                decompflagPTile=decompflagPTile,                                   &
                staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
                tileFilePath='./data/', regDecompPTile=decomptile, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return

     call ESMF_GridAddItem(grid2D, staggerLoc=ESMF_STAGGERLOC_CENTER, itemflag=ESMF_GRIDITEM_MASK, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return

     ! create an fv3 mesh from the grid
     mesh0 = ESMF_MeshCreate(grid=grid2D, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     call ESMF_MeshGet(mesh0, elementCount=ncnt, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return

     ! create a field on mesh0 and use it to retrieve the distgrid cmeps run
     doffield = ESMF_FieldCreate(mesh=mesh0, typekind=ESMF_TYPEKIND_I4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     call ESMF_FieldGet(field=doffield, farrayPtr=dof, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     call ESMF_FieldRead(field=doffield, fileName='data/dof.atm.nc', variableName='dof', rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     print *,dof(1:10)

     ! create a distgrid to match the dof
     distGrid = ESMF_DistGridCreate(arbSeqIndexList=dof, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     ! recreate the mesh using the above distGrid (this may actually be the same distgrid created by cmeps?)
     atmMesh = ESMF_MeshCreate(mesh=mesh0, elementDistgrid=distGrid, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return

    allocate(meshmask(ncnt))
    elemMaskArray = ESMF_ArrayCreate(Distgrid, farrayPtr=meshmask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(atmMesh, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

     maskfield = ESMF_FieldCreate(mesh=atmMesh, typekind=ESMF_TYPEKIND_I4, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     call ESMF_FieldGet(field=maskfield, farrayPtr=lmask, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return
     ! now get the land mask from the file
     call ESMF_FieldRead(field=maskfield, fileName='data/dststatus.meshmask.atm.nc', variableName='mask', rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return

     meshmask(1:ncnt) = lmask(1:ncnt)
     call ESMF_MeshSet(mesh=atmMesh, elementMask=meshmask, rc=rc)
     if (chkerr(rc,__LINE__,u_FILE_u)) return

      call ESMF_ArrayWrite(elemMaskArray, 'test.nc', variableName = 'mask', &
           overwrite=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
   end subroutine setupAtmMesh
end module aMesh

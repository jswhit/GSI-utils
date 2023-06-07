program taperens
!$$$  main program documentation block
!
! program:  taperens               recenter
!
! prgmmr: whitaker         org: esrl/psd               date: 2009-02-23
!
! abstract:  Read ensemble from netcdf history files,
!            compute ens mean, optionally recenter around specified
!            ens mean, taper ens perturbations at top of model,
!            and write out results.
!
! program history log:
!   2023-6-1  Initial version.
!
! usage:
!   input files:
!
!   output files:
!
! attributes:
!   language: f95
!
!
!$$$

  use module_ncio, only: open_dataset, create_dataset, read_attribute, &
                         Dataset, Dimension, close_dataset, has_attr, has_var, &
                         read_vardata, write_attribute, write_vardata, &
                         get_dim, quantize_data

  implicit none

  include "mpif.h"

  logical:: quantize,recenterens,lexist,hybgain,clip,cliptracers,tracer
  character(len=500) filenamein,filenameout,filename_meanin, &
    filename_fgmean,filename_3dvar
  character(len=3) charnanal
  integer iret,mype,mype1,npe,nanals,ierr,k
  integer:: latb,lonb,levs,nbits,nvar,ndims
  real(4),allocatable, dimension(:) :: ak,bk
  real(4),allocatable, dimension(:,:,:) :: values_3d, values_3d_i,&
                           values_3d_meani,values_3d_mean, taper_vert, &
                           values_3d_fgmean, values_3d_3dvar
  real(4),allocatable, dimension(:,:) :: values_2d, values_2d_i,&
         values_2d_fgmean,values_2d_3dvar,values_2d_mean,values_2d_meani
  real(4) :: compress_err,ak_top,ak_bot,alpha,beta
  real(8) :: rnanals

  type(Dataset) :: dseti,dseto,dsetmi,dsetmo,dsetfgm,dset3dvar,dsetmo_new
  type(Dimension) :: londim,latdim,levdim

  namelist /setup/ ak_top, ak_bot, alpha, beta, filenamein, &
                   filenameout, filename_meanin, nanals, &
                   filename_fgmean, filename_3dvar, cliptracers

! Initialize mpi
  call MPI_Init(ierr)

! mype is process number, npe is total number of processes.
  call MPI_Comm_rank(MPI_COMM_WORLD,mype,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,npe,ierr)

  if (mype==0) call w3tagb('TAPERENS',2011,0319,0055,'NP25')

! defaults for namelist
  alpha=0 ! no hybrid gain
  beta=1
  cliptracers = .false.
  filename_fgmean = ""
  filename_3dvar = ""
  filenamein = ""
  filenameout = ""
  filename_meanin = ""
  nanals = 0
  ak_bot = -1
  ak_top = -1

! read namelist on each rank
  inquire(file='taperens.nml', exist=lexist)
  if ( lexist ) then
    open(file='taperens.nml', unit=92, status='old', &
         form='formatted', action='read', access='sequential')
    read(92,nml=setup)
    close(92)
  else
    if (mype .eq. 0) write(6,*) 'taperens.nml does not exist and should, ABORT!'
    call mpi_abort(mpi_comm_world,99,iret)
    stop
  endif

  if (mype==0) then
     write(6,*)'TAPERENS:  PROCESS ',nanals,' ENSEMBLE MEMBERS'
     write(6,*)'filenamein=',trim(filenamein)
     if (trim(filename_meanin) .ne. "") then
        write(6,*)'recenter around ',trim(filename_meanin)
     else
        write(6,*)'no recentering'
     endif
     if (alpha > 0) then
        write(6,*) 'hybrid gain'
        write(6,*)'mean background forecast from ',trim(filename_fgmean)
        write(6,*)'3dvar analysis from ',trim(filename_3dvar)
     endif
     write(6,*)'filenameout=',trim(filenameout)
     write(6,*)'ak_bot,ak_top',ak_bot,ak_top
     write(6,*)'alpha,beta',alpha,beta
  endif

  if (npe .ne. nanals) then
     if (mype .eq. 0) write(6,'(2(a,i4))')'***FATAL ERROR***  npe not equal to nanals.  npe = ',npe,' != nanals = ',nanals
     call mpi_abort(mpi_comm_world,99,iret)
     stop
  endif

  if (trim(filename_meanin) .ne. "") then
     if (alpha > 0) then
        if (mype .eq. 0) write(6,*)'***FATAL ERROR***  both alpha and filename_meanin specified'
        call mpi_abort(mpi_comm_world,99,iret)
        stop
     endif
     dsetmi  = open_dataset(trim(filename_meanin))
     recenterens = .true.
  else
     recenterens = .false.
  endif
  if (alpha > 1.e-5) then
     hybgain = .true.
     dsetfgm  = open_dataset(trim(filename_fgmean))
     dset3dvar  = open_dataset(trim(filename_3dvar))
  else
     hybgain = .false.
  endif
  if (mype .eq. 0) print *,'recenterens,hybgain',recenterens,hybgain

  clip = tiny(alpha)
  rnanals = 1./nanals

  mype1 = mype+1
  write(charnanal,'(i3.3)') mype1
  dseti  = open_dataset(trim(filenamein)//"_mem"//charnanal)
  londim = get_dim(dsetmi,'grid_xt'); lonb = londim%len
  latdim = get_dim(dsetmi,'grid_yt'); latb = latdim%len
  levdim = get_dim(dsetmi,'pfull');   levs = levdim%len
  dseto  = create_dataset(trim(filenameout)//"_mem"//charnanal, dseti, copy_vardata=.true.)
  if (mype .eq. 0) then
     if (hybgain) then
         dsetmo_new  = create_dataset(trim(filenameout)//"_ensmean", dseti, copy_vardata=.true.)
         dsetmo  = create_dataset(trim(filenameout)//"_ensmean.orig", dseti, copy_vardata=.true.)
     else
         dsetmo  = create_dataset(trim(filenameout)//"_ensmean", dseti, copy_vardata=.true.)
     endif
  endif

  allocate(taper_vert(lonb,latb,levs))
  taper_vert = 1.0
  allocate(values_3d_mean(lonb,latb,levs))
  allocate(values_3d_fgmean(lonb,latb,levs))
  allocate(values_3d_3dvar(lonb,latb,levs))
  allocate(values_3d_meani(lonb,latb,levs))
  allocate(values_3d(lonb,latb,levs))
  allocate(values_3d_i(lonb,latb,levs))
  allocate(values_2d_mean(lonb,latb))
  allocate(values_2d_fgmean(lonb,latb))
  allocate(values_2d_3dvar(lonb,latb))
  allocate(values_2d_meani(lonb,latb))
  allocate(values_2d(lonb,latb))
  allocate(values_2d_i(lonb,latb))
  allocate(ak(levs+1),bk(levs+1))

  call read_attribute(dseti, 'ak', ak)
  call read_attribute(dseti, 'bk', bk)

  ! vertical taper function for ens perts
  if (ak_bot > 0) then
     if (mype .eq. 0) print *,'profile to taper ens perts (k,ak,bk,taper):'
     do k=1,levs
        if (k < levs/2 .and. (ak(k) <= ak_bot .and. ak(k) >= ak_top)) then
           taper_vert(:,:,k)= log(ak(k) - ak_top)/log(ak_bot - ak_top)
        else if (bk(k) .eq. 0. .and. ak(k) < ak_top) then
           taper_vert(:,:,k) = 0.
        endif
        if (mype .eq. 0) then
          print *,k,ak(k),bk(k),taper_vert(1,1,k)
        endif
     enddo
  endif

  ! loop over all variables
  do nvar=1,dseti%nvars
     ndims = dseti%variables(nvar)%ndims
     if (trim(dseti%variables(nvar)%name) == 'spfh' .or. &
         trim(dseti%variables(nvar)%name) == 'o3mr' .or. &
         trim(dseti%variables(nvar)%name) == 'clwmr' .or. &
         trim(dseti%variables(nvar)%name) == 'icmr') then
           tracer = .true.
     else
           tracer = .false.
     endif
     if (trim(dseti%variables(nvar)%name) == 'pressfc') then
        ! pressfc is the only 2d var that varies across ensemble
        if (mype == 0) print *,'processing ',trim(dseti%variables(nvar)%name)
        call read_vardata(dseti,trim(dseti%variables(nvar)%name),values_2d_i)
        call mpi_allreduce(values_2d_i,values_2d_mean,lonb*latb,mpi_real4,mpi_sum,mpi_comm_world,iret)
        values_2d_mean = values_2d_mean*rnanals
        print *,mype,'ps mean',minval(values_2d_mean),maxval(values_2d_mean)
        if (hybgain) then ! hybrid gain (blend 3dvar and EnKF mean incr)
           call read_vardata(dsetfgm,trim(dseti%variables(nvar)%name),values_2d_fgmean)
           call read_vardata(dset3dvar,trim(dseti%variables(nvar)%name),values_2d_3dvar)
           values_2d_meani = (1. - alpha - beta)*values_2d_fgmean + &
                             alpha*values_2d_3dvar + beta*values_2d_mean
        else if (recenterens) then ! recenter around specified mean
           call read_vardata(dsetmi,trim(dseti%variables(nvar)%name),values_2d_meani)
        else
           values_2d_meani = values_2d_mean ! don't change mean
        endif
        values_2d = values_2d_i - values_2d_mean + values_2d_meani
        if (has_attr(dseti, 'nbits', trim(dseti%variables(nvar)%name))) then
           call read_attribute(dseti, 'nbits', nbits, &
                trim(dseti%variables(nvar)%name))
           quantize = .true.
           if (nbits < 1) quantize = .false.
        else
           quantize = .false.
        endif
        if (recenterens) then
           if (quantize) then
              values_2d_i = values_2d
              call quantize_data(values_2d_i, values_2d, nbits, compress_err)
              call write_attribute(dseto,&
              'max_abs_compression_error',compress_err,trim(dseti%variables(nvar)%name))
           endif
           call write_vardata(dseto,trim(dseti%variables(nvar)%name),values_2d)
        end if
        ! write ens mean
        if (mype == 0) then
           if (quantize) then
             values_2d_i = values_2d_mean
             call quantize_data(values_2d_i, values_2d, nbits, compress_err)
             call write_attribute(dseto,&
             'max_abs_compression_error',compress_err,trim(dseti%variables(nvar)%name))
           else
             values_2d = values_2d_mean
           endif
           call write_vardata(dsetmo,trim(dseti%variables(nvar)%name),values_2d)
           ! for hyb gain, write new mean file (blended increments)
           if (hybgain) then
              if (quantize) then
                values_2d_i = values_2d_meani
                call quantize_data(values_2d_i, values_2d, nbits, compress_err)
                call write_attribute(dseto,&
                'max_abs_compression_error',compress_err,trim(dseti%variables(nvar)%name))
              else
                values_2d = values_2d_meani
              endif
              call write_vardata(dsetmo_new,trim(dseti%variables(nvar)%name),values_2d)
           endif
        endif
     else if (ndims == 4) then ! 3d variables
        if (mype == 0) print *,'processing ',trim(dseti%variables(nvar)%name)
        call read_vardata(dseti,trim(dseti%variables(nvar)%name),values_3d_i)
        call mpi_allreduce(values_3d_i,values_3d_mean,lonb*latb*levs,mpi_real4,mpi_sum,mpi_comm_world,iret)
        values_3d_mean = values_3d_mean*rnanals
        if (mype .eq. 0) then
          print *,trim(dseti%variables(nvar)%name),' mean',minval(values_3d_mean),maxval(values_3d_mean)
        endif
        if (hybgain) then ! hybrid gain (blend 3dvar and EnKF mean incr)
           call read_vardata(dsetfgm,trim(dseti%variables(nvar)%name),values_3d_fgmean)
           call read_vardata(dset3dvar,trim(dseti%variables(nvar)%name),values_3d_3dvar)
           values_3d_meani = (1. - alpha - beta)*values_3d_fgmean + &
                             alpha*values_3d_3dvar + beta*values_3d_mean
        else if (recenterens) then ! recenter around specified mean
           call read_vardata(dsetmi,trim(dseti%variables(nvar)%name),values_3d_meani)
        else
           values_3d_meani = values_3d_mean ! don't change mean
        endif
        values_3d = (values_3d_i - values_3d_mean)*taper_vert + values_3d_meani
        if (has_attr(dseti, 'nbits', trim(dseti%variables(nvar)%name))) then
           call read_attribute(dseti, 'nbits', nbits, &
                trim(dseti%variables(nvar)%name))
           quantize = .true.
           if (nbits < 1) quantize = .false.
        else
           quantize = .false.
        endif
        if (quantize) then
           values_3d_i = values_3d
           call quantize_data(values_3d_i, values_3d, nbits, compress_err)
           call write_attribute(dseto,&
           'max_abs_compression_error',compress_err,trim(dseti%variables(nvar)%name))
        endif
        if (tracer) then
          if (cliptracers)  where (values_3d < clip) values_3d = clip
          if (mype == 0) &
          print *,'clipping ',trim(dseti%variables(nvar)%name)
        endif
        call write_vardata(dseto,trim(dseti%variables(nvar)%name),values_3d)
        ! write out ensmean
        if (mype .eq. 0) then
           if (quantize) then
             values_3d_i = values_3d_mean
             call quantize_data(values_3d_i, values_3d, nbits, compress_err)
             call write_attribute(dsetmo,&
             'max_abs_compression_error',compress_err,trim(dseti%variables(nvar)%name))
           else
             values_3d = values_3d_mean
           endif
           call write_vardata(dsetmo,trim(dseti%variables(nvar)%name),values_3d)
           ! for hyb gain, write new mean file (blended increments)
           if (hybgain) then
              if (quantize) then
                values_3d_i = values_3d_meani
                call quantize_data(values_3d_i, values_3d, nbits, compress_err)
                call write_attribute(dseto,&
                'max_abs_compression_error',compress_err,trim(dseti%variables(nvar)%name))
              else
                values_3d = values_3d_meani
              endif
              call write_vardata(dsetmo_new,trim(dseti%variables(nvar)%name),values_3d)
           endif
        endif
     endif 
  enddo  ! nvars

  deallocate(values_3d,values_3d_i,values_3d_mean,values_3d_meani,taper_vert)
  deallocate(values_2d,values_2d_i,values_2d_mean,values_2d_meani)
  deallocate(values_3d_fgmean,values_2d_fgmean,values_3d_3dvar,values_2d_3dvar)
  call close_dataset(dseti)
  call close_dataset(dseto)
  if (recenterens) call close_dataset(dsetmi)
  if (hybgain) then
    call close_dataset(dsetfgm)
    call close_dataset(dset3dvar)
  endif
  if (mype .eq. 0) then
      call close_dataset(dsetmo)
      if (hybgain) call close_dataset(dsetmo_new)
  endif

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  if (mype==0) call w3tage('TAPERENS')

  call MPI_Finalize(ierr)
  if (mype .eq. 0 .and. ierr .ne. 0) then
     print *, 'MPI_Finalize error status = ',ierr
  end if

END program taperens

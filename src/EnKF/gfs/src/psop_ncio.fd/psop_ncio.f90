program psop

 ! EnKF surface pressure forward operator.

 use module_ncio, only: open_dataset, create_dataset, read_attribute, &
                        Dataset, Dimension, close_dataset, has_attr, has_var, &
                        read_vardata, write_attribute, write_vardata, get_dim
 use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
                              nc_diag_write, nc_diag_data2d
 implicit none
 include  'mpif.h'
 type(Dataset) :: dset
 type(Dimension) :: londim,latdim,levdim
 character(len=120) filenamein,obsfile,filename,diag_file
 character(len=10) datestring
 integer iret,nlats,nlons,nlevs,ntrac,ntrunc,ierr,nanals,nfhr,nobstot,&
         nproc,numproc,nob,nanal,j,iunit,iunitsig,fhmin,fhmax,fhout,ntimes,&
         nchar,nreal,ii,nn,nlevt,ntime,np,k,nobsh,izob,iunit_nml,iunito,idate
 integer mpi_status(mpi_status_size),ianldate
 real dxob,dyob,dtob,zerr,anal_obt,anal_obp,rlapse,&
      delz_const,ensmean_ob,bias,preduce,palt,zthresh,zdiff,altob,errfact
 real cp,rd,rv,kap,kapr,kap1,fv,pi,grav,deg2rad,rad2deg
 character(len=2) charfhr
 character(len=3) charnanal
 character(len=19) sid
 real, dimension(:), allocatable :: glats, glatspluspoles, ak, bk
 real, dimension(:), allocatable :: oblocx,oblocy,ob,zob,obtime,stdev,&
                                    anal_obz,stdevorig,anal_ob,biasob
 real, dimension(:,:), allocatable :: psg,zsg,analzs,anal_ob2
 real, dimension(:,:,:), allocatable :: tempg,tvg,qg,psig,pslg,&
       analpress,analtemp,analps,zg
 integer, allocatable, dimension(:) :: stattype,iuseob
 character(len=8), allocatable :: statid(:)
 namelist /nam_psop/nlevt,fhmin,fhmax,fhout,datestring,rlapse,&
                    nanals,obsfile,zthresh,errfact,delz_const

! Initialize mpi
 call MPI_Init(ierr)

! nproc is process number, numproc is total number of processes.
 call MPI_Comm_rank(MPI_COMM_WORLD,nproc,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,numproc,ierr)

 pi      = 4.*atan(1.0)
 rad2deg = 180./pi
 deg2rad = 1./rad2deg
 rd     = 2.8705e+2
 rv     = 4.6150e+2
 fv     = rv/rd-1.    ! used in virtual temperature equation 
 cp     = 1004.67
 grav   = 9.80665
 kap = rd/cp
 kapr = cp/rd
 kap1 = kap + 1.

 nchar=1; nreal=19
 iunit = 7
 iunitsig = 22
 iunit_nml = 912
 iunito = 9
 rlapse = 0.0065

 ! read namelist from file on all processors.
 zthresh = 9.9e31
 delz_const = 0.001 ! factor for adjusting ob error based on diff between station and model height
 nlevt = 2 ! use temp at level just above 1st level for pressure reduction.
 errfact = 1.0
 open(iunit_nml,file='psop.nml',form='formatted')
 read(iunit_nml,nam_psop)
 close(iunit_nml)
 if (nproc .eq. 0) write(6,nam_psop)


 if (numproc .lt. nanals) then
    print *,numproc,nanals
    print *,'warning, numproc too small!'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Finalize(ierr)
    stop
 end if

 ntimes = 1+((fhmax-fhmin)/fhout)

 nanal = nproc + 1

 read(datestring,'(i10)') idate

!==> read in obs data (on root process).
 
 if (nproc .eq. 0) then
 
 print *,trim(adjustl(obsfile))
 open(149,form="formatted",file=trim(adjustl(obsfile)))
 nobstot = 0
 do 
   read(149,9801,err=199,end=199) sid
   nobstot = nobstot + 1
 enddo
 199 continue

 end if
 
 call MPI_Bcast(nobstot,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_Status,ierr)

!==> allocate some arrays for obs and obs metadata.
  
 allocate(statid(nobstot))
 allocate(anal_ob(nobstot))
 allocate(biasob(nobstot))
 allocate(anal_ob2(nanals+1,nobstot))
 allocate(anal_obz(nobstot))
 allocate(stattype(nobstot))
 allocate(iuseob(nobstot))
 allocate(oblocx(nobstot))
 allocate(oblocy(nobstot))
 allocate(ob(nobstot))
 allocate(zob(nobstot))
 allocate(obtime(nobstot))
 allocate(stdev(nobstot))
 allocate(stdevorig(nobstot))
 
 if (nproc .eq. 0) then
 
 rewind(149)
 nobsh = 0
 biasob = 0.
 do nob=1,nobstot
!02907    181   23.44  64.33     5  -3.00   958.1   0.0000 1.00
      read(149,9801) statid(nob),stattype(nob),oblocx(nob),oblocy(nob),&
           izob,obtime(nob),ob(nob),bias,stdevorig(nob)
      if (bias .lt. 1.e20) biasob(nob) = bias
      stdev(nob) = errfact*stdevorig(nob)
      zob(nob) = izob
 9801 format(a8,1x,i3,1x,f7.2,1x,f6.2,1x,i5,1x,f6.2,&
             f7.1,1x,f8.4,1x,f4.2)
      if (oblocy(nob) .lt. 0.) nobsh = nobsh + 1
      if (oblocx(nob) .lt. 0.) oblocx(nob) = oblocx(nob) + 360.
 enddo
 print *, nobstot,' total obs'
 print *, nobsh,' total obs in SH'
 close(149)

 end if
 
 call MPI_Bcast(oblocx,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(oblocy,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(obtime,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(ob,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(zob,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(stdevorig,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(stdev,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(stattype,nobstot,MPI_INTEGER,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 call MPI_Bcast(biasob,nobstot,MPI_REAL,0, &
               MPI_COMM_WORLD,MPI_Status,ierr)
 
!==> convert ob location to radians.
  
 oblocx = deg2rad*oblocx
 oblocy = deg2rad*oblocy

 ntime = 0
 do nfhr=fhmin,fhmax,fhout
 ntime = ntime + 1

 write(charfhr,'(i2.2)') nfhr
 write(charnanal,'(i3.3)') nanal
 if (nanal .eq. nanals+1) then
    filenamein = "sfg_"//datestring//"_fhr"//charfhr//"_ensmean"
 else
    filenamein = "sfg_"//datestring//"_fhr"//charfhr//"_mem"//charnanal
 end if


 dset = open_dataset(trim(filenamein),errcode=iret)
 londim = get_dim(dset,'grid_xt'); nlons = londim%len
 latdim = get_dim(dset,'grid_yt'); nlats = latdim%len
 levdim = get_dim(dset,'pfull');   nlevs = levdim%len
 !print *,'filenamein,nlons,nlats,nlevs= ',trim(adjustl(filenamein)),nlons,nlats,nlevs

 if (iret .ne. 0) then
    print *,'error reading file ',iret,trim(filenamein)
    stop
 end if

 if (nfhr .eq. fhmin) then
    allocate(psg(nlons,nlats))
    allocate(zsg(nlons,nlats))
    allocate(tempg(nlons,nlats,nlevs))
    allocate(tvg(nlons,nlats,nlevs))
    allocate(qg(nlons,nlats,nlevs))
    allocate(zg(nlons,nlats,nlevs))
    allocate(pslg(nlons,nlats,nlevs))
    allocate(psig(nlons,nlats,nlevs+1))
    allocate(glats(nlats))
    allocate(glatspluspoles(nlats+2))
    allocate(ak(nlevs+1),bk(nlevs+1))
    allocate(analps(nlons+1,nlats+2,ntimes))
    allocate(analtemp(nlons+1,nlats+2,ntimes))
    allocate(analpress(nlons+1,nlats+2,ntimes))
    allocate(analzs(nlons+1,nlats+2))
 end if

 call read_vardata(dset,'tmp',tempg)
 call read_vardata(dset,'spfh',qg)
 tvg = tempg(:,:,nlevs:1:-1) * ( 1.0 + fv*qg(:,:,nlevs:1:-1) ) ! convert T to Tv, flip vertical
 call read_vardata(dset,'pressfc',psg)
 call read_vardata(dset,'hgtsfc',zsg)
 !call read_vardata(dset, 'grid_xt', glons)
 call read_vardata(dset, 'grid_yt', glats)
 call read_attribute(dset, 'ak', ak)
 call read_attribute(dset, 'bk', bk)
 call close_dataset(dset)
 do k=1,nlevs+1
    psig(:,:,k) = ak(k)+bk(k)*psg
 enddo
 do k=1,nlevs
   ! layer pressure from Phillips vertical interpolation.
   pslg(:,:,nlevs-k+1) = ((psig(:,:,k)**kap1-psig(:,:,k+1)**kap1)/&
                 (kap1*(psig(:,:,k)-psig(:,:,k+1))))**kapr
 enddo
 psg = psg/100. ! convert to mb (units of obs)
 pslg = pslg/100.
 psig = psig/100.
 !if (nproc .eq. 0 .and. nfhr .eq. fhmin) then
 !   print *,'psg min/max ',minval(psg),maxval(psg)
 !   do k=1,nlevs
 !      print *,k,minval(psig(:,:,k)),maxval(psig(:,:,k)),minval(pslg(:,:,k)),maxval(pslg(:,:,k)),minval(tvg(:,:,k)),maxval(tvg(:,:,k))
 !   enddo
 !endif

 ! add wraparound and pole points.
 call addpolewrap(psg,analps(:,:,ntime),nlons,nlats)
 call addpolewrap(zsg,analzs,nlons,nlats)
 call addpolewrap(tvg(:,:,nlevt),analtemp(:,:,ntime),nlons,nlats)
 call addpolewrap(pslg(:,:,nlevt),analpress(:,:,ntime),nlons,nlats)

 enddo ! nfhr

 !==> 0.5*pi-latitudes with poles included (used in bilinear interp routine).
 glats = deg2rad*glats
 glatspluspoles(1) = 0.
 glatspluspoles(nlats+2) = pi
 do j=2,nlats+1
   glatspluspoles(j) = 0.5*pi - glats(j-1)
 enddo

 !==> perform Benjamin and Miller reduction for each ob, compute ob priors.
 !    also do gross qc checks.
 nn = 0
 do nob=1,nobstot
    ! make sure ob location > -90 degrees.
    if (oblocy(nob) .gt. 0.5*pi+1.e-6 .or. oblocy(nob) .lt. -0.5*pi-1.e-6) then
       if (nproc .eq. 0) then
       print *,'WARNING: ob located outside domain',oblocy(nob)
       print *,'the ob latitude will be clipped to the nearest pole'
       end if
    end if
    if (oblocy(nob) .lt. -0.5*pi) oblocy(nob) = -0.5*pi
    if (oblocy(nob) .gt. 0.5*pi) oblocy(nob) = 0.5*pi
    ! longitudes are evenly spaced
    dxob = (0.5*float(nlons)*oblocx(nob)/pi)+1.
    ! gaussian latitudes are not.
    j=1
    dyob = 1.
    do 
       if (glatspluspoles(j) .ge. 0.5*pi-oblocy(nob)) exit
       j = j + 1
    enddo 
    if (j .gt. 1) dyob = float(j-1) + (0.5*pi-oblocy(nob)-glatspluspoles(j-1))/(glatspluspoles(j)-glatspluspoles(j-1))
    dtob = 1.+((real(fhmin)+obtime(nob))/real(fhout))
    call lintrp3(analtemp,anal_obt,&
                 dxob,dyob,dtob,nlons+1,nlats+2,ntimes)
    call lintrp3(analpress,anal_obp,&
                 dxob,dyob,dtob,nlons+1,nlats+2,ntimes)
    call lintrp2(analzs,anal_obz(nob),&
                 dxob,dyob,nlons+1,nlats+2)
    
    ! this is ob prior.
    call lintrp3(analps,anal_ob(nob),&
                 dxob,dyob,dtob,nlons+1,nlats+2,ntimes)
    ! adjust Hx to (perturbed) station height
    anal_ob(nob) = &
    preduce(anal_ob(nob),anal_obp,anal_obt,zob(nob),anal_obz(nob),rlapse,grav,rd)

    zdiff = zob(nob)-anal_obz(nob)
    ! adjust ob error based on diff between station and model height.
    ! (GSI uses 0.005 for delz_const?)
    stdev(nob) = stdev(nob) + delz_const*abs(zdiff)
    if ((obtime(nob) .ge. -3. .and. &
        obtime(nob) .le. 3.) .and. abs(zdiff) .lt. zthresh) then
        iuseob(nob) = 1
    else
        iuseob(nob) = 0
        nn = nn + 1
    end if
! gross error check.
    if (iuseob(nob) .eq. 1) then
       altob = palt(ob(nob),zob(nob),grav,rd)
       if (altob .lt. 850. .or. altob .gt. 1090.) then
          if (nproc .eq. numproc) print *,'failed gross error check',rad2deg*oblocx(nob),rad2deg*oblocy(nob),obtime(nob),ob(nob),zob(nob),anal_obz(nob),altob
          iuseob(nob)=-1
          nn = nn + 1
       end if
    end if   
 enddo
 if (nn .ne. 0) print *,nanal,nn,' failed gross qc check'

! distribute the ob error to all processors.
! first, gather back on last proc.
 if (nproc .eq. nanals) then
    do np=0,nanals-1
       call MPI_Recv(anal_ob2(np+1,:),nobstot,MPI_REAL,np, &
                     1,MPI_COMM_WORLD,MPI_Status,ierr)
    enddo
    do nob=1,nobstot
       ! replace hx computed from ens mean with ens mean hx.
       ensmean_ob = sum(anal_ob2(1:nanals,nob))/float(nanals)
       anal_ob2(nanals+1,nob) = ensmean_ob
       ! recenter hx ensemble about hx computed from ensemble mean.
       !do nanal=1,nanals
       !   anal_ob2(nanal,nob) = anal_ob2(nanal,nob) - ensmean_ob + anal_ob(nob)
       !enddo
       !anal_ob2(nanals+1,nob) = anal_ob(nob)
    enddo
 else
    !print *,'nproc',nproc,'min/max anal_ob',minval(anal_ob),maxval(anal_ob)
    call MPI_Send(anal_ob,nobstot,MPI_REAL,nanals, &
                  1,MPI_COMM_WORLD,ierr)
 end if

! now push back out to all other procs.
 call MPI_Bcast(anal_ob2,nobstot*(nanals+1),MPI_REAL,nanals, &
               MPI_COMM_WORLD,MPI_Status,ierr)

 call MPI_Barrier(MPI_COMM_WORLD,ierr)

 read(datestring,'(i10)') ianldate
 if (nanal .eq. nanals+1) then
    diag_file = "diag_conv_ps_ges."//datestring//"_ensmean.nc4"
 else
    diag_file = "diag_conv_ps_ges."//datestring//"_mem"//charnanal//".nc4"
 end if

 call nc_diag_init(diag_file, append=.false.)
 do nob=1,nobstot
 call nc_diag_header("date_time",ianldate )
 call nc_diag_metadata("Station_ID",              statid(nob)            )
 call nc_diag_metadata("Observation_Class",       "     ps"              )
 call nc_diag_metadata("Observation_Type",        stattype(nob)          )
 call nc_diag_metadata("Observation_Subtype",     0                      )
 call nc_diag_metadata("Latitude",                rad2deg*oblocy(nob)    )
 call nc_diag_metadata("Longitude",               rad2deg*oblocx(nob)    )
 call nc_diag_metadata("Station_Elevation",       zob(nob)               )
 call nc_diag_metadata("Pressure",                ob(nob)                )
 call nc_diag_metadata("Height",                  anal_obz(nob)          )
 call nc_diag_metadata("Time",                    obtime(nob)            )
 call nc_diag_metadata("Analysis_Use_Flag",       iuseob(nob)            )
 call nc_diag_metadata("Errinv_Input",            stdevorig(nob)         )
 call nc_diag_metadata("Errinv_Adjust",           stdev(nob)             )
 call nc_diag_metadata("Errinv_Final",            stdev(nob)             )
 call nc_diag_metadata("Observation",                   ob(nob)          )
 call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   ob(nob)-(anal_ob2(nanal,nob)+biasob(nob))   )
 call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", ob(nob)-anal_ob2(nanal,nob)   )
 ! for Jacobian, need index of ps in control vector
 !st_ind = sum(levels(1:ns3d)) + ps_ind
 !end_ind = sum(levels(1:ns3d)) + ps_ind
 !call nc_diag_header("jac_nnz", 1)
 !call nc_diag_header("jac_nind", 1)
 !call nc_diag_data2d("Observation_Operator_Jacobian_stind", st_ind)
 !call nc_diag_data2d("Observation_Operator_Jacobian_endind", st_ind)
 !call nc_diag_data2d("Observation_Operator_Jacobian_val", 1.)
 enddo
 call nc_diag_write

 if (nanal .eq. nanals+1) then
    open(9,form='formatted',file='psobs_prior.txt')
    do nob=1,nobstot
       if (stdev(nob) .gt. 99.99) stdev(nob) = 99.99
       write(9,9802) stattype(nob),rad2deg*oblocx(nob),rad2deg*oblocy(nob),&
               nint(zob(nob)),nint(anal_obz(nob)),obtime(nob),ob(nob),&
               anal_ob2(nanal,nob)+biasob(nob),stdevorig(nob),stdev(nob),iuseob(nob)
    enddo
    9802 format(i3,1x,f7.2,1x,f6.2,1x,i5,1x,i5,1x,f6.2,1x,f7.1,1x,&
                  f7.1,1x,f5.2,1x,f5.2,1x,i2)
    close(9)
 !else
 !   filename = 'psobs_mem'//charnanal//'.txt'
 !   print *,'write out ',trim(filename)
 !   open(9,form='formatted',file=filename)
 !   do nob=1,nobstot
 !      if (stdev(nob) .gt. 99.99) stdev(nob) = 99.99
 !      write(9,9802) stattype(nob),rad2deg*oblocx(nob),rad2deg*oblocy(nob),&
 !              nint(zob(nob)),nint(anal_obz(nob)),obtime(nob),ob(nob),&
 !              anal_ob2(nanal,nob)+biasob(nob),stdevorig(nob),stdev(nob),iuseob(nob)
 !   enddo
 !   close(9)
 end if
 close(iunito)

 call MPI_Barrier(MPI_COMM_WORLD,ierr)
 call MPI_Finalize(ierr)

end program psop

subroutine addpolewrap(fin,fout,nx,ny)
! add pole and wrap-around points to lon,lat array.
 integer j,nx,ny
 real fin(nx,ny),fout(nx+1,ny+2)
 do j=2,ny+1
    fout(1:nx,j) = fin(:,j-1)
 enddo
 fout(:,1) = sum(fin(:,1))/float(nx)
 fout(:,ny+2) = sum(fin(:,ny))/float(nx)
 fout(nx+1,:) = fout(1,:)
end subroutine addpolewrap

real function preduce(ps,tpress,t,zmodel,zob,rlapse,grav,rd)
! compute MAPS pressure reduction from model to station elevation
! See Benjamin and Miller (1990, MWR, p. 2100)
! uses 'effective' surface temperature extrapolated
! from virtual temp (tv) at tpress mb
! using std atmosphere (ICAO) lapse rate.
! ps - surface pressure to reduce.
! t - virtual temp. at pressure tpress.
! zmodel - model orographic height.
! zob - station height
   implicit none
   real, intent(in) :: t,tpress,zmodel,zob,ps,rlapse,grav,rd
   real t0,alpha
   alpha = rd*rlapse/grav
   t0 = t*(ps/tpress)**alpha ! eqn 4 from B&M
   preduce = ps*((t0 + rlapse*(zob-zmodel))/t0)**(1./alpha) ! eqn 1 from B&M
end function preduce

real function palt(ps,zs,grav,rd)
! compute QNH altimeter setting (in mb) given ps (in mb), zs (in m).
! see WAF, Sept 1998, p 833-850 (Pauley) section 2c
   real, parameter :: rlapse = 0.0065
   real, parameter :: t0 = 288.15
   real, parameter :: p0 = 1013.25
   real alpha
   real,intent(in) :: ps,zs,grav,rd
   alpha = rd*rlapse/grav
   palt = ps*(1.+(rlapse*zs/t0)*(p0/ps)**alpha)**(1./alpha)
end function palt

subroutine lintrp2(f,g,dx,dy,nx,ny)
                                                                      
! subprogram:    lintrp2      linear interpolation in two dimensions.
!
!
!   input argument list:
!     f        - input interpolator
!     dx,dy    - input x,y -coords of interpolation point (grid units)
!     nx,ny    - x,y-dimensions of interpolator grid
!
!   output argument list:
!     g        - output interpolatee
 
 integer, intent(in) :: nx,ny
 real, intent(in) :: f(nx,ny),dx,dy
 real, intent(out) :: g
 integer ix,iy,ixp,iyp
 real delx,dely
  

 ix=dx
 ix=max(1,min(ix,nx))
 iy=dy
 iy=max(1,min(iy,ny))
 ixp=ix+1
 ixp=min(ixp,nx)
 iyp=iy+1
 iyp=min(iyp,ny)
 delx=dx-ix
 dely=dy-iy
 g=f(ix,iy)*(1.-delx)*(1.-dely) &
  +f(ixp,iy)*delx*(1.-dely) &
  +f(ix,iyp)*(1.-delx)*dely &
  +f(ixp,iyp)*delx*dely

end subroutine lintrp2

subroutine lintrp3(f,g,dx,dy,dz,nx,ny,nz)
!                .      .    .                                        
! subprogram:    lintrp3      linear interpolation in three dimensions.
!
!   input argument list:
!     f        - input interpolator
!     dx,dy,gz    - input x,y,z -coords of interpolation point (grid units)
!     nx,ny,nz    - x,y,z-dimensions of interpolator grid
!
!     x is longitude, y is latitude and z is ln(sigma).
!     x is assumed to be regularly spaced, y and z aren't.
!     if z outside range of z(1) to z(nz), g=-9.9e31 is returned.
!
!   output argument list:
!     g        - output interpolatee
!
 integer, intent(in) :: nx,ny,nz
 real, intent(in) :: f(nx,ny,nz),dx,dy,dz
 real, intent(out) :: g
 integer ix,iy,ixp,iyp,iz,izp
 real delx,dely,delz

 ix=dx
 ix=max(1,min(ix,nx))
 iy=dy
 iy=max(1,min(iy,ny))
 iz=dz
 iz=max(1,min(iz,nz))
 ixp=ix+1
 ixp=min(ixp,nx)
 iyp=iy+1
 iyp=min(iyp,ny)
 izp=iz+1
 izp=min(izp,nz)
 delx=dx-ix
 dely=dy-iy
 delz=dz-iz
 g=f(ix,iy,iz)*(1.-delx)*(1.-dely)*(1.-delz)+ &
   f(ixp,iy,iz)*delx*(1.-dely)*(1.-delz)+ &
   f(ix,iyp,iz)*(1.-delx)*dely*(1.-delz)+ &
   f(ixp,iyp,iz)*delx*dely*(1.-delz)+ &
   f(ix,iy,izp)*(1.-delx)*(1.-dely)*delz+ &
   f(ixp,iy,izp)*delx*(1.-dely)*delz+ &
   f(ix,iyp,izp)*(1.-delx)*dely*delz+ &
   f(ixp,iyp,izp)*delx*dely*delz

end subroutine lintrp3

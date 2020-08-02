module m_gen_acqui
use m_string
use m_setup
use m_arrayop
use m_model, only:m

    !abbreviation:
    ! s,src : source
    ! r,rcv : receiver
    ! o : origin
    ! d : increment, spacing
    ! n : number of
    ! z : depth
    ! x : inline
    ! y : crossline
    ! pos : position
    ! comp : component

    private
    public :: gen_acqui_init, gen_acqui_nshot, gen_acqui_ntr, gen_acqui_nt, gen_acqui_dt

    character(:),allocatable :: acqui_type
    real :: osz, osx, osy !origin of src position
    real :: orz, orx, ory !origin of rcv position
    real :: dsz, dsx, dsy !src position sampling
    real :: drz, drx, dry !rcv position sampling
    type(t_string),dimension(:),allocatable :: scomp, rcomp
    integer :: nscomp, nrcomp !number of source & receiver components

    integer :: gen_acqui_nshot !number of shotgathers
    integer :: gen_acqui_ntr   !number of traces

    integer gen_acqui_nt !number of time steps
    real gen_acqui_dt !time interval

    contains
    
    subroutine gen_acqui_init
        character(:),allocatable :: tmp

        acqui_type=setup_get_char('ACQUI_TYPE',default='spread')

        tmp=setup_get_char('SOURCE_ORIGIN','OSRC');   read(tmp,*) osz, osx, osy
        tmp=setup_get_char('RECEIVER_ORIGIN','ORCV'); read(tmp,*) orz, orx, ory

        tmp=setup_get_char('SOURCE_SPACING','DSRC');  read(tmp,*) dsz, dsx, dsy
        tmp=setup_get_char('RECEIVER_SPACING','DRCV');read(tmp,*) drz, drx, dry
        
        ns=setup_get_int('NUMBER_SOURCE','NSRC')
        nr=setup_get_int('NUMBER_RECEIVER','NRCV')

        if(.not.m%is_cubic) then
            osy=0; ory=0
            dsy=0; dry=0
            nsy=1; nry=1
        endif

        scomp=partition(setup_get_char('SOURCE_COMPONENT',  'SRC_COMP',default='P'))
        rcomp=partition(setup_get_char('RECEIVER_COMPONENT','RCV_COMP',default='P'))
        nscomp=size(scomp)
        nrcomp=size(rcomp)

        gen_acqui_nshot=nscomp*ns
        gen_acqui_ntr  =nrcomp*nr
        
        gen_acqui_nt=setup_get_int('TIME_STEP','NT')
        gen_acqui_dt=setup_get_real('TIME_INTERVAL','DT')
        
    end subroutine
    
    subroutine gen_acqui_shotgather(ishot, &
                                    sc,sz,sx,sy, &
                                    rc,rz,rx,ry, &
                                    nt,dt)
        integer,intent(in) :: ishot
        character(4),intent(out)    :: sc !src comp
        real,intent(out)            :: sz,sx,sy !src pos
        character(4),dimension(gen_acqui_ntr),intent(out)   :: rc !rcv comp
        real,dimension(gen_acqui_ntr),intent(out)           :: rz,rx,ry !rcv pos

        character(:),allocatable :: sshot

        !src side
        sc = scomp( int(ishot/ns)+1 )%string

        ish = ishot- int(ishot/ns)*ns
        sz = osz + (ish-1)*dsz
        sx = osx + (ish-1)*dsx
        sy = osy + (ish-1)*dsy

        !rcv side
        do ic=0,nrcomp-1
            do ir=0,nr-1
                rc(1+ir+ic*nr) = rcomp(ic+1)%string
                rz(1+ir+ic*nr) = orz + ir*drz
                rx(1+ir+ic*nr) = orz + ir*drz
                ry(1+ir+ic*nr) = orz + ir*drz
            enddo
        enddo

        ! if(acqui_type=='spread') then
        !     call gen_acqui_spread
        ! elseif(acqui_type=='streamer') then
        !     call gen_acqui_streamer
        ! elseif(acqui_type=='irregularOBN') then
        !     call gen_acqui_irregularOBN
        ! elseif(acqui_type=='spread3D') then
        !     call gen_acqui_spread3D
        ! else
        !     call hud('Sorry, only spread, streamer, or irregular OBN acquisition geometries are implemented now')
        !     stop
        ! endif
        select case (acqui_type)
            !case ('spread')
            case ('streamer')
                rz=rz+sz
                rx=rx+sx
                ry=ry+sy
        end select

        !time
        nt=gen_acqui_nt
        dt=gen_acqui_dt

        !check source positions
        sshot='Shot# '//num2str(ishot)//' : '
        if(sz<0.) then
            sz = m%dz
            call hud(sshot//'sz < 0 above top of model! Shift sz = dz.',mpiworld%iproc)
        endif
        if(sz>(m%nz-1)*m%dz) then
            sz = (m%nz-2)*m%dz
            call hud(sshot//'sz below bottom of model! Shift sz = (nz-2)*dz.',mpiworld%iproc)
        endif
        if(sx<0.) then
            sx = m%dx
            call hud(sshot//'sx < 0 outside left bnd of model! Shift sx = dx.',mpiworld%iproc)
        endif
        if(sx>(m%nx-1)*m%dx) then
            sx = (m%nx-2)*m%dx
            call hud(sshot//'sx beyond right bnd of model! Shift sx = (nx-2)*dx.',mpiworld%iproc)
        endif
        if(sy<0.) then
            sy = m%dy
            call hud(sshot//'sy < 0 outside front end of model! Shift sy = dy.',mpiworld%iproc)
        endif
        if(sy>(m%ny-1)*m%dy) then
            sy = (m%ny-2)*m%dy
            call hud(sshot//'sy beyond rear end of model! Shift sy = (ny-2)*dy.',mpiworld%iproc)
        endif

        !check receiver positions
        sshot='Shot# '//num2str(ishot)//' has receiver(s) '
        if(any(rz<0.)) then
            where(rz<0.) rz=m%dz
            call hud(sshot//'rz < 0 above top of model! Shift rz = dz.',mpiworld%iproc)
        endif
        if(any(rz>(m%nz-1)*m%dz)) then
            where(rz>(m%nz-1)*m%dz) rz=(m%nz-2)*dz
            call hud(sshot//'rz below bottom of model! Shift rz = (nz-2)*dz.',mpiworld%iproc)
        endif
        if(any(rx<0.)) then
            where(rx<0.) rx=m%dx
            call hud(sshot//'rx < 0 outside left bnd of model! Shift rx = dx.',mpiworld%iproc)
        endif
        if(any(rx>(m%nx-1)*m%dx)) then
            where(rx>(m%nx-1)*m%dx) rx=(m%nx-2)*m%dx
            call hud(sshot//'rx beyond right bnd of model! Shift rx = (nx-2)*dx.',mpiworld%iproc)
        endif
        if(any(ry<0.)) then
            where(ry<0.) ry=m%dy
            call hud(sshot//'ry outside front end of model! Shift ry = dy.',mpiworld%iproc)
        endif
        if(any(ry>(m%ny-1)*m%dy)) then
            where(ry>(m%ny-1)*m%dy) ry=(m%ny-2)*m%dy
            call hud(sshot//'ry beyond rear end of model! Shift sy = (ny-2)*dy.',mpiworld%iproc)
        endif

    end subroutine

    ! subroutine gen_acqui_spread
        
    !     acqui%nsrc=ns*nscomp
    !     if(allocated(acqui%src))deallocate(acqui%src)
    !     allocate(acqui%src(acqui%nsrc))
        
    !     do k=1,nscomp

    !         sx=fsx; sy=fsy; sz=fsz
    !         do i=1,ns
    !             ii=i+(k-1)*ns
    !             acqui%src(ii)%x=sx;  sx=sx+dsx
    !             acqui%src(ii)%y=sy;  sy=sy+dsy
    !             acqui%src(ii)%z=sz;  sz=sz+dsz
    !             acqui%src(ii)%comp=src_comp(k)
                
    !             acqui%src(ii)%nrcv=nr*nrcomp
    !             if(allocated(acqui%src(ii)%rcv))deallocate(acqui%src(ii)%rcv)
    !             allocate(acqui%src(ii)%rcv(acqui%src(ii)%nrcv))
                
    !             do l=1,nrcomp

    !                 rx=frx; ry=fry; rz=frz
    !                 do j=1,nr
    !                     jj=j+(l-1)*nr
    !                     acqui%src(ii)%rcv(jj)%x=rx;  rx=rx+drx
    !                     acqui%src(ii)%rcv(jj)%y=ry;  ry=ry+dry
    !                     acqui%src(ii)%rcv(jj)%z=rz;  rz=rz+drz
    !                     acqui%src(ii)%rcv(jj)%comp=rcv_comp(l)
    !                 enddo

    !             enddo

    !         enddo

    !     enddo
        
    ! end subroutine
    
    ! subroutine gen_acqui_streamer
    !     real :: fsx,fsy,fsz,       lsx,lsy,lsz,       dsx,dsy,dsz
    !     real :: foffx,foffy,foffz, loffx,loffy,loffz, doffx,doffy,doffz
    !     integer :: ns, noff
        
    !     read(src_pos,*) fsz,  fsx,  fsy,   lsz,  lsx,  lsy,   ns
    !     read(rcv_pos,*) foffz,foffx,foffy, loffz,loffx,loffy, noff
        
    !     if(.not.m%is_cubic) then
    !         fsy=0; lsy=0;
    !         foffy=0; loffy=0;
    !     endif
        
    !     dsx=(lsx-fsx)/(ns-1);  doffx=(loffx-foffx)/(noff-1)
    !     dsy=(lsy-fsy)/(ns-1);  doffy=(loffy-foffy)/(noff-1)
    !     dsz=(lsz-fsz)/(ns-1);  doffz=(loffz-foffz)/(noff-1)
        
    !     acqui%nsrc=ns*nscomp
    !     if(allocated(acqui%src))deallocate(acqui%src)
    !     allocate(acqui%src(acqui%nsrc))
        
    !     do k=1,nscomp

    !         sx=fsx; sy=fsy; sz=fsz
    !         do i=1,ns
    !             ii=i+(k-1)*ns
    !             acqui%src(ii)%x=sx;  sx=sx+dsx
    !             acqui%src(ii)%y=sy;  sy=sy+dsy
    !             acqui%src(ii)%z=sz;  sz=sz+dsz
    !             acqui%src(ii)%comp=src_comp(k)
                
    !             acqui%src(ii)%nrcv=noff*nrcomp
    !             if(allocated(acqui%src(ii)%rcv))deallocate(acqui%src(ii)%rcv)
    !             allocate(acqui%src(ii)%rcv(acqui%src(ii)%nrcv))
                
    !             do l=1,nrcomp

    !                 rx=sx+foffx; ry=sy+foffy; rz=sz+foffz
    !                 do j=1,noff
    !                     jj=j+(l-1)*nr
    !                     acqui%src(ii)%rcv(jj)%x=rx;  rx=rx+drx
    !                     acqui%src(ii)%rcv(jj)%y=ry;  ry=ry+dry
    !                     acqui%src(ii)%rcv(jj)%z=rz;  rz=rz+drz
    !                     acqui%src(ii)%rcv(jj)%comp=rcv_comp(l)
    !                 enddo

    !             enddo

    !         enddo

    !     enddo
        
    ! end subroutine

    ! subroutine gen_acqui_irregularOBN
    !     integer :: ns, nr
    !     character(80) :: text

    !     !read sources
    !     open(13,file=src_pos,action='read')

    !         !count number of sources
    !         ns=0
    !         do
    !             read (13,*,iostat=msg) z,x,y
    !             if(msg/=0) exit
    !             ns=ns+1
    !         enddo
    !         if(mpiworld%is_master) write(*,*) 'Will read',ns,'source positions.'
            
    !         acqui%nsrc=ns*nscomp
    !         if(allocated(acqui%src))deallocate(acqui%src)
    !         allocate(acqui%src(acqui%nsrc))

    !     !read source positions
    !     rewind(13)
    !         i=1
    !         do
    !             read (13,*,iostat=msg) z,x,y
    !             if(msg/=0) exit
    !             acqui%src(i)%z=z
    !             acqui%src(i)%x=x
    !             acqui%src(i)%y=y
    !             acqui%src(i)%comp=src_comp(1)
    !             i=i+1
    !         enddo

    !     !duplicate source positions for other components
    !         do k=2,nscomp
    !             do i=1,ns
    !                 ii=i+(k-1)*ns
    !                 acqui%src(ii)%x=acqui%src(i)%x
    !                 acqui%src(ii)%y=acqui%src(i)%y
    !                 acqui%src(ii)%z=acqui%src(i)%z
    !                 acqui%src(ii)%comp=src_comp(k)
    !             enddo
    !         enddo

    !     close(13)


    !     !read receivers
    !     open(15,file=rcv_pos,action='read')

    !         !count number of receivers
    !         nr=0
    !         do
    !             read (15,*,iostat=msg) z,x,y
    !             if(msg/=0) exit
    !             nr=nr+1
    !         end do
    !         if(mpiworld%is_master) write(*,*) 'Will read',nr,'receiver positions.'

    !         acqui%src(:)%nrcv=nr*nrcomp
    !         do i=1,acqui%nsrc
    !             if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
    !             allocate(acqui%src(i)%rcv(acqui%src(i)%nrcv))
    !         enddo

    !     !read receiver positions
    !     rewind(15)
    !         j=1
    !         do
    !             read (15,*,iostat=msg) z,x,y
    !             if(msg/=0) exit
    !                 acqui%src(1)%rcv(j)%z=z
    !                 acqui%src(1)%rcv(j)%x=x
    !                 acqui%src(1)%rcv(j)%y=y
    !                 acqui%src(1)%rcv(j)%comp=rcv_comp(1)
    !                 j=j+1
    !         end do

    !     !duplicate receiver positions for other components
    !         do l=2,nrcomp
    !             do j=1,nr
    !                 jj=j+(l-1)*nr
    !                 acqui%src(1)%rcv(jj)%z=acqui%src(1)%rcv(j)%z
    !                 acqui%src(1)%rcv(jj)%x=acqui%src(1)%rcv(j)%x
    !                 acqui%src(1)%rcv(jj)%y=acqui%src(1)%rcv(j)%y
    !                 acqui%src(1)%rcv(jj)%comp=rcv_comp(l)
    !             enddo
    !         enddo

    !         do k=2,nscomp
    !             do i=1,ns
    !                 ii=i+(k-1)*ns
    !                 acqui%src(ii)%rcv(:)%z=acqui%src(i)%rcv(:)%z
    !                 acqui%src(ii)%rcv(:)%x=acqui%src(i)%rcv(:)%x
    !                 acqui%src(ii)%rcv(:)%y=acqui%src(i)%rcv(:)%y
    !                 acqui%src(ii)%rcv(:)%comp=rcv_comp(l)
    !             enddo
    !         enddo

    !     close(15)

    ! end subroutine

    ! subroutine gen_acqui_3Dspread
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !  Quadrilateral Geometry  !
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !
    !     !  IL direction ->
    !     ! XL  P1---------------P2
    !     ! dir -------------------
    !     ! |   -------------------
    !     ! v   P3---------------P4
    !     !
    !     !4 anchor points, P1-4, and number of points in IL and XL directions
    !     !should be read from setup (SOURCE_POSITION, RECEIVER_POSITION)
    !     !other points are interpolated from P1-4 and nil, nxl

    !     real :: fsx,fsy,fsz, frx,fry,frz
    !     real :: lsx,lsy,lsz, lrx,lry,lrz
        
    !     read(src_pos,*)  sz1,sx1,sy1, sz2,sx2,sy2, sz3,sx3,sy3, sz4,sx4,sy4, nsil, nsxl
    !     read(rcv_pos,*)  rz1,rx1,ry1, rz2,rx2,ry2, rz3,rx3,ry3, rz4,rx4,ry4, nril, nrxl
        
    !     if(.not.m%is_cubic) then
    !         sx3=sx1; sy3=sy1; sz3=sz1;
    !         sx4=sx2; sy4=sy2; sz4=sz2;  nsxl=1

    !         rx3=rx1; ry3=ry1; rz3=rz1;
    !         rx4=rx2; ry4=ry2; rz4=rz2;  nrxl=1
    !     endif

    !     ns=nsil*nsxl
    !     nr=nril*nrxl

    !     acqui%nsrc=ns*nscomp
    !     if(allocated(acqui%src))deallocate(acqui%src)
    !     allocate(acqui%src(acqui%nsrc))

    !     do ixl=1,nsxl

    !         fsx = ( sx1*(nsxl-ixl) + sx3*(ixl-1) )/(nsxl-1)
    !         fsy = ( sy1*(nsxl-ixl) + sy3*(ixl-1) )/(nsxl-1)
    !         fsz = ( sz1*(nsxl-ixl) + sz3*(ixl-1) )/(nsxl-1)

    !         lsx = ( sx2*(nsxl-ixl) + sx4*(ixl-1) )/(nsxl-1)
    !         lsy = ( sy2*(nsxl-ixl) + sy4*(ixl-1) )/(nsxl-1)
    !         lsz = ( sz2*(nsxl-ixl) + sz4*(ixl-1) )/(nsxl-1)

    !     do iil=1,nsil

    !         sx  = ( fsx*(nsil-iil) + lsx*(iil-1) )/(nsil-1)
    !         sy  = ( fsy*(nsil-iil) + lsy*(iil-1) )/(nsil-1)
    !         sz  = ( fsz*(nsil-iil) + lsz*(iil-1) )/(nsil-1)

    !         i   = iil + (ixl-1)*nsil

    !         acqui%src(i)%x=sx
    !         acqui%src(i)%y=sy
    !         acqui%src(i)%z=sz
    !         acqui%src(i)%comp=src_comp(1)
            

    !         acqui%src(i)%nrcv=nr*nrcomp
    !         if(allocated(acqui%src(i)%rcv))deallocate(acqui%src(i)%rcv)
    !         allocate(acqui%src(i)%rcv(acqui%src(i)%nrcv))
            
    !         do jxl=1,nrxl

    !             frx = ( rx1*(nrxl-jxl) + rx3*(jxl-1) )/(nrxl-1)
    !             fry = ( ry1*(nrxl-jxl) + ry3*(jxl-1) )/(nrxl-1)
    !             frz = ( rz1*(nrxl-jxl) + rz3*(jxl-1) )/(nrxl-1)

    !             lrx = ( rx2*(nrxl-jxl) + rx4*(jxl-1) )/(nrxl-1)
    !             lry = ( ry2*(nrxl-jxl) + ry4*(jxl-1) )/(nrxl-1)
    !             lrz = ( rz2*(nrxl-jxl) + rz4*(jxl-1) )/(nrxl-1)

    !         do jil=1,nril

    !             rx  = ( frx*(nril-jil) + lrx*(jil-1) )/(nril-1)
    !             ry  = ( fry*(nril-jil) + lry*(jil-1) )/(nril-1)
    !             rz  = ( frz*(nril-jil) + lrz*(jil-1) )/(nril-1)

    !             j   = jil + (jxl-1)*nril

    !             acqui%src(i)%rcv(j)%x=rx
    !             acqui%src(i)%rcv(j)%y=ry
    !             acqui%src(i)%rcv(j)%z=rz
    !             acqui%src(i)%rcv(j)%z=rcv_comp(1)
                
    !         enddo

    !         enddo

    !     enddo

    !     enddo

    !     !duplicate positions for other components
    !     do l=2,nrcomp
    !         do j=1,nr
    !             jj=j+(l-1)*nr
    !             acqui%src(1)%rcv(jj)%z=acqui%src(1)%rcv(j)%z
    !             acqui%src(1)%rcv(jj)%x=acqui%src(1)%rcv(j)%x
    !             acqui%src(1)%rcv(jj)%y=acqui%src(1)%rcv(j)%y
    !         enddo
    !     enddo

    !     do k=2,nscomp
    !         do i=1,ns
    !             ii=i+(k-1)*ns
    !             acqui%src(ii)%rcv(:)%z=acqui%src(i)%rcv(:)%z
    !             acqui%src(ii)%rcv(:)%x=acqui%src(i)%rcv(:)%x
    !             acqui%src(ii)%rcv(:)%y=acqui%src(i)%rcv(:)%y
    !         enddo
    !     enddo

    ! end subroutine
     
end

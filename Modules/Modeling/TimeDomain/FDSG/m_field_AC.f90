module m_field
use m_mpienv
use m_arrayop
use m_sysio
use m_model, only:m
use m_shot, only:shot
use m_computebox, only: cb

    private fdcoeff_o4,fdcoeff_o8,c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z,c4x,c4y,c4z
    private k_x,k_y,k_z,npower,rcoef
    public
    
    !FD coeff
    real,dimension(2),parameter :: fdcoeff_o4 = [1.125,      -1./24.]
    real,dimension(4),parameter :: fdcoeff_o8 = [1225./1024, -245./3072., 49./5120., -5./7168]
    
    real :: c1x, c1y, c1z
    real :: c2x, c2y, c2z
    real :: c3x, c3y, c3z
    real :: c4x, c4y, c4z
    
    !CPML
    real,parameter :: k_x = 1.
    real,parameter :: k_y = 1.
    real,parameter :: k_z = 1.
    real,parameter :: npower=2.
    real,parameter :: rcoef=0.001
    
    !local models
    real,dimension(:),allocatable ::  buox, buoy, buoz, kpa, invkpa
    real,dimension(:,:,:),pointer :: pbuox,pbuoy,pbuoz,pkpa,pinvkpa !aliases
    
    !fields
    type t_field
        !flat arrays
        real,dimension(:),allocatable :: vx,vy,vz,p
        real,dimension(:),allocatable :: prev_vx,prev_vy,prev_vz,prev_p !for time derivatives
        real,dimension(:),allocatable :: cpml_dvx_dx,cpml_dvy_dy,cpml_dvz_dz, cpml_dp_dz,cpml_dp_dx,cpml_dp_dy !for cpml
    end type
    
    !hicks interpolation
    logical :: if_hicks
    
    !gradient info
    character(*),parameter :: gradient_info='kpa-rho'
    
    contains
    
    !========= use before propagation =================
    subroutine field_print_info
        !modeling method
        call hud('WaveEq : Time-domain ACoustic isotropic modeling')
        call hud('Staggered-grid Finite-difference method')
        call hud('O(x4,t2) stencil')
        
        !stencil constant
        if(mpiworld%is_master) then
            write(*,*) 'Coeff:',fdcoeff_o4
        endif
        c1x=fdcoeff_o4(1)/m%dx; c1y=fdcoeff_o4(1)/m%dy; c1z=fdcoeff_o4(1)/m%dz
        c2x=fdcoeff_o4(2)/m%dx; c2y=fdcoeff_o4(2)/m%dy; c2z=fdcoeff_o4(2)/m%dz
        
        !hicks point interpolation
        if_hicks=get_setup_logical('IF_HICKS',default=.true.)
        
    end subroutine
    
    subroutine check_model  !not required
        
    end
    
    subroutine check_discretization
        !grid dispersion condition
        if (5.*m%dmax > cb%velmin/shot%src%fpeak/2.) then  !FDTDo4 rule
            write(*,*) 'WARNING: Shot# '//shot%cindex//' can have grid dispersion!'
            write(*,*) 'Shot# '//shot%cindex//' 5*dx, velmin, fpeak:',5.*m%dmax, cb%velmin,shot%src%fpeak
        endif
        
        !CFL condition
        if (m%is_cubic) then
            cfl = cb%velmax*shot%src%dt/m%dmin / (sqrt(3.)*(sum(abs(fdcoeff_o4))))! ~0.494
        else
            cfl = cb%velmax*shot%src%dt/m%dmin / (sqrt(2.)*(sum(abs(fdcoeff_o4))))! ~0.606
        end if
        if(mpiworld%is_master) write(*,*) 'CFL value:',CFL
        
        if(cfl>1.) then
            write(*,*) 'ERROR: CFL > 1 on shot# '//shot%cindex//'!'
            write(*,*) 'Shot# '//shot%cindex//' velmax, dt, dx:',cb%velmax,shot%src%dt,m%dmin
            stop
        endif
        
    end subroutine
    
    !========= use inside propagation =================
    
    subroutine init_field_localmodel
        call alloc(buox,cb%n)
        call alloc(buoy,cb%n)
        call alloc(buoz,cb%n)
        call alloc(kpa, cb%n)
        call alloc(invkpa,cb%n)
        
        call inflate_alias(buox,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pbuox)
        call inflate_alias(buoy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pbuoy)
        call inflate_alias(buoz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pbuoz)
        call inflate_alias(kpa, [cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pkpa)
        call inflate_alias(invkpa,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pinvkpa)
        
        pkpa=cb%vp*cb%vp*cb%rho
        invkpa=1./kpa
        
        pbuoz(cb%ifz,:,:)=1./cb%rho(cb%ifz,:,:)
        pbuox(:,cb%ifx,:)=1./cb%rho(:,cb%ifx,:)
        pbuoy(:,:,cb%ify)=1./cb%rho(:,:,cb%ify)
        
        do iz=cb%ifz+1,cb%ilz
            pbuoz(iz,:,:)=0.5/cb%rho(iz,:,:)+0.5/cb%rho(iz-1,:,:)
        enddo
        
        do ix=cb%ifx+1,cb%ilx
            pbuox(:,ix,:)=0.5/cb%rho(:,ix,:)+0.5/cb%rho(:,ix-1,:)
        enddo

        do iy=cb%ify+1,cb%ily
            pbuoy(:,:,iy)=0.5/cb%rho(:,:,iy)+0.5/cb%rho(:,:,iy-1)
        enddo
        
    end subroutine
    
    subroutine init_field(f,if_save_previous)
        type(t_field) :: f
        logical,optional :: if_save_previous
        
        call alloc(f%vx,cb%n)
        call alloc(f%vy,cb%n)
        call alloc(f%vz,cb%n)
        call alloc(f%p,cb%n)
        
        if(present(if_save_previous)) then
        if(if_save_previous) then
            !wavefield at previous incident timestep for gradient computation
            call alloc(f%prev_vx,cb%n)
            call alloc(f%prev_vy,cb%n)
            call alloc(f%prev_vz,cb%n)
            call alloc(f%prev_p,cb%n)
        endif
        endif
        
        call alloc(f%cpml_dvx_dx,cb%n)
        call alloc(f%cpml_dvy_dy,cb%n)
        call alloc(f%cpml_dvz_dz,cb%n)
        call alloc(f%cpml_dp_dx,cb%n)
        call alloc(f%cpml_dp_dy,cb%n)
        call alloc(f%cpml_dp_dz,cb%n)
        
    end subroutine
    
    subroutine check_field(f,name)
        type(t_field) :: f
        character(*) :: name
        
        if(mpiworld%is_master) write(*,*) name//' sample values:',minval(f%p),maxval(f%p)
        if(any(f%p-1.==f%p)) then
            write(*,*) 'ERROR: '//name//' values become +-Infinity !!'
            stop
        endif
        if(any(isnan(f%p))) then
            write(*,*) 'ERROR: '//name//' values become NaN !!'
            stop
        endif
        
    end subroutine
    
    subroutine field_cpml_reinitialize(f)
        type(t_field) :: f
        f%cpml_dvx_dx=0.
        f%cpml_dvy_dy=0.
        f%cpml_dvz_dz=0.
        f%cpml_dp_dx=0.
        f%cpml_dp_dy=0.
        f%cpml_dp_dz=0.
    end subroutine
    
    subroutine field_save_previous(f)
        type(t_field) :: f
        f%prev_vx=f%vx
        f%prev_vy=f%vy
        f%prev_vz=f%vz
        f%prev_p=f%p
    end subroutine
    
    subroutine write_field(iunit,f,if_difference)
        type(t_field) :: f
        logical,optional :: if_difference
        if(present(if_difference)) then
            if(if_difference) then
                write(iunit) f%prev_vz-f%vz
            else
                write(iunit) f%vz
            endif
        else
            write(iunit) f%vz
        endif
    end subroutine
    
    !========= forward propagation =================
    !WE: du_dt = MDu
    !u=[vx vy vz p]^T, p=(sxx+syy+szz)/3=sxx+syy+szz
    !  [diag3(b)  0  ]    [0   0   0   dx]
    !M=[   0     kpa ], D=|0   0   0   dy|
    !  (b=1/rho)          |0   0   0   dz|
    !                     [dx  dy  dz  0 ]
    !
    !Discretization (staggered grid in space and time):
    !  (grid index)     (real index)
    !  s(iz,ix,iy):=   s[iz,ix,iy]^it+0.5,it+1.5,...
    ! vx(iz,ix,iy):=  vx[iz,ix-0.5,iy]^it,it+1,...
    ! vy(iz,ix,iy):=  vy[iz,ix,iy-0.5]^it,it+1,...
    ! vz(iz,ix,iy):=  vy[iz-0.5,ix,iy]^it,it+1,...
    
    !add RHS to v^it
    subroutine put_velocities(time_dir,it,w,f)
        integer :: time_dir,it
        real :: w
        type(t_field) :: f
        
        real,dimension(:,:,:),pointer :: pvx,pvy,pvz
        
        call inflate_alias(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvx)
        call inflate_alias(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvy)
        call inflate_alias(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvz)
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (2)
                pvx(ifz:ilz,ifx:ilx,ify:ily) = pvx(ifz:ilz,ifx:ilx,ify:ily) + source_term *pbuox(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
                case (3)
                pvy(ifz:ilz,ifx:ilx,ify:ily) = pvy(ifz:ilz,ifx:ilx,ify:ily) + source_term *pbuoy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
                case (4)
                pvz(ifz:ilz,ifx:ilx,ify:ily) = pvz(ifz:ilz,ifx:ilx,ify:ily) + source_term *pbuoz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff
                
            end select
            
        else
            select case (shot%src%icomp)
                case (2) !horizontal x force on vx[iz,ix-0.5,iy]
                pvx(iz,ix,iy) = pvx(iz,ix,iy) + source_term*pbuox(iz,ix,iy)
                
                case (3) !horizontal y force on vy[iz,ix,iy-0.5]
                pvy(iz,ix,iy) = pvy(iz,ix,iy) + source_term*pbuoy(iz,ix,iy)
                
                case (4) !vertical force     on vz[iz-0.5,ix,iy]
                pvz(iz,ix,iy) = pvz(iz,ix,iy) + source_term*pbuoz(iz,ix,iy)
                
            end select
            
        endif
        
    end subroutine
    
    !v^it -> v^it+1 by FD of s^it+0.5
    subroutine update_velocities(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        ifz=bl(1)+2
        ilz=bl(2)-1
        ifx=bl(3)+2
        ilx=bl(4)-1
        if(m%is_cubic) then !3D
            ify=bl(5)+2
            ily=bl(6)-1
        else !2D
            ify=1
            ily=1
        endif
        
        if(m%if_freesurface) ifz=max(ifz,1)
        
        dir_dt=time_dir*shot%src%dt
        
        dp_dx=0.;dp_dy=0.;dp_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,&
        !$omp         izm2_ix_iy,izm1_ix_iy,iz_ix_iy,izp1_ix_iy,&
        !$omp         iz_ixm2_iy,iz_ixm1_iy,iz_ixp1_iy,&
        !$omp         iz_ix_iym2,iz_ix_iym1,iz_ix_iyp1,&
        !$omp         dp_dx,dp_dy,dp_dz)
        !$omp do schedule(dynamic)
        do iy=ify,ily
        do ix=ifx,ilx
        do iz=ifz,ilz
            i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(iy-cb%ify)*nz*nx+1
            
            izm2_ix_iy=i-2  !iz-2,ix,iy
            izm1_ix_iy=i-1  !iz-1,ix,iy
            iz_ix_iy  =i    !iz,ix,iy
            izp1_ix_iy=i+1  !iz+1,ix,iy
            
            iz_ixm2_iy=i-2*nz  !iz,ix-2,iy
            iz_ixm1_iy=i  -nz  !iz,ix-1,iy
            iz_ixp1_iy=i  +nz  !iz,ix+1,iy
            
            iz_ix_iym2=i-2*nz*nx  !iz,ix,iy-2
            iz_ix_iym1=i  -nz*nx  !iz,ix,iy-1
            iz_ix_iyp1=i  +nz*nx  !iz,ix,iy+1
            
            dp_dx= c1x*(f%p(iz_ix_iy)-f%p(iz_ixm1_iy)) +c2x*(f%p(iz_ixp1_iy)-f%p(iz_ixm2_iy))
            if(m%is_cubic) &
            dp_dy= c1y*(f%p(iz_ix_iy)-f%p(iz_ix_iym1)) +c2y*(f%p(iz_ix_iyp1)-f%p(iz_ix_iym2))
            dp_dz= c1z*(f%p(iz_ix_iy)-f%p(izm1_ix_iy)) +c2z*(f%p(izp1_ix_iy)-f%p(izm2_ix_iy))
            
            !cpml
            f%cpml_dp_dx(iz_ix_iy)= cb%a_x(ix)*dp_dx + cb%b_x(ix)*f%cpml_dp_dx(iz_ix_iy)
            dp_dx=dp_dx*k_x + f%cpml_dp_dx(iz_ix_iy)

            f%cpml_dp_dy(iz_ix_iy)= cb%a_y(iy)*dp_dy + cb%b_y(iy)*f%cpml_dp_dy(iz_ix_iy) 
            dp_dy=dp_dy*k_y + f%cpml_dp_dy(iz_ix_iy)
            
            f%cpml_dp_dz(iz_ix_iy)= cb%a_z(iz)*dp_dz + cb%b_z(iz)*f%cpml_dp_dz(iz_ix_iy)
            dp_dz=dp_dz*k_z + f%cpml_dp_dz(iz_ix_iy)
            
            !velocity
            f%vx(iz_ix_iy)=f%vx(iz_ix_iy) + dir_dt*buox(iz_ix_iy)*dp_dx
            f%vy(iz_ix_iy)=f%vy(iz_ix_iy) + dir_dt*buoy(iz_ix_iy)*dp_dy
            f%vz(iz_ix_iy)=f%vz(iz_ix_iy) + dir_dt*buoz(iz_ix_iy)*dp_dz
        enddo
        enddo
        enddo
        !$omp enddo 
        !$omp end parallel
        
        
        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,iy] level
        !so symmetric mirroring: vz[0.5]=vz[1.5], ie. vz(1,ix,iy)=vz(2,ix,iy) -> dp(1,ix,iy)=0.
        if(m%if_freesurface) then
            !$omp parallel default (shared)&
            !$omp private(ix,iy,i)
            !$omp do schedule(dynamic)
            do iy=ify,ily
            do ix=ifx,ilx
                i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy
                
                f%vz(i)=f%vz(i+1)
            enddo
            enddo
            !$omp enddo
            !$omp end parallel
        endif
        
    end subroutine
    
    !add RHS to s^it+0.5
    subroutine put_stresses(time_dir,it,w,f)
        integer :: time_dir
        real :: w
        type(t_field) :: f
        
        real,dimension(:,:,:),pointer :: pp
        
        call inflate_alias(f%p,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pp)
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
        
        source_term=time_dir*w
        
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                pp(ifz:ilz,ifx:ilx,ify:ily) = pp(ifz:ilz,ifx:ilx,ify:ily) + source_term *shot%src%interp_coeff
            end select
        
        else
            select case (shot%src%icomp)
                case (1) !explosion on s[iz,ix,iy]
                pp(iz,ix,iy) = pp(iz,ix,iy) + source_term
            end select
            
        endif
        
    end subroutine
    
    !s^it+0.5 -> s^it+1.5 by FD of v^it+1
    subroutine update_stresses(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        nz=cb%nz
        nx=cb%nx
        ny=cb%ny
        
        ifz=bl(1)+1
        ilz=bl(2)-2
        ifx=bl(3)+1
        ilx=bl(4)-2
        if(m%is_cubic) then !3D
            ify=bl(5)+1
            ily=bl(6)-2
        else !2D
            ify=1
            ily=1
        endif
        
        if(m%if_freesurface) ifz=max(ifz,1)
        
        dir_dt=time_dir*shot%src%dt
        
        dvx_dx=0.;dvy_dy=0.;dvz_dz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,&
        !$omp         izm1_ix_iy,iz_ix_iy,izp1_ix_iy,izp2_ix_iy,&
        !$omp         iz_ixm1_iy,iz_ixp1_iy,iz_ixp2_iy,&
        !$omp         iz_ix_iym1,iz_ix_iyp1,iz_ix_iyp2,&
        !$omp         dvx_dx,dvy_dy,dvz_dz)
        !$omp do schedule(dynamic)
        do iy=ify,ily
        do ix=ifx,ilx
        do iz=ifz,ilz
            i=(iz-cb%ifz)+(ix-cb%ifx)*nz+(iy-cb%ify)*nz*nx+1
            
            izm1_ix_iy=i-1  !iz-1,ix,iy
            iz_ix_iy  =i    !iz,ix,iy
            izp1_ix_iy=i+1  !iz+1,ix,iy
            izp2_ix_iy=i+2  !iz+2,ix,iy
            
            iz_ixm1_iy=i  -nz  !iz,ix-1,iy
            iz_ixp1_iy=i  +nz  !iz,ix+1,iy
            iz_ixp2_iy=i +2*nz !iz,ix+2,iy
            
            iz_ix_iym1=i  -nz*nx  !iz,ix,iy-1
            iz_ix_iyp1=i  +nz*nx  !iz,ix,iy+1
            iz_ix_iyp2=i+2*nz*nx  !iz,ix,iy+2
                
                
            dvx_dx= c1x*(f%vx(iz_ixp1_iy)-f%vx(iz_ix_iy))  +c2x*(f%vx(iz_ixp2_iy)-f%vx(iz_ixm1_iy))
            if(m%is_cubic) &
            dvy_dy= c1y*(f%vy(iz_ix_iyp1)-f%vy(iz_ix_iy))  +c2y*(f%vy(iz_ix_iyp2)-f%vy(iz_ix_iym1))
            dvz_dz= c1z*(f%vz(izp1_ix_iy)-f%vz(iz_ix_iy))  +c2z*(f%vz(izp2_ix_iy)-f%vz(izm1_ix_iy))
            
            !cpml
            f%cpml_dvx_dx(iz_ix_iy)=cb%b_x(ix)*f%cpml_dvx_dx(iz_ix_iy)+cb%a_x(ix)*dvx_dx
            dvx_dx=dvx_dx*k_x + f%cpml_dvx_dx(iz_ix_iy)
            
            f%cpml_dvy_dy(iz_ix_iy)=cb%b_y(iy)*f%cpml_dvy_dy(iz_ix_iy)+cb%a_y(iy)*dvy_dy
            dvy_dy=dvy_dy*k_y + f%cpml_dvy_dy(iz_ix_iy)
            
            f%cpml_dvz_dz(iz_ix_iy)=cb%b_z(iz)*f%cpml_dvz_dz(iz_ix_iy)+cb%a_z(iz)*dvz_dz
            dvz_dz=dvz_dz*k_z + f%cpml_dvz_dz(iz_ix_iy)
            
            !pressure
            f%p(iz_ix_iy) = f%p(iz_ix_iy) + dir_dt * kpa(iz_ix_iy)*(dvx_dx+dvy_dy+dvz_dz)
        enddo
        enddo
        enddo
        !$omp enddo 
        !$omp end parallel
        
        
        !apply free surface boundary condition if needed
        !free surface is located at [1,ix,iy] level
        !so explicit boundary condition: p(1,ix,iy)=0
        !and antisymmetric mirroring: p(0,ix,iy)=-p(2,ix,iy) -> vz(2,ix,iy)=vz(1,ix,iy)
        if(m%if_freesurface) then
            !$omp parallel default (shared)&
            !$omp private(ix,iy,i)
            !$omp do schedule(dynamic)
            do iy=ify,ily
            do ix=ifx,ilx
                i=(1-cb%ifz) + (ix-cb%ifx)*nz + (iy-cb%ify)*nz*nx +1 !iz=1,ix,iy 
                
                f%p(i)=0.
                
                f%p(i-1)=-f%p(i+1)
            enddo
            enddo
            !$omp enddo
            !$omp end parallel
        endif
        
    end subroutine
    
    !get v^it+1 or s^it+1.5
    subroutine get_field(f,seismo)
        type(t_field) :: f
        real,dimension(*) :: seismo
        
        real,dimension(:,:,:),pointer :: pp,pvx,pvy,pvz
        
        call inflate_alias(f%p,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pp)
        call inflate_alias(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvx)
        call inflate_alias(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvy)
        call inflate_alias(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvz)
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1)
                    seismo(ircv)=sum(  pp(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (2)
                    seismo(ircv)=sum( pvx(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (3)
                    seismo(ircv)=sum( pvy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                    case (4)
                    seismo(ircv)=sum( pvz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff )
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !p[iz,ix,iy]
                    seismo(ircv)= pp(iz,ix,iy)
                    case (2) !vx[iz,ix-0.5,iy]
                    seismo(ircv)=pvx(iz,ix,iy)
                    case (3) !vy[iz,ix,iy-0.5]
                    seismo(ircv)=pvy(iz,ix,iy)
                    case (4) !vz[iz-0.5,ix,iy]
                    seismo(ircv)=pvz(iz,ix,iy)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    
    !========= adjoint propagation =================
    !WE:        du_dt   = MDu   + src
    !->        Ndu_dt   = Du    + Nsrc, N=M^-1
    !Adjoint:  Ndv_dt^T = D^Tv  + adjsrc
    !->        -dv_dt   =-MDv   + Mr
    !     (reverse time) (negative derivatives)
    !
    !Discrete form (staggered grid in space and time, 2D as example):
    ! N  [vx^it+1  ] [vx^it    ]   [ 0     0    dx(-)][vx^it+1  ]
    !----|vz^it+1  |-|vz^it    | = | 0     0    dz(-)||vz^it+1  |  +Nsrc
    ! dt [ s^it+1.5] [ s^it+0.5]   [dx(+) dz(+)  0   ][ s^it+0.5]
    !dx(-):=a1(s(ix)-s(ixm1))+a2(s(ixp1)-s(ixm2))  (O(x4))
    !dx(+):=a1(v(ixp1)-v(ix))+a2(v(ixp2)-v(ixm1))  (O(x4))
    !Adjoint:
    ! N  [vx^it    ] [vx^it+1  ]   [ 0       0      dx(+)^T][vx^it+1  ]
    !----|vz^it    |-|vz^it+1  | = | 0       0      dz(+)^T||vz^it+1  |  +Nsrc
    ! dt [ s^it+0.5] [ s^it+1.5]   [dx(-)^T dz(-)^T  0     ][ s^it+0.5]
    !dx(+)^T=a1(s(ixm1)-s(ix))+a2(s(ixm2)-s(ixp1))=-dx(-)
    !dx(-)^T=a1(v(ix)-v(ixp1))+a2(v(ixm1)-v(ixp2))=-dx(+)
    !-dx,-dz can be regarded as -dt, thus allowing to use same code with flip of dt sign.
    !transposes of time derivatives centered on v^it+0.5 and s^it+1
    !so v^it+1   - v^it     => v^it - v^it+1
    !   s^it+1.5 - s^it+0.5 => s^it+0.5 - s^it+1.5
    !
    !In each time step:
    !d=RGAsrc, src:source term, A:put source, G:propagator, R:sampling at receivers
    !d=[vx vy vz p]^T, R=[I  0   0 ], A=[diag3(b) 0], src=[fx fy fz p]^T
    !                    |0 2/3 1/3|    [   0     1]
    !                    [   0    1]
    !Adjoint: src=A^T N G^T M R^T d
    !M R^T:put adjoint sources, 
    !G^T:  adjoint propagator,
    !A^T N:get adjoint fields
    
    
    !add RHS to s^it+1.5
    subroutine put_stresses_adjoint(time_dir,it,adjsource,f)
        integer :: time_dir
        real,dimension(*) :: adjsource
        type(t_field) :: f
        
        real,dimension(:,:,:),pointer :: pp
        
        call inflate_alias(f%p,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pp)
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            tmp=adjsource(ircv)
        
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    pp(ifz:ilz,ifx:ilx,ify:ily) = pp(ifz:ilz,ifx:ilx,ify:ily) &
                        +tmp* pkpa(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (1) !pressure adjsource
                    !p[iz,ix,iy]
                    pp(iz,ix,iy) = pp(iz,ix,iy)  &
                        +tmp* pkpa(iz,ix,iy)
                end select
                
            endif
            
        enddo
        
    end subroutine
    
    !s^it+1.5 -> s^it+0.5 by FD^T of v^it+1
    subroutine update_stresses_adjoint(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        !same code with -dt
        call update_stresses(time_dir,it,f,bl)
        
    end subroutine
    
    !add RHS to v^it+1
    subroutine put_velocities_adjoint(time_dir,it,adjsource,f)
        integer :: time_dir
        real,dimension(*) :: adjsource
        type(t_field) :: f
        
        real,dimension(:,:,:),pointer :: pvx,pvy,pvz
        
        call inflate_alias(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvx)
        call inflate_alias(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvy)
        call inflate_alias(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvz)
        
        do ircv=1,shot%nrcv
            ifz=shot%rcv(ircv)%ifz; iz=shot%rcv(ircv)%iz; ilz=shot%rcv(ircv)%ilz
            ifx=shot%rcv(ircv)%ifx; ix=shot%rcv(ircv)%ix; ilx=shot%rcv(ircv)%ilx
            ify=shot%rcv(ircv)%ify; iy=shot%rcv(ircv)%iy; ily=shot%rcv(ircv)%ily
            
            tmp=adjsource(ircv)
            
            if(if_hicks) then
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    pvx(ifz:ilz,ifx:ilx,ify:ily) = pvx(ifz:ilz,ifx:ilx,ify:ily) + tmp*pbuox(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                    
                    case (3) !horizontal y adjsource
                    pvy(ifz:ilz,ifx:ilx,ify:ily) = pvy(ifz:ilz,ifx:ilx,ify:ily) + tmp*pbuoy(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                    
                    case (4) !horizontal z adjsource
                    pvz(ifz:ilz,ifx:ilx,ify:ily) = pvz(ifz:ilz,ifx:ilx,ify:ily) + tmp*pbuoz(ifz:ilz,ifx:ilx,ify:ily) *shot%rcv(ircv)%interp_coeff
                end select
                
            else
                select case (shot%rcv(ircv)%icomp)
                    case (2) !horizontal x adjsource
                    !vx[ix-0.5,iy,iz]
                    pvx(iz,ix,iy) = pvx(iz,ix,iy) + tmp*pbuox(iz,ix,iy)
                    
                    case (3) !horizontal y adjsource
                    !vy[ix,iy-0.5,iz]
                    pvy(iz,ix,iy) = pvy(iz,ix,iy) + tmp*pbuoy(iz,ix,iy)
                    
                    case (4) !horizontal z adjsource
                    !vz[ix,iy,iz-0.5]
                    pvz(iz,ix,iy) = pvz(iz,ix,iy) + tmp*pbuoz(iz,ix,iy)
                end select
                
            endif
        
        enddo
        
    end subroutine
    
    !v^it+1 -> v^it by FD^T of s^it+0.5
    subroutine update_velocities_adjoint(time_dir,it,f,bl)
        integer :: time_dir, it
        type(t_field) :: f
        integer,dimension(6) :: bl
        
        !same code with -dt
        call update_velocities(time_dir,it,f,bl)
        
    end subroutine
    
    !get v^it or s^it+0.5
    subroutine get_field_adjoint(f,w)
        type(t_field) :: f
        real w
        
        real,dimension(:,:,:),pointer :: pp,pvx,pvy,pvz
        
        call inflate_alias(f%p,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pp)
        call inflate_alias(f%vx,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvx)
        call inflate_alias(f%vy,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvy)
        call inflate_alias(f%vz,[cb%ifz,cb%ilz],[cb%ifx,cb%ilx],[cb%ify,cb%ily],pvz)
        
        ifz=shot%src%ifz; iz=shot%src%iz; ilz=shot%src%ilz
        ifx=shot%src%ifx; ix=shot%src%ix; ilx=shot%src%ilx
        ify=shot%src%ify; iy=shot%src%iy; ily=shot%src%ily
            
        if(if_hicks) then
            select case (shot%src%icomp)
                case (1)
                w = sum( pinvkpa(ifz:ilz,ifx:ilx,ify:ily)*pp(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (2)
                w = sum( pvx(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (3)
                w = sum( pvy(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
                case (4)
                w = sum( pvz(ifz:ilz,ifx:ilx,ify:ily) *shot%src%interp_coeff )
                
            end select
            
        else
            select case (shot%src%icomp)
                case (1) !p[iz,ix,iy]
                w=pinvkpa(iz,ix,iy)*pp(iz,ix,iy)
                
                case (2) !vx[iz,ix-0.5,iy]
                w=pvx(iz,ix,iy)
                
                case (3) !vy[iz,ix,iy-0.5]
                w=pvy(iz,ix,iy)
                
                case (4) !vz[iz-0.5,ix,iy]
                w=pvz(iz,ix,iy)
            end select
            
        endif
        
    end subroutine
    
    !========= for wavefield correlation ===================
    !gkpa = -1/kpa2 dp_dt x adjp
    !dp_dt^it+1 ~= ^it+1.5 - p^it+0.5
    ! adjp^it+1 is not available, use adjp^it+0.5 instead (->very accurate gkpa)
    
    !grho = dvx_dt x adjvx + dvy_dt x adjvy + dvz_dt x adjvz
    !dv_dt^it+0.5 ~= v^it+1 - v^it
    ! adjv^it+0.5 is not available, use adjv^it instead (->grho is reasonably accurate)
    !use v[i+1]-v[i] to approximate v[i+0.5], so is adjv
    
    subroutine field_correlation_stresses(it,sf,rf,sb,rb,corr)
        type(t_field) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(*) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),1)
        ilz=min(sb(2),rb(2),cb%mz)
        ifx=max(sb(3),rb(3),1)
        ilx=min(sb(4),rb(4),cb%mx)
        ify=max(sb(5),rb(5),1)
        ily=min(sb(6),rb(6),cb%my)
        if(.not.m%is_cubic) then
            ify=1; ily=1
        endif
        
        
        dsp=0.
         rp=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,dsp,rp)
        !$omp do schedule(dynamic)
        do iy=ify,ily
        do ix=ifx,ilx
        do iz=ifz,ilz
            i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
            j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !corr has no boundary layers
            
            dsp=sf%prev_p(i)-sf%p(i)
             rp=rf%p(i)*2.    !rp=rf%prev_p(i)+rf%p(i)
            
            corr(j)=corr(j) + 0.5*dsp*rp
            
        end do
        end do
        end do
        !$omp end do 
        !$omp end parallel
        
    end subroutine
    
    subroutine field_correlation_velocities(it,sf,rf,sb,rb,corr)
        type(t_field) :: sf,rf
        integer,dimension(6) :: sb,rb
        real,dimension(*) :: corr
        
        !nonzero only when sf touches rf
        ifz=max(sb(1),rb(1),1)
        ilz=min(sb(2),rb(2),cb%mz-1)
        ifx=max(sb(3),rb(3),1)
        ilx=min(sb(4),rb(4),cb%mx-1)
        ify=max(sb(5),rb(5),1)
        ily=min(sb(6),rb(6),cb%my-1)
        if(.not.m%is_cubic) then
            ify=1; ily=1
        endif
        
        
        dsvx=0.; dsvy=0.; dsvz=0.
         rvx=0.;  rvy=0.;  rvz=0.
        
        !$omp parallel default (shared)&
        !$omp private(iz,ix,iy,i,j,&
        !$omp         iz_ix_iy,izp1_ix_iy,&
        !$omp         iz_ixp1_iy,&
        !$omp         iz_ix_iyp1,&
        !$omp         dsvx,dsvy,dsvz,&
        !$omp          rvx, rvy, rvz)
        !$omp do schedule(dynamic)
        do iy=ify,ily
        do ix=ifx,ilx
        do iz=ifz,ilz
            i=(iz-cb%ifz)+(ix-cb%ifx)*cb%nz+(iy-cb%ify)*cb%nz*cb%nx+1 !field has boundary layers
            j=(iz-1)     +(ix-1)     *cb%mz+(iy-1)     *cb%mz*cb%mx+1 !corr has no boundary layers
            
            iz_ix_iy  =i          !iz,ix,iy
            izp1_ix_iy=i+1        !iz+1,ix,iy
            iz_ixp1_iy=i  +cb%nz     !iz,ix+1,iy
            iz_ix_iyp1=i  +cb%nz*cb%nx  !iz,ix,iy+1
            
            dsvx=(sf%prev_vx(iz_ixp1_iy)+sf%prev_vx(iz_ix_iy)) - (sf%vx(iz_ixp1_iy)+sf%vx(iz_ix_iy))
            if(m%is_cubic) &
            dsvy=(sf%prev_vy(iz_ix_iyp1)+sf%prev_vy(iz_ix_iy)) - (sf%vy(iz_ix_iyp1)+sf%vy(iz_ix_iy))
            dsvz=(sf%prev_vz(izp1_ix_iy)+sf%prev_vz(iz_ix_iy)) - (sf%vz(izp1_ix_iy)+sf%vz(iz_ix_iy))

            rvx=rf%vx(iz_ixp1_iy)+rf%vx(iz_ix_iy)
            if(m%is_cubic) &
            rvy=rf%vy(iz_ix_iyp1)+rf%vy(iz_ix_iy)
            rvz=rf%vz(izp1_ix_iy)+rf%vz(iz_ix_iy)
            
            corr(j)=corr(j) + 0.25*( dsvx*rvx + dsvy*rvy + dsvz*rvz )
        
        enddo
        enddo
        enddo
        !$omp end do
        !$omp end parallel
        
    end subroutine
    
    subroutine field_correlation_scaling(grad)
        real,dimension(cb%mz,cb%mx,cb%my,2) :: grad
        
        grad(:,:,:,1)=grad(:,:,:,1) * (-pinvkpa(1:cb%mz,1:cb%mx,1:cb%my)*pinvkpa(1:cb%mz,1:cb%mx,1:cb%my))
        
        !set unit of gkpa to be [m3], grho to be [m5/s2]
        !such that after multiplied by (kpa_max-kpa_min) or (rho_max-rho_min) (will be done in m_parameterization)
        !the unit of parameter update is [Nm], same as Lagrangian
        !and the unit of gradient scaling factor is [1/N/m] (in m_scaling)
        !(parameters are unitless)
        grad=grad*m%cell_size
        
!        if(m%rho(1,1,1)>10.) then !density is in [kg/m3]
!            grad=grad*1e6
!            call hud('Gradient has been x1e6 to avoid potential fp underflow')
!            call hud('(If rho is in [kg/m3], gkpa magnitude can be ~1e-36)')
!        endif
        
    end subroutine
    
end

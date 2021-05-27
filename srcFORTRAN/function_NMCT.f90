module param_mod
    implicit none
    integer,        parameter :: rkind            = 8
    integer,        parameter :: MAX_FILENAME_LEN = 200

    real(rkind),    parameter :: pi = 4.0_rkind * atan(1.0_rkind)   ! 双精度pi值
    complex(rkind), parameter :: ci = cmplx(0.0, 1.0, kind=rkind)   ! 复数
end module

module cheb_mod
    use param_mod
    implicit none

    interface Cheb
        module procedure ChebReal
        module procedure ChebComplex
    end interface

    interface Convolution
        module procedure ConvolutionReal
        module procedure ConvolutionComplex
    end interface

contains

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Chebyshev正变换。已知原函数在CGL节点上的取值，求其Chebyshev变换后的展开系数.
    !> 给定原函数在n个Chebyshev-Gauss-Lobbatto节点z(1:n)上的取值fx(1:n)，
    !> 求n阶Chebyshev变换后的展开系数fk(0:n)。
    !
    !> @param[in]  z(1:n)   n个Chebyshev-Gauss-Lobbatto节点
    !> @param[in] fx(1:n)   原函数在给定 z(1:n) 处的值
    !> @param[out]  fk(0:n) 各阶 Chebyshev 多项式对应的展开系数
    !> @return none
    !---------------------------------------------------------------------------  
    subroutine ChebReal(fk, fx, z)
        implicit none
        real(rkind),	intent(in)    :: z(:)
        real(rkind),    intent(inout) :: fk(:), fx(:)
        real(rkind),	allocatable   :: T(:,:)
        integer                       :: m, i
        m = size(fx,1)
        allocate(T(m, m))             

        T = 0.0_rkind
        fk= 0.0_rkind
        
        do i=0, m-1
            T(:,i+1) = cos(i * acos(z))	  
        enddo   

        fx(1) = fx(1) * 0.5_rkind
        fx(m) = fx(m) * 0.5_rkind

        fk = matmul( transpose(T), fx ) * 2.0_rkind / (m - 1)
        
        fk(1) = fk(1) * 0.5_rkind
        fk(m) = fk(m) * 0.5_rkind
        
    end subroutine ChebReal

    subroutine ChebComplex(fk, fx, z)
        implicit none
        real(rkind),    intent(in)    :: z(:)
        complex(rkind), intent(inout) :: fk(:), fx(:)
        real(rkind),    allocatable   :: T(:, :)
        integer                       :: m, i

        m = size(fx, 1)
        allocate ( T(m, m) )

        T  = 0.0_rkind
        fk = 0.0_rkind

        do i = 0, m - 1
            T(:, i + 1) = cos( i * acos(z) )
        end do

        fx(1) = fx(1) * 0.5_rkind
        fx(m) = fx(m) * 0.5_rkind

        fk = matmul( transpose(T), fx ) * 2.0_rkind / (m - 1)

        fk(1) = fk(1) * 0.5_rkind
        fk(m) = fk(m) * 0.5_rkind

    end subroutine ChebComplex

   !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Chebyshev反变换。已知一组函数的Chebyshev变换系数，求取这组函数在指定位置处的函数近似值.
    !> 给定d个函数关于m阶Chebyshev变换的展开系数fk(0:m, 1:d)，以及给定的自变量离散点z(1:n)，
    !> 求取原函数在这些离散点上的函数近似值fx(1:n, 1:d)。
    !
    !> @param[in]  z(1:n)   表示原函数要在这些离散点上求取近似值
    !> @param[in]  fk(0:m, 1:d)   fk(:, j) 表示各阶 Chebyshev 多项式对应的展开系数（列向量）
    !> @param[out] fx(1:n, 1:d)   fx(:, j) 表示原函数在给定 z(1:n) 处的近似值
    !> @return none
    !---------------------------------------------------------------------------  
    subroutine InvChebMatrix(fx, fk, z)
        implicit none
        real(rkind),    intent(in)  :: z(:)
        complex(rkind), intent(in)  :: fk(:, :)
        complex(rkind), intent(out) :: fx(:, :)
        integer                     :: i
        real(rkind)                 :: T(size(z), size(fk, 1))

        T  = 0.0_rkind
        fx = 0.0_rkind

        do i = 0, size(fk, 1) - 1
            T(:, i+1) = cos( i * acos( z(:) ) )
        end do

        fx = matmul(T, fk)
    end subroutine InvChebMatrix

    subroutine DerivationMatrix(D, m)
        implicit none
        integer,      intent(in) :: m
        real(rkind),  intent(out):: D(m, m)
        integer                  :: i, j
    
        D = 0.0_rkind
        do i = 1, m
        do j = i + 1, m
            if ( mod((i+j), 2) == 1 ) D(i, j) = 2.0_rkind * j - 2.0_rkind
        end do
        end do
        D(1, :) = D(1, :) * 0.5_rkind
    
    end subroutine DerivationMatrix
    
    subroutine ConvolutionReal(Co, v)
        implicit none
        real(rkind), intent(in)   :: v  (:)
        real(rkind), intent(out)  :: Co (:, :)
        integer                   :: i, j, k, n
    
        Co = 0.0_rkind
    
		n  = size(v)
		do  i = 1, n
		do  k = 1, n
		     j= k-i+1
		     if (1 <= j .and. j <= n) then
			    Co(k,i) = Co(k,i) + v(j) 
		     endif
		     j = i-k+1
		     if (j <=n .and. j >= 1) then
				Co(k,i) = Co(k,i) + v(j) 
		     endif
		     if (k >1) then 
				j = i+k-1   
			    if (1 <= j .and. j <= n) then
					Co(k,i) = Co(k,i) + v(j)
				endif 
		     endif   
		     Co(k,i) = Co(k,i) * 0.5_rkind
		enddo
		enddo
    
    end subroutine ConvolutionReal
    
    subroutine ConvolutionComplex(Co,v)
        implicit none
        complex(rkind), intent(in)   ::v(:)
        complex(rkind), intent(out)  ::Co(:,:)
        integer                      ::i, j, k, n  
      
        Co = 0.0_rkind + ci * 0.0_rkind
    
		n  = size(v)
		do  i = 1, n
		do  k = 1, n
		     j= k-i+1
		     if (1 <= j .and. j <= n) then
			    Co(k,i) = Co(k,i) + v(j) 
		     endif
		     j = i-k+1
		     if (j <=n .and. j >= 1) then
				Co(k,i) = Co(k,i) + v(j) 
		     endif
		     if (k >1) then 
				j = i+k-1   
			    if (1 <= j .and. j <= n) then
					Co(k,i) = Co(k,i) + v(j)
				endif 
		     endif   
		     Co(k,i) = Co(k,i) * 0.5_rkind
		enddo
		enddo
    
    end subroutine ConvolutionComplex

end module cheb_mod

module util_mod
    use param_mod
    implicit none

contains

    subroutine assert (cond, msg)
        implicit none
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: msg

        if (.not. cond) then
            write(*, *) "ERROR : ", msg
            stop
        end if
    end subroutine assert
end module util_mod

module nmct_mod
    use param_mod
    use cheb_mod
    implicit none

contains

    subroutine ReadEnvParameter(casename, Nw, Nb, cpmax, freq, zs, zr, rmax, dr, hinterface, Hbottom,&
        dz, Lowerboundary, tlmin, tlmax, rhow, rhob, alphaw, alphab, cw, cb, data_file)
        use util_mod
        implicit none
        character(len=MAX_FILENAME_LEN), intent(out) :: casename
		character(len=MAX_FILENAME_LEN), intent(in)  :: data_file
        integer,                         intent(out) :: Nw
        integer,                         intent(out) :: Nb
        integer,                         intent(out) :: Lowerboundary
        real(rkind),                     intent(out) :: cpmax
        real(rkind),                     intent(out) :: dr
        real(rkind),                     intent(out) :: zs
        real(rkind),                     intent(out) :: zr
        real(rkind),                     intent(out) :: rmax
        real(rkind),                     intent(out) :: freq
        real(rkind),                     intent(out) :: hinterface
        real(rkind),                     intent(out) :: Hbottom
        real(rkind),                     intent(out) :: dz
        real(rkind),                     intent(out) :: tlmin
        real(rkind),                     intent(out) :: tlmax
        real(rkind), allocatable,        intent(out) :: rhow(:), alphaw(:), cw(:)
        real(rkind), allocatable,        intent(out) :: rhob(:), alphab(:), cb(:)
        real(rkind), allocatable, dimension(:)       :: temp_alphaw, depw, temp_rhow, temp_cw
        real(rkind), allocatable, dimension(:)       :: temp_alphab, depb, temp_rhob, temp_cb
        integer										 :: n_w,n_b,i
                                        
        open(unit=1, status='unknown', file=data_file) 
    
        read (1, *) casename
        read (1, *) Nw
        read (1, *) Nb
        read (1, *) cpmax
        read (1, *) freq
        read (1, *) zs
        read (1, *) zr
        read (1, *) rmax
        read (1, *) dr
        read (1, *) hinterface
        read (1, *) Hbottom
        read (1, *) dz
        read (1, *) Lowerboundary
        read (1, *) tlmin
        read (1, *) tlmax
        read (1, *) n_w
        read (1, *) n_b
        
        !read the param_mod of ocean
        allocate(depw(n_w), temp_cw(n_w), temp_rhow(n_w), temp_alphaw(n_w))
        allocate(depb(n_b), temp_cb(n_b), temp_rhob(n_b), temp_alphab(n_b))

        if (hinterface > 0.0 .and. Hbottom > hinterface) then
            do i = 1, n_w
                read(1, *)depw(i), temp_cw(i), temp_rhow(i), temp_alphaw(i)
            end do

            do i = 1, n_b
                read(1, *)depb(i), temp_cb(i), temp_rhob(i), temp_alphab(i)
            end do
        else
            call assert(.false., "interface must less than Hbottom and greater than 0!")
        end if
        close(1)

        call assert( Nw > 2 .and. Nb > 2, 'Nw and Nb must greater than 2! ' )
        call assert( depw(1) == 0.0 .and. depw(n_w) == hinterface .and. &
                     depb(1) == hinterface .and. depb(n_b) == Hbottom,  &
                     'input sound profile is unsuitable!' )
        call assert( hinterface / dz == floor(hinterface / dz) .and. &
                     Hbottom / dz == floor(Hbottom / dz), &
                     'The input dz unsuitable!' )
					 
        call assert( zs > 0 .and. zs < Hbottom .and. zr > 0 .and. zr < Hbottom,&
											'Nw and Nb must greater than 2! ' )
		
        call assert( rmax / dr == floor(rmax / dr), 'Please reinput the dr and rmax!' )

        call assert( Lowerboundary == 0 .or. Lowerboundary == 1, &
                     'The lower boundary must be rigid (0) or soft (1)!' )

        !Interplating to the CGL points
        allocate(cw(Nw+1), cb(Nb+1), rhow(Nw+1), rhob(Nb+1), alphaw(Nw+1), alphab(Nb+1))

        call Interpolation(depw,temp_cw,cw,temp_rhow,rhow,temp_alphaw,alphaw,Nw)	
        call Interpolation(depb,temp_cb,cb,temp_rhob,rhob,temp_alphab,alphab,Nb)	

        deallocate(depw, depb, temp_cw, temp_cb, temp_rhow, temp_rhob, temp_alphaw, temp_alphab)

    end subroutine ReadEnvParameter

    subroutine Interpolation(dep,b1,b2,c1,c2,d1,d2,N)
        implicit none
        integer,       intent(in)   ::N
        real(rkind),   intent(in)   ::dep(:), b1(:), c1(:), d1(:)
        real(rkind),   intent(out)  ::b2(:), c2(:), d2(:)
        real(rkind)                 ::x(N+1), z(N+1)
        integer                     ::i, j, m
        
        m = size(dep)
		
        do i = 1, N + 1
            x(i) = cos((i - 1) * pi / N)
            z(i) = ((dep(m) + dep(1)) / (dep(m) - dep(1)) - x(i)) * (dep(m) - dep(1)) / 2.0
        end do

        do i=1,N+1	
            do j=1,m-1
                if((z(i) >= dep(j)).and.(z(i) <= dep(j+1))) then
                    b2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * b1(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * b1(j)
                    
                    c2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * c1(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * c1(j)
                    
                    d2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * d1(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * d1(j)				
                    
                endif				
            enddo 
            
        if(z(i) > dep(m)) then
            b2(i) = b1(m)
            c2(i) = c1(m)
            d2(i) = d1(m)							
        end if
            
        if(z(i) < dep(1)) then
            b2(i) = b1(1)
            c2(i) = c1(1)	
            d2(i) = d1(1)					
            endif
        enddo	
        
    end subroutine Interpolation

    subroutine Initialization(Nw, Nb, freq, rmax, dr, zs, rhow, rhob, cw, cb, alphaw, alphab, &
        hinterface, Hbottom, nr, r, rhozs, kw, kb)
        implicit none
        integer,                     intent(in)  :: Nw, Nb
        integer,                     intent(out) :: nr
        real(rkind),                 intent(in)  :: freq
        real(rkind),                 intent(in)  :: rmax
        real(rkind),                 intent(in)  :: dr
        real(rkind),                 intent(in)  :: zs
        real(rkind),                 intent(in)  :: hinterface
        real(rkind),                 intent(in)  :: Hbottom
        real(rkind), dimension(Nw+1),intent(in)  :: rhow, cw, alphaw
        real(rkind), dimension(Nb+1),intent(in)  :: rhob, cb, alphab
        real(rkind),                 intent(out) :: rhozs
        real(rkind),    allocatable, intent(out) :: r(:)
        complex(rkind), allocatable, intent(out) :: kw(:), kb(:)

        real(rkind), dimension(Nw+1)             :: x1, z1
        real(rkind), dimension(Nb+1)             :: x2, z2
        real(rkind)                              :: w
        integer                                  :: i

        w  = 2.0_rkind * pi * freq
        nr = int(rmax / dr)

        allocate ( r(nr) )
        do i = 1, nr
            r(i) = i * dr
        end do
        do i = 1, Nw+1
            x1(i) = cos( (i-1) * pi / Nw )
        end do
        do i = 1, Nb+1
            x2(i) = cos( (i-1) * pi / Nb )
        end do
        
        if ( zs <= hinterface ) then
            z1 = (1.0_rkind - x1) * 0.5_rkind * hinterface;
            do i = 1, Nw
                if ( zs >= z1(i) .and. zs <= z1(i+1) ) then
                    rhozs = (zs - z1(i)) / (z1(i+1) - z1(i)) * rhow(i+1) &
                        + (z1(i+1) - zs) / (z1(i+1) - z1(i)) * rhow(i)
                end if
            end do
        else
            z2 = (1.0_rkind - x2) * 0.5_rkind * (Hbottom - hinterface) + hinterface
            do i = 1, Nb
                if (zs >= z2(i) .and. zs <= z2(i+1)) then
                    rhozs = (zs - z2(i)) / (z2(i+1) - z2(i)) * rhob(i+1) &
                        + (z2(i+1) - zs) / (z2(i+1) - z2(i)) * rhob(i)
                end if
            end do
        end if

        allocate (kw(Nw+1), kb(Nb+1))
        kw = w / cw * (1.0_rkind + ci * alphaw / (40.0_rkind * pi * log10(exp(1.0_rkind))))
        kb = w / cb * (1.0_rkind + ci * alphab / (40.0_rkind * pi * log10(exp(1.0_rkind))))
        
    end subroutine Initialization

    subroutine EigenValueVector(Nw, Nb, hinterface, Hbottom, kw, kb, rhow, rhob, &
        Lowerboundary, kr, eigvectorw, eigvectorb)
        implicit none
        integer,                     intent(in)    :: Nw, Nb
        integer,                     intent(in)    :: Lowerboundary
        real(rkind),                 intent(in)    :: hinterface, Hbottom
        real(rkind),                 intent(in)    :: rhow(Nw+1), rhob(Nb+1)
        complex(rkind),              intent(inout) :: kw(Nw+1), kb(Nb+1)
        complex(rkind), allocatable, intent(out)   :: kr(:)
        complex(rkind), allocatable, intent(out)   :: eigvectorw(:, :), eigvectorb(:, :)
        real(rkind),    dimension(Nw+1)            :: rhow1, rhow2, Pu
        real(rkind),    dimension(Nw+1, Nw+1)      :: D1
        real(rkind),    dimension(Nb+1)            :: rhob1, rhob2, Pd
        real(rkind)                                :: D2(Nb+1, Nb+1)
        real(rkind)                                :: x1(Nw+1), tem11(Nw+1), tem21(Nw+1, Nw+1)
        real(rkind)                                :: x2(Nb+1), tem12(Nb+1), tem22(Nb+1, Nb+1)
        complex(rkind)                             :: A(Nw+1, Nw+1), tem31(Nw+1), tem41(Nw+1, Nw+1)
        complex(rkind)                             :: B(Nb+1, Nb+1), tem32(Nb+1), tem42(Nb+1, Nb+1)
        complex(rkind)                             :: U(Nw+Nb+2, Nw+Nb+2)
        complex(rkind)                             :: bottom_boundary(Nb+1)
        complex(rkind)                             :: L11(Nw+Nb-2, Nw+Nb-2)
        complex(rkind)                             :: L12(Nw+Nb-2, 4)
        complex(rkind)                             :: L21(4, Nw+Nb-2)
        complex(rkind)                             :: L22(4, 4)
        complex(rkind)                             :: L(Nw+Nb-2, Nw+Nb-2)
        complex(rkind)                             :: v2(4, Nw+Nb-2)
        complex(rkind)                             :: VL(Nw+Nb-2)
        complex(rkind)                             :: VR(Nw+Nb-2, Nw+Nb-2)
        complex(rkind)                             :: WORK(2*(Nw+Nb-2))
        real(rkind)                                :: RWORK(2*(Nw+Nb-2))
        integer                                    :: i, info, j(1), IPIV(4)

        allocate(kr(Nw+Nb-2), eigvectorw(Nw+1, Nw+Nb-2), eigvectorb(Nb+1, Nw+Nb-2))

        call DerivationMatrix(D1, Nw+1)
        call DerivationMatrix(D2, Nb+1)

        do i = 1, Nw + 1
            x1(i) = cos( (i-1) * pi / Nw )
            kw(i) = kw(i) ** 2
        end do
        do i = 1, Nb + 1
            x2(i) = cos((i - 1) * pi / Nb)
            kb(i) = kb(i) ** 2
        end do

        rhow2 = rhow
        rhob2 = rhob
        rhow1 = 1.0_rkind / rhow
        rhob1 = 1.0_rkind / rhob

        call Cheb(tem11, rhow2, x1)
        call Convolution(tem21, tem11)
        A = 4.0_rkind / hinterface / hinterface * matmul(tem21, D1)
        call Cheb(tem11, rhow1, x1)
        call Convolution(tem21, tem11)
        A = matmul(A, tem21)
        A = matmul(A, D1)
        call Cheb(tem31, kw, x1)
        call Convolution(tem41, tem31)
        A = A + tem41

        call Cheb(tem12, rhob2, x2)
        call Convolution(tem22, tem12)
        B = 4.0_rkind / (Hbottom - hinterface) / (Hbottom - hinterface) * matmul(tem22, D2)
        call Cheb(tem12, rhob1, x2)
        call Convolution(tem22, tem12)
        B = matmul(B, tem22)
        B = matmul(B, D2)
        call Cheb(tem32, kb, x2)
        call Convolution(tem42, tem32)
        B = B + tem42

        U = 0.0_rkind
        U(1:Nw-1, 1:Nw-1) = A(1:Nw-1, 1:Nw-1)
        U(1:Nw-1, Nw+Nb-1:Nw+Nb) = A(1:Nw-1, Nw:Nw+1)
        U(Nw:Nw+Nb-2, Nw:Nw+Nb-2) = B(1:Nb-1, 1:Nb-1)
        U(Nw:Nw+Nb-2, Nw+Nb+1:Nw+Nb+2) = B(1:Nb-1, Nb:Nb+1)
        !upper boundary
        U(Nw+Nb-1, 1:Nw-1)        = 1.0_rkind
        U(Nw+Nb-1, Nw+Nb-1:Nw+Nb) = 1.0_rkind
        !lower boundary
        do i = 1, Nb + 1
            bottom_boundary(i) = (-1.0_rkind) ** (i-1)
        end do
        if (Lowerboundary == 1)  bottom_boundary = matmul(bottom_boundary, D2)
        U(Nw+Nb+2, Nw:Nw+Nb-2)      = bottom_boundary(1 : Nb-1)
        U(Nw+Nb+2, Nw+Nb+1:Nw+Nb+2) = bottom_boundary(Nb: Nb+1)
        !interface equal
        do i = 1, Nw - 1
        U(Nw+Nb, i) = (-1.0_rkind) ** (i-1)
        end do
        U(Nw+Nb, Nw+Nb-1)        = (-1.0_rkind) ** (Nw - 1)
        U(Nw+Nb, Nw+Nb)          = (-1.0_rkind) ** Nw
        U(Nw+Nb, Nw:Nw+Nb-2)     = -1.0_rkind
        U(Nw+Nb, Nw+Nb+1:Nw+Nb+2)= -1.0_rkind
        !interface derivative
        do i = 1, Nw + 1
        Pu(i) = (-1.0_rkind) ** (i-1)
        end do
        Pd = 1.0_rkind
        Pu = 1.0_rkind / rhow(Nw+1) / hinterface * matmul(Pu, D1)
        Pd = -1.0_rkind / rhob(1) / (Hbottom - hinterface) * matmul(Pd, D2)

        U(Nw+Nb+1, 1:Nw-1)         = Pu(1:Nw-1)
        U(Nw+Nb+1, Nw:Nw+Nb-2)     = Pd(1:Nb-1)
        U(Nw+Nb+1, Nw+Nb-1:Nw+Nb)  = Pu(Nw:Nw+1)
        U(Nw+Nb+1, Nw+Nb+1:Nw+Nb+2)= Pd(Nb:Nb+1)

        L11 = U(1:Nw+Nb-2, 1:Nw+Nb-2)
        L12 = U(1:Nw+Nb-2, Nw+Nb-1:Nw+Nb+2)
        L21 = U(Nw+Nb-1:Nw+Nb+2, 1:Nw+Nb-2)
        L22 = U(Nw+Nb-1:Nw+Nb+2, Nw+Nb-1:Nw+Nb+2)

		call zgesv(4,Nw+Nb-2,L22,4,IPIV,L21,4,info)
		L = L11 - matmul(L12, L21)

        call zgeev('N', 'V', Nw+Nb-2, L, Nw+Nb-2, kr, VL, 1, VR, Nw+Nb-2, WORK, 2*(Nw+Nb-2), RWORK, INFO)

		v2 = -matmul(L21, VR)

        eigvectorw = 0.0_rkind
        eigvectorb = 0.0_rkind
        eigvectorw(1:Nw-1,  :) = VR(1:Nw-1,     :)
        eigvectorw(Nw:Nw+1, :) = v2(1:2,        :)
        eigvectorb(1:Nb-1,  :) = VR(Nw:Nw+Nb-2, :)
        eigvectorb(Nb:Nb+1, :) = v2(3:4,        :)

        kr = sqrt(kr)
        !L and L11 store the sorted eigenvectors respectively. VL stores the sorted eigenvalpsi_zs
        do i = 1, Nw + Nb - 2
            j = maxloc(real(kr))
            VL(i) = kr(j(1))
            L  (1:Nw+1, i) = eigvectorw(:, j(1))
            L11(1:Nb+1, i) = eigvectorb(:, j(1))
            kr(j(1)) = - 1.0_rkind
        end do

        kr = VL
        eigvectorw = L  (1:Nw+1, :)
        eigvectorb = L11(1:Nb+1, :)

    end subroutine EigenValueVector

    subroutine NumofModes(freq, kr, nmodes, cpmax)
        use util_mod
        implicit none
        integer,        intent(out):: nmodes
        real(rkind),    intent(in) :: freq
        real(rkind),    intent(in) :: cpmax
        complex(rkind), intent(in) :: kr(:)
        real(rkind), allocatable   :: cp(:)
        integer                    :: i

        allocate ( cp(size(kr)) )
        cp = 2.0_rkind * pi * freq / real(kr, kind=rkind)
        nmodes = 0
        do i = 1, size(kr)
            if (cp(i) <= cpmax) nmodes = i
        end do

        call assert( nmodes /= 0, 'Incorrect maximum phase speed input!')

        deallocate(cp)
    end subroutine NumofModes

    subroutine GenerateModes(nmodes, dz, hinterface, Hbottom, eigvectorw, eigvectorb, psi, z)
        implicit none
        integer,                         intent(in)  :: nmodes
        real(rkind),                     intent(in)  :: dz, hinterface, Hbottom
        complex(rkind), dimension(:, :), intent(in)  :: eigvectorw, eigvectorb
        real(rkind),    allocatable,     intent(out) :: z(:)
        complex(rkind), allocatable,     intent(out) :: psi(:, :)

        real(rkind),    allocatable, dimension(:)    :: zt1, zt2
        complex(rkind), allocatable, dimension(:, :) :: psi1, psi2
        integer                                      :: i

        allocate ( zt1(nint(hinterface / dz) + 1) )
        allocate ( zt2(nint((Hbottom - hinterface) / dz) + 1) )
        allocate ( z  (nint(Hbottom / dz) + 1) )

        do i = 1, nint(hinterface / dz) + 1
            zt1(i) = (i-1) * dz
        end do
        do i = 1, nint((Hbottom - hinterface) / dz) + 1
            zt2(i) = hinterface + (i - 1) * dz
        end do
        do i = 1, nint(Hbottom / dz) + 1
            z(i) = (i - 1) * dz
        end do
        zt1 = -2.0_rkind / hinterface * zt1 + 1.0_rkind
        zt2 = -2.0_rkind / (Hbottom - hinterface) * zt2 + (Hbottom + hinterface) / (Hbottom - hinterface)

        allocate ( psi1(size(zt1), nmodes), psi2(size(zt2), nmodes), psi(size(z), nmodes) )
        call InvChebMatrix(psi1, eigvectorw(:, 1:nmodes), zt1)
        call InvChebMatrix(psi2, eigvectorb(:, 1:nmodes), zt2)

        psi(1:size(zt1)-1, :)	  = psi1(1:size(zt1)-1, :)
        psi(size(zt1):size(z), :) = psi2

        deallocate(psi1, psi2, zt1, zt2)
    end subroutine GenerateModes

    subroutine Normalization(eigvectorw, eigvectorb, rhow, rhob, hinterface, Hbottom, Nw, Nb, nmodes, psi)
        implicit none
        complex(rkind), dimension(:, :), intent(in)  :: eigvectorw, eigvectorb
        complex(rkind), dimension(:, :), intent(out) :: psi

        integer            :: i, Nw, Nb, nmodes
        real(rkind)        :: hinterface, Hbottom, rhow(Nw), rhob(Nb), x1(Nw-1), x2(Nb-1), &   
                            rw(Nw), rb(Nb), Co11(Nw, Nw), Co22(Nb, Nb), P(Nw), Q(Nb)
        complex(rkind)     :: Co1(Nw, Nw), Co2(Nb, Nb), f1(Nw), f2(Nb), norm(nmodes)

        do i = 1, Nw
            x1(i) = cos( (i-1) * pi / (Nw-1) )
        end do
        do i = 1, Nb
            x2(i) = cos( (i-1) * pi / (Nb-1) )
        end do

        rhow = 1.0_rkind / rhow
        rw = 0.0_rkind
        call Cheb(rw, rhow, x1) ! rhow is modified
        call Convolution(Co11, rw)
        rhob = 1.0_rkind / rhob
        rb = 0.0_rkind
        call Cheb(rb, rhob, x2)
        call Convolution(Co22, rb)

        P = 0.0_rkind
        Q = 0.0_rkind
        do i = 0, Nw-1, 2
            P(i+1) = - 2.0_rkind / (i * i - 1.0_rkind)
        end do
        do i = 0, Nb-1, 2
            Q(i+1) = - 2.0_rkind / (i * i - 1.0_rkind)
        end do

        do i = 1, nmodes
            call Convolution(Co1, eigvectorw(:, i))
            f1 = matmul(Co11, matmul(Co1, eigvectorw(:, i)))

            call Convolution(Co2, eigvectorb(:, i))
            f2 = matmul(Co22, matmul(Co2, eigvectorb(:, i)))

            norm(i) = sqrt(dot_product(P, f1) * hinterface * 0.5_rkind + &
                        dot_product(Q, f2) * (Hbottom - hinterface) * 0.5_rkind)
        end do

        do i = 1, nmodes
            psi(:, i) = psi(:, i) / norm(i)
        end do

    end subroutine Normalization

    subroutine SynthesizeSoundField(nmodes, nr, r, kr, rhozs, zs, dz, psi, tl)
        implicit none
        integer,                  intent(in)    :: nmodes, nr
        real(rkind),              intent(in)    :: r(nr), rhozs, zs, dz
        real(rkind), allocatable, intent(inout) :: tl(:, :)
        complex(rkind),           intent(in)    :: kr(:)
        complex(rkind),           intent(inout) :: psi(:, :)

        complex(rkind), allocatable             :: p(:, :)
        complex(rkind)                          :: bessel(nmodes, nr)
        complex(rkind)                          :: psizs(nmodes, nmodes)
        integer                                 :: i, k, s, IERR1, IERR2
        real(rkind)                             :: CYR, CYI

        do k = 1, nmodes
        do i = 1, nr
            bessel(k, i) = r(i) * kr(k)
            call ZBESH(real(bessel(k, i)), aimag(bessel(k, i)), 0.0_rkind, 1, 1, 1, CYR, CYI, IERR1, IERR2)
            bessel(k, i) = cmplx(CYR, CYI)
        end do
        end do
            
        allocate(p(size(psi, 1), nr), tl(size(psi, 1), nr))
        s = ceiling(zs / dz)
        psizs(:, :)=0.0_rkind
        do k = 1, nmodes
            psizs(k, k) = (zs / dz - s) * psi(s + 1, k) + ( s + 1 - zs / dz) * psi(s, k)
        end do
        
        psi=matmul(psi, psizs)
        p = matmul(psi, bessel)
        p = p * ci * pi / rhozs
        
        tl = - 20.0_rkind * log10(abs(p))
    end subroutine SynthesizeSoundField

    subroutine SaveSoundField(filename, tlmin, tlmax, r, z, tl)
        implicit none
        character(len=MAX_FILENAME_LEN), intent(in) :: filename
        real(rkind),                     intent(in) :: tlmin, tlmax
        real(rkind), dimension(:),       intent(in) :: r, z
        real(rkind), dimension(:, :),    intent(in) :: tl

        open(unit=20, status='unknown', file=filename, access='stream', form='unformatted')
        write(20)  size(z), size(r), tlmin, tlmax, z, r, tl
        close(20)
    end subroutine SaveSoundField

end module nmct_mod
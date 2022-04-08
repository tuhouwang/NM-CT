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

    function ChebReal(fx,z)

        implicit none
        real(rkind),    intent(in)  :: z (:)
        real(rkind),    intent(in)  :: fx(:)   
        real(rkind)                 :: ChebReal(size(z))
        real(rkind)                 :: fx2(size(z))
        real(rkind)                 :: T(size(z),size(z)) 
        integer                     :: m, i
        
        m = size(fx)           
        T = 0.0_rkind
        ChebReal = 0.0_rkind

        do i=0, m-1
            T(:,i+1) = cos(i * acos(z))      
        enddo 
        
        fx2    = fx
        fx2(1) = fx2(1) * 0.5_rkind
        fx2(m) = fx2(m) * 0.5_rkind

        ChebReal    = matmul(transpose(T),fx2)*2.0_rkind/(m-1)       
        ChebReal(1) = ChebReal(1) * 0.5_rkind
        ChebReal(m) = ChebReal(m) * 0.5_rkind
        
    end function

    function ChebComplex(fx,z)
    
        implicit none
        real(rkind),    intent(in)  :: z (:)
        complex(rkind), intent(in)  :: fx(:)   
        complex(rkind)              :: ChebComplex(size(z))
        complex(rkind)              :: fx2(size(z))
        real(rkind)                 :: T(size(z),size(z)) 
        integer                     :: m, i
        
        m = size(fx)           
        T = 0.0_rkind
        ChebComplex = 0.0_rkind

        do i=0, m-1
            T(:,i+1) = cos(i * acos(z))      
        enddo 
        
        fx2    = fx
        fx2(1) = fx2(1) * 0.5_rkind
        fx2(m) = fx2(m) * 0.5_rkind

        ChebComplex    = matmul(transpose(T),fx2)*2.0_rkind/(m-1)       
        ChebComplex(1) = ChebComplex(1) * 0.5_rkind
        ChebComplex(m) = ChebComplex(m) * 0.5_rkind
        
    end function

    function InvChebMatrix(fk,z)
    
        implicit none
        real(rkind),    intent(in)  :: z(:)
        complex(rkind), intent(in)  :: fk(:, :)
        complex(rkind)              :: InvChebMatrix(size(z),size(fk,2))
        real(rkind)                 :: T(size(z), size(fk,1))
        integer                     :: i

        T  = 0.0_rkind

        do i = 0, size(fk, 1) - 1
            T(:, i+1) = cos( i * acos( z(:) ) )
        end do

        InvChebMatrix = matmul(T, fk)
        
    end function

    function DerivationMatrix(m)
    
        implicit none
        integer,      intent(in) :: m
        real(rkind)              :: DerivationMatrix(m, m)
        integer                  :: i, j
    
        DerivationMatrix = 0.0_rkind
        do i = 1, m
        do j = i + 1, m
            if ( mod((i+j), 2) == 1 ) then
                DerivationMatrix(i, j) = 2.0_rkind * j - 2.0_rkind
            endif
        end do
        end do
        DerivationMatrix(1, :) = DerivationMatrix(1, :) * 0.5_rkind
    
    end function
    
    function ConvolutionReal(v)
    
        implicit none
        real(rkind), intent(in)   :: v(:)
        real(rkind)               :: ConvolutionReal(size(v),size(v))
        integer                   :: i, j, k, n
    
        ConvolutionReal = 0.0_rkind
    
        n  = size(v)
        do  i = 1, n
        do  k = 1, n
             j = k-i+1
             if (1 <= j .and. j <= n) then
                ConvolutionReal(k,i) = ConvolutionReal(k,i) + v(j)
             endif
             j = i-k+1
             if (j <=n .and. j >= 1) then
                ConvolutionReal(k,i) = ConvolutionReal(k,i) + v(j)
             endif
             if (k > 1) then 
                j = i+k-1   
                if (1 <= j .and. j <= n) then
                    ConvolutionReal(k,i) = ConvolutionReal(k,i) + v(j)
                endif 
             endif   
             ConvolutionReal(k,i) = ConvolutionReal(k,i) * 0.5_rkind
        enddo
        enddo
    
    end function
    
    function ConvolutionComplex(v)
    
        implicit none
        complex(rkind), intent(in)   :: v(:)
        complex(rkind)               :: ConvolutionComplex(size(v),size(v))
        integer                      :: i, j, k, n
    
        ConvolutionComplex = 0.0_rkind
    
        n  = size(v)
        do  i = 1, n
        do  k = 1, n
             j = k-i+1
             if (1 <= j .and. j <= n) then
                ConvolutionComplex(k,i) = ConvolutionComplex(k,i) + v(j)
             endif
             j = i-k+1
             if (j <=n .and. j >= 1) then
                ConvolutionComplex(k,i) = ConvolutionComplex(k,i) + v(j)
             endif
             if (k > 1) then 
                j = i+k-1   
                if (1 <= j .and. j <= n) then
                    ConvolutionComplex(k,i) = ConvolutionComplex(k,i) + v(j)
                endif 
             endif   
             ConvolutionComplex(k,i) = ConvolutionComplex(k,i) * 0.5_rkind
        enddo
        enddo
    
    end function

end module

module util_mod
    
    implicit none

contains

    subroutine assert (cond,msg)
        implicit none
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: msg

        if (.not. cond) then
            write(*, *) 'ERROR : ', msg
            stop
        end if
    end subroutine
    
end module

module nmct_mod

    use param_mod
    use util_mod
    use cheb_mod
    implicit none

contains

    subroutine ReadEnvParameter(casename,Nw,Nb,cpmax,freq,zs,zr,rmax,dr,hinterface,Hb,dz,&
        Lowerboundary,tlmin,tlmax,rhow,rhob,rhoh,alphaw,alphab,alphah,cw,cb,ch,data_file)
        
        implicit none
        character(len=MAX_FILENAME_LEN), intent(out) :: casename
        character(len=MAX_FILENAME_LEN), intent(in)  :: data_file
        character(len=1),                intent(out) :: Lowerboundary
        integer,                         intent(out) :: Nw
        integer,                         intent(out) :: Nb
        real(rkind),                     intent(out) :: cpmax
        real(rkind),                     intent(out) :: dr
        real(rkind),                     intent(out) :: zs
        real(rkind),                     intent(out) :: zr
        real(rkind),                     intent(out) :: rmax
        real(rkind),                     intent(out) :: freq
        real(rkind),                     intent(out) :: hinterface
        real(rkind),                     intent(out) :: Hb
        real(rkind),                     intent(out) :: dz
        real(rkind),                     intent(out) :: tlmin
        real(rkind),                     intent(out) :: tlmax
        real(rkind),                     intent(out) :: rhoh,    alphah,    ch        
        real(rkind), allocatable,        intent(out) :: rhow(:), alphaw(:), cw(:)
        real(rkind), allocatable,        intent(out) :: rhob(:), alphab(:), cb(:)
        real(rkind), allocatable, dimension(:)       :: temp_alphaw, depw, temp_rhow, temp_cw
        real(rkind), allocatable, dimension(:)       :: temp_alphab, depb, temp_rhob, temp_cb
        integer                                      :: n_w,n_b,i
                                        
        open(unit=1, status='unknown', file=data_file) 
    
        read (1, '(A200)') casename
        read (1, *) Nw
        read (1, *) Nb
        read (1, *) cpmax
        read (1, *) freq
        read (1, *) zs
        read (1, *) zr
        read (1, *) rmax
        read (1, *) dr
        read (1, *) hinterface
        read (1, *) Hb
        read (1, *) dz
        read (1, *) tlmin
        read (1, *) tlmax
        read (1, *) n_w
        read (1, *) n_b

        !read the param_mod of ocean
        allocate(depw(n_w), temp_cw(n_w), temp_rhow(n_w), temp_alphaw(n_w))
        allocate(depb(n_b), temp_cb(n_b), temp_rhob(n_b), temp_alphab(n_b))

        if (hinterface > 0.0 .and. Hb > hinterface) then
            do i = 1, n_w
                read(1, *)depw(i), temp_cw(i), temp_rhow(i), temp_alphaw(i)
            end do

            do i = 1, n_b
                read(1, *)depb(i), temp_cb(i), temp_rhob(i), temp_alphab(i)
            end do
        else
            call assert(.false., 'interface must less than Hb and greater than 0!')
        end if
        
        read (1, *) Lowerboundary
        
        call assert(Lowerboundary == 'V' .or. Lowerboundary == 'R' &
               .or. Lowerboundary == 'A', &
             'Error! The lower boundary must be vaccum, rigid or halfspace!') 
             
        if (Lowerboundary == 'A') then
                read(1, *) ch, rhoh, alphah    
        endif
        
        close(1)

        call assert(Nw > 2 .and. Nb > 2, 'Nw and Nb must greater than 2!')
        
        call assert(depw(1) == 0.0 .and. depw(n_w) == hinterface .and. &
                    depb(1) == hinterface .and. depb(n_b) == Hb, &
                     'Error! input sound profile is unsuitable !')
                     
        call assert(zs > 0 .and. zs < Hb .and. zr > 0 .and. zr < Hb,&
                    'zs and zr must be greater than 0 and less than H!')
        
        call assert(rmax / dr == floor(rmax / dr), 'Please reinput the dr and rmax!')

        !Interplating to the CGL points
        allocate(cw(Nw+1),cb(Nb+1),rhow(Nw+1),rhob(Nb+1),alphaw(Nw+1),alphab(Nb+1))

        call Interpolation(depw,temp_cw,cw,temp_rhow,rhow,temp_alphaw,alphaw,Nw)    
        call Interpolation(depb,temp_cb,cb,temp_rhob,rhob,temp_alphab,alphab,Nb)    

        deallocate(depw,depb,temp_cw,temp_cb,temp_rhow,temp_rhob,temp_alphaw,temp_alphab)

    end subroutine

    subroutine Interpolation(dep,b1,b2,c1,c2,d1,d2,N)
        implicit none
        integer,       intent(in)   :: N
        real(rkind),   intent(in)   :: dep(:)
        real(rkind),   intent(in)   :: b1(:), c1(:), d1(:)
        real(rkind),   intent(out)  :: b2(:), c2(:), d2(:)
        real(rkind)                 :: x(N+1), z(N+1)
        integer                     :: i, j, m
        
        m = size(dep)
        
        do i = 1, N+1
            x(i) = cos((i - 1) * pi / N)
            z(i) = ((dep(m) + dep(1)) / (dep(m) - dep(1)) - x(i)) &
            * (dep(m) - dep(1)) / 2.0
        end do

        do i=1, N+1    
            do j=1, m-1
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
        
    end subroutine

    subroutine Interpolation_SSP(dep,c,dep2,c2)
        implicit none
        real   (rkind),intent(in)  :: dep(:), dep2(:)
        complex(rkind),intent(in)  :: c(:)
        complex(rkind),intent(out) :: c2(:)
        integer                    :: n_b, n_b2
        integer                    :: i, j
        n_b  = size(dep)
        n_b2 = size(dep2)
        
        call assert(n_b  == size(c) ,  'The dimensions of the interpolation array do not match!') 
        call assert(n_b2 == size(c2),  'The dimensions of the interpolation array do not match!') 
            
        do i = 1, n_b2
            do j = 1, n_b-1
                if((dep2(i) >= dep(j)) .and. (dep2(i) <= dep(j+1))) then
                    c2(i) = (dep2(i) - dep(j)) / (dep(j+1) - dep(j)) * c(j+1)&
                            +(dep(j+1)-dep2(i))/ (dep(j+1) - dep(j)) * c(j)
                endif
            enddo 
            
            if(dep2(i) > dep(n_b)) then
                c2(i) = c(n_b)               
            end if
            
            if(dep2(i) < dep(1)) then
                c2(i) = c(1)           
            endif
        enddo    

    end subroutine

    subroutine Initialization(Nw,Nb,freq,rmax,dr,zs,rhow,rhob,cw,cb,ch,alphaw, &
                              alphab,alphah,hinterface,Hb,nr,r,rhozs,kw,kb,kh)
        implicit none
        integer,                     intent(in)  :: Nw, Nb
        integer,                     intent(out) :: nr
        real(rkind),                 intent(in)  :: freq
        real(rkind),                 intent(in)  :: rmax
        real(rkind),                 intent(in)  :: dr
        real(rkind),                 intent(in)  :: zs
        real(rkind),                 intent(in)  :: hinterface
        real(rkind),                 intent(in)  :: Hb
        real(rkind), dimension(Nw+1),intent(in)  :: rhow, cw, alphaw
        real(rkind), dimension(Nb+1),intent(in)  :: rhob, cb, alphab
        real(rkind),                 intent(in)  :: ch, alphah
        real(rkind),                 intent(out) :: rhozs
        real(rkind),    allocatable, intent(out) :: r(:)
        complex(rkind), allocatable, intent(out) :: kw(:), kb(:)
        complex(rkind),              intent(out) :: kh
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
            z2 = (1.0_rkind - x2) * 0.5_rkind * (Hb - hinterface) + hinterface
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
        kh = w / ch * (1.0_rkind + ci * alphah / (40.0_rkind * pi * log10(exp(1.0_rkind))))
        
    end subroutine

    subroutine EigenValueVector(Nw,Nb,hinterface,Hb,kw,kb,kh,rhow,rhob,rhoh,&
        alphah,Lowerboundary,kr,eigvectorw,eigvectorb)
        implicit none
        integer,                     intent(in)     :: Nw, Nb
        character(len=1),            intent(in)     :: Lowerboundary
        real(rkind),                 intent(in)     :: hinterface, Hb
        real(rkind),                 intent(in)     :: rhow(Nw+1), rhob(Nb+1)
        real(rkind),                 intent(in)     :: rhoh
        real(rkind),                 intent(in)     :: alphah     
        complex(rkind),              intent(inout)  :: kw(Nw+1), kb(Nb+1)
        complex(rkind),              intent(in)     :: kh        
        complex(rkind),allocatable,  intent(out)    :: kr(:)
        complex(rkind),allocatable,  intent(out)    :: eigvectorw(:,:)
        complex(rkind),allocatable,  intent(out)    :: eigvectorb(:,:)
        real(rkind),    dimension(Nw+1)             :: Pu
        real(rkind),    dimension(Nb+1)             :: Pd        
        real(rkind),    dimension(Nw+1, Nw+1)       :: D1
        real(rkind),    dimension(Nb+1, Nb+1)       :: D2        
        real(rkind),    dimension(Nw+1)             :: x1
        real(rkind),    dimension(Nb+1)             :: x2    
        complex(rkind), dimension(Nw+1, Nw+1)       :: A
        complex(rkind), dimension(Nb+1, Nb+1)       :: B               
        complex(rkind), dimension(Nw+Nb+2, Nw+Nb+2) :: U
        complex(rkind), dimension(Nb+1)             :: Hb_boundary
                                                    
        complex(rkind), allocatable, dimension(:,:) :: V, W, E
        complex(rkind), allocatable, dimension(:,:) :: A_, B_        
        complex(rkind), allocatable, dimension(:,:) :: L11
        complex(rkind), allocatable, dimension(:,:) :: L12
        complex(rkind), allocatable, dimension(:,:) :: L21
        complex(rkind), allocatable, dimension(:,:) :: L22
        complex(rkind), allocatable, dimension(:,:) :: L
        complex(rkind), allocatable, dimension(:,:) :: v2
        complex(rkind), allocatable, dimension(:,:) :: GEVL, VR        
        complex(rkind), allocatable, dimension(:)   :: VL
        complex(rkind), allocatable, dimension(:)   :: WORK
        real(rkind)   , allocatable, dimension(:)   :: RWORK, IPIV                                                    
        complex(rkind), allocatable, dimension(:)   :: ALPHA
        complex(rkind), allocatable, dimension(:)   :: BETA                                                            
        integer                                     :: j(1), INFO 
        integer                                     :: i

        D1 = DerivationMatrix(Nw+1)
        D2 = DerivationMatrix(Nb+1)

        do i = 1, Nw + 1
            x1(i) = cos((i - 1) * pi / Nw)
            kw(i) = kw(i) ** 2            
        end do

        do i = 1, Nb + 1
            x2(i) = cos((i - 1) * pi / Nb)
            kb(i) = kb(i) ** 2        
        end do

        A = 4.0_rkind / hinterface ** 2 * matmul(Convolution(Cheb(rhow, x1)), D1)
        A = matmul(A, Convolution(Cheb(1.0_rkind / rhow, x1)))
        A = matmul(A, D1) + Convolution(Cheb(kw, x1))

        B = 4.0_rkind / (Hb - hinterface) ** 2 * matmul(Convolution(Cheb(rhob, x2)), D2)
        B = matmul(B, Convolution(Cheb(1.0_rkind / rhob, x2)))
        B = matmul(B, D2) + Convolution(Cheb(kb, x2))

        U = 0.0_rkind
        
        !for the second interface boundary
        do i = 1, Nw + 1
            Pu(i) = (-1.0_rkind) ** (i - 1)
        end do        
        Pu =  1.0_rkind / rhow(Nw+1) / hinterface * matmul(Pu, D1)
        Pd =  1.0_rkind
        Pd = -1.0_rkind / rhob(1) / (Hb - hinterface) * matmul(Pd, D2)
        !for the lower boundary
        do i = 1, Nb + 1
            Hb_boundary(i) = (-1.0_rkind) ** (i - 1)
        end do
       
        !Which type of lower boundary condition is used? 
        if (Lowerboundary == 'A') then
            U(1:Nw+1, 1:Nw+1)             = A
            U(Nw+2:Nw+Nb+2, Nw+2:Nw+Nb+2) = B
            allocate(V(Nw+Nb+2,Nw+Nb+2),W(Nw+Nb+2,Nw+Nb+2),E(Nw+Nb+2,Nw+Nb+2))
            V = 0.0_rkind
            W = 0.0_rkind
            E = 0.0_rkind
            do i = 1, Nw+Nb+2
                U(i, i) = U(i, i) - kh ** 2
                W(i, i) = 1.0_rkind
                E(i, i) = 1.0_rkind
            enddo
            !upper boundary
            U(Nw, 1:Nw+1) = 1.0_rkind
            W(Nw, Nw)     = 0.0_rkind            
            !first interface boundary
            do i = 1, Nw + 1
                U(Nw+1, i)   = (-1.0_rkind) ** (i - 1)
            end do
            U(Nw+1, Nw+2:Nw+Nb+2) = -1.0_rkind
            W(Nw+1, Nw+1)         =  0.0_rkind            
            !second interface boundary
            U(Nw+Nb+1, 1:Nw+1)       = Pu
            U(Nw+Nb+1, Nw+2:Nw+Nb+2) = Pd
            W(Nw+Nb+1, Nw+Nb+1)      = 0.0_rkind 
            !the most important lower boundary 
            U(Nw+Nb+2, Nw+2:Nw+Nb+2) = 2 * ci * rhoh / rhow(Nw+1) / (hinterface - Hb) * matmul(Hb_boundary, D2)            
            V(Nw+Nb+2, Nw+2:Nw+Nb+2) = Hb_boundary
            W(Nw+Nb+2, Nw+Nb+2)      = 0.0_rkind
 
            allocate(A_(2*(Nw+Nb+2),2*(Nw+Nb+2)), B_(2*(Nw+Nb+2),2*(Nw+Nb+2)))
            A_ = 0.0_rkind
            B_ = 0.0_rkind
            A_(1:Nw+Nb+2,       1:   Nw+Nb+2)  = -V
            A_(1:Nw+Nb+2, Nw+Nb+3:2*(Nw+Nb+2)) = -U
            A_(Nw+Nb+3:2*(Nw+Nb+2), 1:Nw+Nb+2) =  E
            B_(1:Nw+Nb+2,           1:Nw+Nb+2) =  W
            B_(Nw+Nb+3:2*(Nw+Nb+2), Nw+Nb+3:2*(Nw+Nb+2)) = E
            deallocate(V, W, E)
            allocate(ALPHA(2*(Nw+Nb+2)), BETA(2*(Nw+Nb+2)), GEVL(2*(Nw+Nb+2),2*(Nw+Nb+2)))
            allocate(VR(2*(Nw+Nb+2),2*(Nw+Nb+2)), WORK(4*(Nw+Nb+2)), RWORK(16*(Nw+Nb+2)))
           
            call zggev('N', 'V', 2*(Nw+Nb+2), A_, 2*(Nw+Nb+2), B_, 2*(Nw+Nb+2),&
                        ALPHA, BETA, GEVL, 2*(Nw+Nb+2), VR, 2*(Nw+Nb+2), WORK, &
                        4*(Nw+Nb+2), RWORK, INFO)                       
            ALPHA = ALPHA / BETA         
            
            do i = 1, 2*(Nw+Nb+2)
                if(real(ALPHA(i)) >= 0.0_rkind) then
                    ALPHA(  INFO+1) = sqrt(kh ** 2 - ALPHA(i) ** 2)
                    !GEVL temporarily stores the right eigenvectors 
                    GEVL(:, INFO+1) = VR(:, i)
                    INFO =  INFO+1
                endif
            end do
                        
            if(INFO > 1) then
                allocate(kr(INFO), eigvectorw(Nw+1,INFO), eigvectorb(Nb+1,INFO))
                ALPHA(INFO+1:2*(Nw+Nb+2)) = -9999.0_rkind
                !VR temporarily stores the sorted eigenvectors 
                do i = 1, INFO
                    j = maxloc(real(ALPHA))
                    kr(i)       = ALPHA(j(1))
                    VR(:, i)    = GEVL(:, j(1))
                    ALPHA(j(1)) = -9999.0_rkind
                end do                

                eigvectorw = VR(Nw+Nb+3:2*Nw+Nb+3, 1:INFO)
                eigvectorb = VR(2*Nw+Nb+4:2*(Nw+Nb+2), 1:INFO)
            else
                stop 'Error! No suitable modes!'
            endif           
            
            deallocate(ALPHA, BETA, GEVL, VR, WORK, RWORK)
            
        else
            U(1:Nw-1,       1:Nw-1)        = A(1:Nw-1, 1:Nw-1)
            U(1:Nw-1, Nw+Nb-1:Nw+Nb)       = A(1:Nw-1, Nw:Nw+1)
            U(Nw:Nw+Nb-2,  Nw:Nw+Nb-2)     = B(1:Nb-1, 1:Nb-1)
            U(Nw:Nw+Nb-2, Nw+Nb+1:Nw+Nb+2) = B(1:Nb-1, Nb:Nb+1)
            !upper boundary
            U(Nw+Nb-1, 1:Nw-1)        = 1.0_rkind
            U(Nw+Nb-1, Nw+Nb-1:Nw+Nb) = 1.0_rkind
            !lower boundary
            if (Lowerboundary == 'R')  Hb_boundary = matmul(Hb_boundary, D2)
            U(Nw+Nb+2, Nw:Nw+Nb-2)      = Hb_boundary(1 : Nb-1)
            U(Nw+Nb+2, Nw+Nb+1:Nw+Nb+2) = Hb_boundary(Nb: Nb+1)
            !first interface boundary
            do i = 1, Nw-1
                U(Nw+Nb, i) = (-1.0_rkind) ** (i - 1)
            end do
            U(Nw+Nb, Nw+Nb-1)         = (-1.0_rkind) **(Nw - 1)
            U(Nw+Nb, Nw+Nb)           = (-1.0_rkind) ** Nw
            U(Nw+Nb, Nw:Nw+Nb-2)      =  -1.0_rkind
            U(Nw+Nb, Nw+Nb+1:Nw+Nb+2) =  -1.0_rkind
            !second interface boundary
            U(Nw+Nb+1,  1:Nw-1)         = Pu(1 :Nw-1)
            U(Nw+Nb+1, Nw:Nw+Nb-2)      = Pd(1 :Nb-1)
            U(Nw+Nb+1, Nw+Nb-1:Nw+Nb)   = Pu(Nw:Nw+1)
            U(Nw+Nb+1, Nw+Nb+1:Nw+Nb+2) = Pd(Nb:Nb+1)

            allocate(L11(Nw+Nb-2,Nw+Nb-2), L12(Nw+Nb-2,4), L21(4,Nw+Nb-2),&
                     L22(4,4), L(Nw+Nb-2,Nw+Nb-2), v2(4, Nw+Nb-2), IPIV(4)) 
                                        
            L11 = U(1:Nw+Nb-2, 1:Nw+Nb-2)
            L12 = U(1:Nw+Nb-2, Nw+Nb-1:Nw+Nb+2)
            L21 = U(Nw+Nb-1:Nw+Nb+2, 1:Nw+Nb-2)
            L22 = U(Nw+Nb-1:Nw+Nb+2, Nw+Nb-1:Nw+Nb+2)

            allocate(kr(Nw+Nb-2), eigvectorw(Nw+1,Nw+Nb-2), eigvectorb(Nb+1,Nw+Nb-2))

            call zgesv(4, Nw+Nb-2, L22, 4, IPIV, L21, 4, INFO)
            
            L = L11 - matmul(L12, L21)
            
            allocate(VL(Nw+Nb-2), VR(Nw+Nb-2, Nw+Nb-2), WORK(2*(Nw+Nb-2)), RWORK(2*(Nw+Nb-2)))
            call zgeev('N', 'V', Nw+Nb-2, L, Nw+Nb-2, kr, VL, 1, VR, Nw+Nb-2, WORK, 2*(Nw+Nb-2), RWORK, INFO)

            v2 = -matmul(L21, VR)
            kr = sqrt(kr)  
            eigvectorw = 0.0_rkind
            eigvectorb = 0.0_rkind
            eigvectorw(1:Nw-1,  :) = VR(1:Nw-1,     :)
            eigvectorw(Nw:Nw+1, :) = v2(1:2,        :)
            eigvectorb(1:Nb-1,  :) = VR(Nw:Nw+Nb-2, :)
            eigvectorb(Nb:Nb+1, :) = v2(3:4,        :) 
            
            !L and L11 temporarily store the sorted eigenvectors respectively 
            !VL temporarily stores the sorted eigenvalues
            do i = 1, Nw+Nb-2
                j = maxloc(real(kr))
                VL(i) = kr(j(1))
                L  (1:Nw+1, i) = eigvectorw(:, j(1))
                L11(1:Nb+1, i) = eigvectorb(:, j(1))
                kr(j(1)) = -9999.0_rkind
            end do

            kr = VL
            eigvectorw = L  (1:Nw+1, :)
            eigvectorb = L11(1:Nb+1, :)            

            deallocate(L11, L12, L21, L22, L, v2, IPIV)           
        endif

    end subroutine

    subroutine NumofModes(freq,kr,nmodes,cpmax)
       
        implicit none
        real(rkind),    intent(in) :: freq
        real(rkind),    intent(in) :: cpmax
        complex(rkind), intent(in) :: kr(:)
        integer,        intent(out):: nmodes        
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
        
    end subroutine

    subroutine GenerateModes(nmodes,dz,hinterface,Hb,eigvectorw,eigvectorb,psi,z)

        implicit none
        integer,                         intent(in)  :: nmodes
        real(rkind),                     intent(in)  :: dz
        real(rkind),    allocatable,     intent(out) :: z(:)    
        real(rkind),                     intent(in)  :: hinterface, Hb        
        complex(rkind), dimension(:, :), intent(in)  :: eigvectorw, eigvectorb        
        complex(rkind), allocatable,     intent(out) :: psi(:, :)
        real(rkind),    allocatable, dimension(:)    :: xt1, xt2, zt1, zt2, zt3
        complex(rkind), allocatable, dimension(:, :) :: psi1, psi2, psit
        integer                                      :: i

        allocate (zt1(ceiling(hinterface / dz)))
        allocate (zt2(ceiling((Hb - hinterface) / dz)))
        allocate (z  (ceiling(Hb / dz) + 1) )

        do i = 1, size(zt1)
            zt1(i) = (i - 1) * dz
        end do
        do i = 1, size(zt2)
            zt2(i) = hinterface + (i - 1) * dz
        end do
        do i = 1, size(z) - 1
            z(i) = (i - 1) * dz
        end do
        z(size(z)) = Hb

        xt1 = -2.0_rkind / hinterface * zt1 + 1.0_rkind
        xt2 = -2.0_rkind / (Hb - hinterface) * zt2 + (Hb + hinterface) / (Hb - hinterface)

        allocate(psi1(size(xt1), nmodes), psi2(size(xt2), nmodes), psi(size(z), nmodes))
        psi1 = InvChebMatrix(eigvectorw(:, 1:nmodes), xt1)
        psi2 = InvChebMatrix(eigvectorb(:, 1:nmodes), xt2)    
            
        allocate(zt3(size(zt1)+size(zt2)),psit(size(zt1)+size(zt2),nmodes))
        zt3(1:size(zt1)) = zt1
        zt3(size(zt1)+1:size(zt1)+size(zt2)) = zt2
        psit(1:size(zt1),:) = psi1
        psit(size(zt1)+1:size(zt1)+size(zt2),:) = psi2    
        
        do i =1, nmodes
            call Interpolation_SSP(zt3,psit(:,i),z,psi(:,i))
        end do

        deallocate(psi1, psi2, psit, xt1, xt2, zt1, zt2, zt3)

    end subroutine

    subroutine Normalization(eigvectorw,eigvectorb,rhow,rhob,rhoh,kr,kh,Lowerboundary,&
                             hinterface,Hb,Nw,Nb,nmodes,psi)

        implicit none
        character(len=1),                intent(in)  :: Lowerboundary
        complex(rkind), dimension(:, :), intent(in)  :: eigvectorw, eigvectorb
        integer,                         intent(in)  :: Nw, Nb
        integer,                         intent(in)  :: nmodes
        real(rkind),                     intent(in)  :: rhoh
        real(rkind),                     intent(in)  :: hinterface, Hb
        complex(rkind), dimension(:, :), intent(out) :: psi        
        real(rkind),                     intent(out) :: rhow(Nw+1)        
        real(rkind),                     intent(out) :: rhob(Nb+1)
        complex(rkind),                  intent(in)  :: kh
        complex(rkind),                  intent(in)  :: kr(nmodes)      
        real(rkind)                                  :: x1(Nw+1), P(Nw+1)
        real(rkind)                                  :: x2(Nb+1), Q(Nb+1)        
        real(rkind)                                  :: Co11(Nw+1, Nw+1), Co22(Nb+1, Nb+1)    
        complex(rkind)                               :: Co1 (Nw+1, Nw+1), Co2 (Nb+1, Nb+1)
        complex(rkind)                               :: f1(Nw+1), f2(Nb+1)
        complex(rkind)                               :: norm
        integer                                      :: i
        
        do i = 1, Nw + 1
            x1(i) = cos( (i-1) * pi / Nw)        
        end do
        do i = 1, Nb + 1
            x2(i) = cos( (i-1) * pi / Nb)
        end do

        Co11 = Convolution(Cheb(1.0_rkind / rhow, x1))
        Co22 = Convolution(Cheb(1.0_rkind / rhob, x2))

        P = 0.0_rkind
        Q = 0.0_rkind
        do i = 0, Nw-1, 2
            P(i+1) = - 2.0_rkind / (i * i - 1.0_rkind)
        end do
        do i = 0, Nb-1, 2
            Q(i+1) = - 2.0_rkind / (i * i - 1.0_rkind)
        end do

        do i = 1, nmodes
            Co1 = Convolution(eigvectorw(:, i))

            f1  = matmul(Co11, matmul(Co1, eigvectorw(:, i)))

            Co2 = Convolution(eigvectorb(:, i))
            f2  = matmul(Co22, matmul(Co2, eigvectorb(:, i)))

            norm = dot_product(P, f1) * hinterface * 0.5_rkind + &
                   dot_product(Q, f2) * (Hb - hinterface) * 0.5_rkind

            if(Lowerboundary == 'A') then
                norm = norm + 0.5_rkind / rhoh * psi(size(psi, 1), i) ** 2  / sqrt(kh ** 2 - kr(i) ** 2)
            end if

            psi(:, i) = psi(:, i) / sqrt(norm)
        end do

    end subroutine

    subroutine SynthesizeSoundField(nmodes,nr,r,kr,rhozs,zs,dz,psi,tl)

        implicit none
        integer,                  intent(in)    :: nmodes, nr
        real(rkind),              intent(in)    :: r(nr), rhozs, zs, dz
        real(rkind), allocatable, intent(inout) :: tl(:, :)
        complex(rkind),           intent(in)    :: kr(:)
        complex(rkind),           intent(inout) :: psi(:, :)
        complex(rkind), allocatable             :: p(:, :)
        complex(rkind)                          :: bessel(nmodes, nr)
        complex(rkind)                          :: psizs(nmodes, nmodes)
        real(rkind)                             :: CYR, CYI
        integer                                 :: IERR1, IERR2        
        integer                                 :: i, k, s

        do k = 1, nmodes
        do i = 1, nr
            bessel(k, i) = r(i) * kr(k)
            call ZBESH(real(bessel(k, i)), aimag(bessel(k, i)), 0.0_rkind, 1, 1, 1, CYR, CYI, IERR1, IERR2)
            bessel(k, i) = cmplx(CYR, CYI)
        end do
        end do
            
        allocate(p(size(psi, 1), nr), tl(size(psi, 1), nr))
        s = ceiling(zs / dz)
        
        psizs(:, :) = 0.0_rkind
        do k = 1, nmodes
            psizs(k, k) = (zs / dz - s) * psi(s + 1, k) + ( s + 1 - zs / dz) * psi(s, k)
        end do
        
        psi = matmul(psi, psizs)
        p   = matmul(psi, bessel)
        p   = p * ci * pi / rhozs
        
        tl  = - 20.0_rkind * log10(abs(p))
        
    end subroutine

    subroutine SaveSoundField(filename,casename,tlmin,tlmax,r,z,tl)

        implicit none
        character(len=MAX_FILENAME_LEN), intent(in) :: filename
        character(len=MAX_FILENAME_LEN), intent(in) :: casename
        real(rkind),                     intent(in) :: tlmin, tlmax
        real(rkind), dimension(:),       intent(in) :: r, z
        real(rkind), dimension(:, :),    intent(in) :: tl

        open(unit=20, status='unknown', file=filename, access='stream', form='unformatted')
        write(20)  casename, size(z), size(r), tlmin, tlmax, z, r, tl
        close(20)
        
    end subroutine

end module
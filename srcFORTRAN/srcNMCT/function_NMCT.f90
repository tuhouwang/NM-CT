module param_mod

    implicit none
    integer,        parameter :: rkind            = 8
    integer,        parameter :: MAX_FILENAME_LEN = 200
    real(rkind),    parameter :: pi = 4.0_rkind * atan(1.0_rkind)
    complex(rkind), parameter :: ci = cmplx(0.0, 1.0, kind=rkind)
    integer, 		parameter :: LL	= 5
	integer, 		parameter :: CC	= 2500

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
        integer,        intent(in) :: m
        real(rkind)                :: DerivationMatrix(m, m)
        integer                    :: i, j
    
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
    
    interface Interpolation_zs
        module procedure Interpolation_zs_R
        module procedure Interpolation_zs_C
    end interface

contains

    function Interpolation_zs_R(z,zs,v)
        implicit none
        real   (rkind),intent(in) :: z(:), zs
        real   (rkind),intent(in) :: v(:)
        real   (rkind)            :: Interpolation_zs_R
        integer                   :: j
        Interpolation_zs_R = 0.0_rkind
        
        call assert(zs <= z(size(z)) .and. zs >= z(1), 'Error! zs is not in the range of z, interpolation failed!')
        
        call assert(size(z)==size(v),'Error! the dimensions of the interpolation array do not match!')
            
        do j=1, size(z)-1
            if (zs >= z(j) .and. zs <= z(j+1)) then
                Interpolation_zs_R=(zs-z(j))/(z(j+1)-z(j))*v(j+1)+(z(j+1)-zs)/(z(j+1)-z(j))*v(j)
            endif				
        end do 

    end function

    function Interpolation_zs_C(z,zs,v)
        implicit none
        real   (rkind),intent(in) :: z(:), zs
        complex(rkind),intent(in) :: v(:)
        complex(rkind)            :: Interpolation_zs_C
        integer                   :: j
        Interpolation_zs_C = 0.0_rkind
        
        call assert(zs <= z(size(z)) .and. zs >= z(1), 'Error! zs is not in the range of z, interpolation failed!')
        
        call assert(size(z)==size(v),'Error! the dimensions of the interpolation array do not match!')
        
        do j=1, size(z)-1
            if (zs >= z(j) .and. zs <= z(j+1)) then
                Interpolation_zs_C=(zs-z(j))/(z(j+1)-z(j))*v(j+1)+(z(j+1)-zs)/(z(j+1)-z(j))*v(j)
            endif				
        end do 

    end function

    subroutine ReadEnvParameter(casename,Layers,Ns,cpmax,freq,zs,zr,rmax,dr,hinterface,dz,&
        tlmin,tlmax,dep,c,rho,alpha,ch,rhoh,alphah,Lowerboundary,data_file)
        
        implicit none
        character(len=MAX_FILENAME_LEN), intent(out) :: casename
        character(len=MAX_FILENAME_LEN), intent(in)  :: data_file
        character(len=1),                intent(out) :: Lowerboundary
        integer ,                        intent(out) :: Layers
        integer ,      allocatable,      intent(out) :: Ns(:)
        real(rkind),                     intent(out) :: cpmax
        real(rkind),                     intent(out) :: freq        
        real(rkind),                     intent(out) :: zs
        real(rkind),                     intent(out) :: zr
        real(rkind),                     intent(out) :: rmax
        real(rkind),                     intent(out) :: dr
        real(rkind),   allocatable,      intent(out) :: hinterface(:)
        real(rkind),                     intent(out) :: dz
        real(rkind),                     intent(out) :: tlmin
        real(rkind),                     intent(out) :: tlmax
        real(rkind),dimension(LL,CC),    intent(out) :: dep
        real(rkind),dimension(LL,CC),    intent(out) :: c        
        real(rkind),dimension(LL,CC),    intent(out) :: rho    
        real(rkind),dimension(LL,CC),    intent(out) :: alpha        
        real(rkind),                     intent(out) :: ch        
        real(rkind),                     intent(out) :: rhoh 
        real(rkind),                     intent(out) :: alphah  
        integer,      allocatable			         :: nprofile(:)
        integer                                      :: i, j
        dep   = 0.0_rkind
        c     = 9999.0_rkind
        rho   = 0.0_rkind
        alpha = 0.0_rkind
                                        
        open(unit=1, status='unknown', file=data_file) 
    
        read (1, '(A200)') casename
        read (1, *) Layers
        
        call assert(Layers <= LL, 'Error! Layers must less than or equal to 5!')
      
        allocate(Ns(Layers), hinterface(Layers), nprofile(Layers)) 
        
        do i=1, Layers
            read(1, *) Ns(i)
        enddo       
        read (1, *) cpmax
        read (1, *) freq
        read (1, *) zs
        read (1, *) zr
        read (1, *) rmax
        read (1, *) dr
        
        call assert(rmax/dr - int(rmax/dr) == 0, 'Please reinput the dr and rmax!')

        do i=1, Layers
            read(1, *) hinterface(i)
        enddo
        
        do i=2, Layers
            call assert(hinterface(i) >= hinterface(i-1), 'Error! h must greater than 0 and less than H!')
        enddo
        
        call assert(zs > 0 .and. zs < hinterface(Layers) .and. zr > 0 .and. zr < hinterface(Layers),&
                    'zs and zr must be greater than 0 and less than H!')
                    
        read (1, *) dz        
        read (1, *) tlmin
        read (1, *) tlmax
        call assert(tlmin < tlmax, 'tlmin must less than tlmax!')
        
        call assert(hinterface(Layers)/dz-floor(hinterface(Layers)/dz) == 0, 'Error! The input dz unsuitable!')      
                
        do i=1, Layers
            read(1, *) nprofile(i)
            if (nprofile(i) > CC) then
               stop 'Error! nprofiles must less than or equal to 2500!'	
            end if		   
        enddo
        
        do i=1, Layers
            do j=1, nprofile(i)
                read(1, *) dep(i,j), c(i,j), rho(i,j), alpha(i,j)
            enddo	
        enddo
        
        do i=2, Layers
            call assert(hinterface(i) == dep(i,nprofile(i)) .and. hinterface(i-1) == dep(i,1), &
                        'Error! input sound profile is unsuitable !')			
        enddo
        
        do i=1, Layers
            call ChebInterpolation(nprofile(i),dep(i,:),c(i,:),rho(i,:),alpha(i,:),Ns(i))	
        enddo	
        
        read (1, *) Lowerboundary
        
        call assert(Lowerboundary == 'V' .or. Lowerboundary == 'R' .or. Lowerboundary == 'A',&
                    'Error! The lower boundary must be vaccum, rigid or half-space!') 
             
        if (Lowerboundary == 'A') then
            read(1, *) ch, rhoh, alphah
        else
            ch     = 0.0_rkind
            rhoh   = 0.0_rkind
            alphah = 0.0_rkind
        endif
        
        close(1)

    end subroutine

    subroutine ChebInterpolation(m,dep,c,rho,alpha,N)
        implicit none
        integer,       intent(in)    :: N
        integer,       intent(in)    :: m
        real(rkind),   intent(inout) :: dep(N+1)
        real(rkind),   intent(inout) :: c(N+1)
        real(rkind),   intent(inout) :: rho(N+1)
        real(rkind),   intent(inout) :: alpha(N+1)        
        real(rkind)                  :: x(N+1)
        real(rkind)                  :: z(N+1)
        real(rkind)                  :: b2(N+1)
        real(rkind)                  :: c2(N+1)
        real(rkind)                  :: d2(N+1)
        integer                      :: i,j
        
        do i = 1, N+1
            x(i) = cos((i - 1) * pi / N)
            z(i) = ((dep(m) + dep(1)) / (dep(m) - dep(1)) - x(i)) * (dep(m) - dep(1)) * 0.5_rkind
        end do

        do i=1, N+1	
            do j=1, m-1
                if((z(i) >= dep(j)) .and. (z(i) <= dep(j+1))) then
                    b2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * c(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * c(j)
                    
                    c2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * rho(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * rho(j)
                    
                    d2(i) = (z(i) - dep(j)) / (dep(j+1) - dep(j)) * alpha(j+1)&
                            +(dep(j+1)-z(i))/ (dep(j+1) - dep(j)) * alpha(j)				              
                end if				
            end do 
            
            if(z(i) > dep(m)) then
                b2(i) = c(m)
                c2(i) = rho(m)
                d2(i) = alpha(m)							
            end if
                
            if(z(i) < dep(1)) then
                b2(i) = c(1)
                c2(i) = rho(1)	
                d2(i) = alpha(1)					
            endif
        end do
        
        dep   = z
        c     = b2
        rho   = c2
        alpha = d2
        
    end subroutine

    subroutine Interpolation_SSP(dep,c,dep2,c2)
        implicit none
        real(rkind),   intent(in) :: dep(:)
        real(rkind),   intent(in) :: dep2(:)
        complex(rkind),intent(in) :: c(:)
        complex(rkind),intent(out):: c2(:)
        integer                   :: n_b
        integer                   :: n_b2        
        integer                   :: i, j
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

    subroutine Initialization(Layers,Ns,freq,rmax,dr,zs,nr,r,rhozs,hinterface,ki,&
                                                    dep,c,rho,alpha,ch,alphah,kh)
        implicit none
        integer,                        intent(in)  :: Layers
        integer,                        intent(in)  :: Ns(Layers)
        real(rkind),                    intent(in)  :: freq
        real(rkind),                    intent(in)  :: rmax
        real(rkind),                    intent(in)  :: dr
        real(rkind),                    intent(in)  :: zs        
        integer,                        intent(out) :: nr
        real(rkind),    allocatable,    intent(out) :: r(:)
        real(rkind),                    intent(out) :: rhozs
        real(rkind),                    intent(in)  :: hinterface(Layers)
        complex(rkind),dimension(LL,CC),intent(out) :: ki
        real(rkind),   dimension(LL,CC),intent(in)  :: dep
        real(rkind),   dimension(LL,CC),intent(in)  :: c        
        real(rkind),   dimension(LL,CC),intent(in)  :: rho    
        real(rkind),   dimension(LL,CC),intent(in)  :: alpha          
        real(rkind),                    intent(in)  :: ch    
        real(rkind),                    intent(in)  :: alphah   
        complex(rkind),                 intent(out) :: kh
        integer                                     :: i, j

        nr = int(rmax / dr)

        allocate (r(nr))
        do i = 1, nr
            r(i) = i * dr
        end do
        
        ki = 2.0_rkind*pi*freq/c *(1.0+ci*alpha /(40.0*pi*log10(exp(1.0))))
        kh = 2.0_rkind*pi*freq/ch*(1.0+ci*alphah/(40.0*pi*log10(exp(1.0))))

        if( zs <= hinterface(1) ) then       
            rhozs = Interpolation_zs(dep(1,1:Ns(1)+1),zs,rho(1,1:Ns(1)+1))
        end if	
        do i = 2, Layers
            if( zs <= hinterface(i) .and. zs >= hinterface(i-1)) then     
                rhozs = Interpolation_zs(dep(i,1:Ns(i)+1),zs,rho(i,1:Ns(i)+1))
            end if
        end do
        
    end subroutine

    subroutine EigenValueVector(Layers,Ns,dep,ki,rho,kh,rhoh,Lowerboundary,kr,eigvector)
        implicit none
        integer,                        intent(in)   :: Layers
        integer,                        intent(in)   :: Ns(Layers)
        real(rkind),   dimension(LL,CC),intent(in)   :: dep
        real(rkind),   dimension(LL,CC),intent(in)   :: rho        
        complex(rkind),dimension(LL,CC),intent(in)   :: ki
        complex(rkind),                 intent(in)   :: kh
        real(rkind),                    intent(in)   :: rhoh        
        character(len=1),               intent(in)   :: Lowerboundary
        complex(rkind), allocatable,    intent(out)  :: kr(:)
        complex(rkind), allocatable,    intent(out)  :: eigvector(:, :)
        complex(rkind), allocatable, dimension(:)    :: k
        real(rkind),    allocatable, dimension(:, :) :: D
        complex(rkind), allocatable, dimension(:, :) :: A
        complex(rkind), allocatable, dimension(:, :) :: B
        complex(rkind), allocatable, dimension(:, :) :: U
        complex(rkind), allocatable, dimension(:, :) :: V
        complex(rkind), allocatable, dimension(:, :) :: W
        complex(rkind), allocatable, dimension(:, :) :: E                
        real(rkind)   , allocatable, dimension(:)    :: x    
        real(rkind)   , allocatable, dimension(:)    :: Pu  
        real(rkind)   , allocatable, dimension(:)    :: Pd                                                       
        complex(rkind), allocatable, dimension(:,:)  :: L11
        complex(rkind), allocatable, dimension(:,:)  :: L12
        complex(rkind), allocatable, dimension(:,:)  :: L21
        complex(rkind), allocatable, dimension(:,:)  :: L22
        complex(rkind), allocatable, dimension(:,:)  :: L
        complex(rkind), allocatable, dimension(:,:)  :: v2                                                     
        complex(rkind), allocatable, dimension(:,:)  :: GEVL       
        complex(rkind), allocatable, dimension(:)    :: VL
        complex(rkind), allocatable, dimension(:,:)  :: VR 
        complex(rkind), allocatable, dimension(:)    :: WORK
        real(rkind)   , allocatable, dimension(:)    :: RWORK
        real(rkind)   , allocatable, dimension(:)    :: IPIV        
        complex(rkind), allocatable, dimension(:)    :: ALPHA
        complex(rkind), allocatable, dimension(:)    :: BETA                                                            
        integer                                      :: j(1) 
        integer                                      :: INFO
        integer                                      :: n      
        integer                                      :: i
        
        allocate(U(sum(Ns+1), sum(Ns+1)))
        U = 0.0_rkind
        if(Lowerboundary == 'A') then
            INFO = 1
            do i = 1, Layers
                allocate(D(Ns(i)+1,Ns(i)+1), A(Ns(i)+1,Ns(i)+1), x(Ns(i)+1), k(Ns(i)+1))
                D = DerivationMatrix(Ns(i)+1)
                do n = 1, Ns(i) + 1
                    x(n) = cos((n - 1) * pi / Ns(i))
                    k(n) = ki(i, n) ** 2            
                end do
                A = 4.0_rkind/(dep(i,Ns(i)+1)-dep(i,1))**2*matmul(Convolution(Cheb(rho(i,1:Ns(i)+1), x)), D)    
                A = matmul(A, Convolution(Cheb(1.0_rkind/rho(i,1:Ns(i)+1), x)))
                A = matmul(A, D) + Convolution(Cheb(k(1:Ns(i)+1), x))
                do n = 1, Ns(i) + 1
                    A(n, n) = A(n, n) - kh ** 2       
                end do            
                U(INFO:INFO+Ns(i), INFO:INFO+Ns(i)) = A
                INFO = INFO + Ns(i) + 1
                deallocate(D,A,x,k)
            end do
            
            allocate(V(sum(Ns+1),sum(Ns+1)), W(sum(Ns+1),sum(Ns+1)), E(sum(Ns+1),sum(Ns+1)))
            V = 0.0_rkind
            W = 0.0_rkind
            E = 0.0_rkind
            do i = 1, sum(Ns+1)
                W(i, i) = 1.0_rkind
                E(i, i) = 1.0_rkind
            end do
            
            ! boundary condition
            INFO = 1
            do i = 1, Layers - 1
                allocate(Pu(Ns(i)+1), Pd(Ns(i+1)+1))
                do n = 0, Ns(i)
                    Pu(n+1) = (-1.0_rkind) ** n
                end do
                Pd = 1.0_rkind

                ! sound pressure is continuous
                U(INFO+Ns(i)-1, INFO:INFO+Ns(i)                  ) = Pu
                U(INFO+Ns(i)-1, INFO+Ns(i)+1:INFO+Ns(i)+Ns(i+1)+1) = -1.0_rkind
                W(INFO+Ns(i)-1, INFO+Ns(i)-1                     ) = 0.0_rkind
                
                Pu = -1.0_rkind / rho(i,Ns(i)+1) / (dep(i,Ns(i)+1) - dep(i,1)) * matmul(Pu, DerivationMatrix(Ns(i)+1))
                Pd =  1.0_rkind / rho(i+1,1) / (dep(i+1,Ns(i)+1) - dep(i+1,1)) * matmul(Pd, DerivationMatrix(Ns(i+1)+1))
                
                ! normal velocity is continuous
                U(INFO+Ns(i),         INFO:INFO+Ns(i))           = Pu
                U(INFO+Ns(i), INFO+Ns(i)+1:INFO+Ns(i)+Ns(i+1)+1) = Pd
                W(INFO+Ns(i), INFO+Ns(i)) = 0.0_rkind
                
                deallocate(Pu, Pd)                
                INFO = INFO + Ns(i) + 1              
            end do            
            
            ! upper boundary
            U(size(U,1)-1,   1:Ns(1)+1) = 1.0_rkind
            W(size(W,1)-1, size(W,1)-1) = 0.0_rkind
            ! the most important lower boundary
            allocate(Pu(Ns(Layers)+1))
            do n = 0, Ns(Layers)
                Pu(n+1) = (-1.0_rkind) ** n
            end do
            
            U(size(U,1), size(U,1)-Ns(Layers):size(U,1)) = -2.0_rkind * ci * rhoh / rho(Layers, Ns(Layers)+1) / &
                          (dep(Layers,Ns(Layers)+1) - dep(Layers, 1)) * matmul(Pu, DerivationMatrix(Ns(Layers)+1))
            V(size(V,1), size(V,1)-Ns(Layers):size(V,1)) = Pu
            W(size(W,1), size(W,1)) = 0.0_rkind
            deallocate(Pu)
            
            allocate(A(2*sum(Ns+1),2*sum(Ns+1)), B(2*sum(Ns+1),2*sum(Ns+1)))
            A = 0.0_rkind
            B = 0.0_rkind
            A(1:sum(Ns+1),       1    :  sum(Ns+1)) = -V
            A(1:sum(Ns+1), sum(Ns+1)+1:2*sum(Ns+1)) = -U
            A(sum(Ns+1)+1:2*sum(Ns+1), 1:sum(Ns+1)) =  E
            B(1:sum(Ns+1),       1    :  sum(Ns+1)) =  W
            B(sum(Ns+1)+1:2*sum(Ns+1),   sum(Ns+1)+1:2*sum(Ns+1)) = E
            deallocate(U, V, W, E)
            
            allocate(ALPHA(2*sum(Ns+1)), BETA(2*sum(Ns+1)), GEVL(2*sum(Ns+1),2*sum(Ns+1)))
            allocate(VR(2*sum(Ns+1),2*sum(Ns+1)), WORK(4*sum(Ns+1)), RWORK(16*sum(Ns+1)))
           
            call zggev('N', 'V', 2*sum(Ns+1), A, 2*sum(Ns+1), B, 2*sum(Ns+1), &
                        ALPHA, BETA, GEVL, 2*sum(Ns+1), VR, 2*sum(Ns+1), WORK,&
                        4*sum(Ns+1), RWORK, INFO)                       
            ALPHA = ALPHA / BETA         
            
            do i = 1, 2*sum(Ns+1)
                if(real(ALPHA(i)) >= 0.0_rkind .and. real(sqrt(kh**2-ALPHA(i)**2)) <= maxval(real(ki))) then
                    ALPHA(  INFO+1) = sqrt(kh ** 2 - ALPHA(i) ** 2)
                    !GEVL temporarily stores the right eigenvectors 
                    GEVL(:, INFO+1) = VR(:, i)
                    INFO =  INFO+1
                endif
            end do
                        
            if(INFO > 1) then
                allocate(kr(INFO), eigvector(sum(Ns+1),INFO))
                eigvector = 0.0_rkind
                ALPHA(INFO+1:2*sum(Ns+1)) = -9999.0_rkind
                !VR temporarily stores the sorted eigenvectors 
                do i = 1, INFO
                    j = maxloc(real(ALPHA))
                    kr(i)       = ALPHA(j(1))
                    VR(:, i)    = GEVL(:, j(1))
                    ALPHA(j(1)) = -9999.0_rkind
                end do

                eigvector = VR(sum(Ns+1)+1:2*sum(Ns+1), 1:INFO) 
            else
                stop 'Error! No suitable modes!'
            endif           
            
            deallocate(A, B, ALPHA, BETA, GEVL, VR, WORK, RWORK)
            
        else
            INFO = 1
            do i = 1, Layers
                allocate(D(Ns(i)+1,Ns(i)+1), A(Ns(i)+1,Ns(i)+1), x(Ns(i)+1), k(Ns(i)+1))
                D = DerivationMatrix(Ns(i)+1)
                do n = 1, Ns(i) + 1
                    x(n) = cos((n - 1) * pi / Ns(i))
                    k(n) = ki(i, n) ** 2            
                end do

                A = 4.0_rkind/(dep(i,Ns(i)+1)-dep(i,1))**2*matmul(Convolution(Cheb(rho(i,1:Ns(i)+1), x)), D)    
                A = matmul(A, Convolution(Cheb(1.0_rkind/rho(i,1:Ns(i)+1), x)))
                A = matmul(A, D) + Convolution(Cheb(k(1:Ns(i)+1), x))
         
                U(INFO:INFO+Ns(i)-2,            INFO:INFO+Ns(i)-2 ) = A(1:Ns(i)-1,    1:Ns(i)-1)
                U(INFO:INFO+Ns(i)-2, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = A(1:Ns(i)-1,Ns(i):Ns(i)+1)
                INFO = INFO + Ns(i) - 1
                deallocate(D, A, x, k)
            end do 
            ! boundary condition
            INFO = 1
            do i = 1, Layers - 1
                allocate(Pu(Ns(i)+1), Pd(Ns(i+1)+1))
                do n = 0, Ns(i)
                    Pu(n+1) = (-1.0_rkind) ** n
                end do
                Pd = 1.0_rkind

                ! sound pressure is continuous
                U(sum(Ns-1)+2*i-1,            INFO:INFO+Ns(i)-2 ) = Pu(1    :Ns(i)-1)
                U(sum(Ns-1)+2*i-1, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = Pu(Ns(i):Ns(i)+1)
                
                Pu = -1.0_rkind / rho(i,Ns(i)+1) / (dep(i,Ns(i)+1) - dep(i,1)) * matmul(Pu, DerivationMatrix(Ns(i)+1))
                Pd =  1.0_rkind / rho(i+1,1) / (dep(i+1,Ns(i)+1) - dep(i+1,1)) * matmul(Pd, DerivationMatrix(Ns(i+1)+1))
                
                ! normal velocity is continuous
                U(sum(Ns-1)+2*i,            INFO:INFO+Ns(i)-2 ) = Pu(1    :Ns(i)-1)
                U(sum(Ns-1)+2*i, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = Pu(Ns(i):Ns(i)+1)
                
                INFO = INFO + Ns(i) - 1
                ! sound pressure is continuous
                U(sum(Ns-1)+2*i-1,            INFO:INFO+Ns(i+1)-2 ) = -1.0_rkind
                U(sum(Ns-1)+2*i-1, sum(Ns-1)+2*i+1:sum(Ns-1)+2*i+2) = -1.0_rkind
                ! normal velocity is continuous
                U(sum(Ns-1)+2*i,            INFO:INFO+Ns(i+1)-2 ) = Pd(      1:Ns(i+1)-1)
                U(sum(Ns-1)+2*i, sum(Ns-1)+2*i+1:sum(Ns-1)+2*i+2) = Pd(Ns(i+1):Ns(i+1)+1)
                          
                deallocate(Pu, Pd)                           
            end do                    
            ! upper boundary, pressure-free boundary
            U(sum(Ns+1)-1, 1          :    Ns(1)-1) = 1.0_rkind
            U(sum(Ns+1)-1, sum(Ns-1)+1:sum(Ns-1)+2) = 1.0_rkind
            ! lower boundary, perfectly free / rigid
            allocate(Pd(Ns(Layers)+1))
            do n = 0, Ns(Layers)
                Pd(n+1) = (-1.0_rkind) ** n
            end do

            if(Lowerboundary == 'R') then
                Pd = matmul(Pd, DerivationMatrix(Ns(Layers)+1))
            end if

            U(sum(Ns+1), sum(Ns-1)-Ns(Layers)+2:sum(Ns-1)) = Pd(1         :Ns(Layers)-1)
            U(sum(Ns+1), sum(Ns+1)-1           :sum(Ns+1)) = Pd(Ns(Layers):Ns(Layers)+1)

            allocate(L11(sum(Ns-1),sum(Ns-1)), L12(sum(Ns-1),2*Layers), L21(2*Layers,sum(Ns-1)),&
                     L22(2*Layers,2*Layers), L(sum(Ns-1),sum(Ns-1)), VL(sum(Ns-1)), kr(sum(Ns-1)),&
                     eigvector(sum(Ns+1),sum(Ns-1)), v2(2*Layers,sum(Ns-1)), IPIV(2*Layers))
            !blocking
            L11 = U(1          :sum(Ns-1), 1          :sum(Ns-1))   
            L12 = U(1          :sum(Ns-1), sum(Ns-1)+1:sum(Ns+1))
            L21 = U(sum(Ns-1)+1:sum(Ns+1), 1          :sum(Ns-1))
            L22 = U(sum(Ns-1)+1:sum(Ns+1), sum(Ns-1)+1:sum(Ns+1))
            
            call zgesv(2*Layers,sum(Ns-1),L22,2*Layers,IPIV,L21,2*Layers,INFO)	
            L = L11 - matmul(L12,L21)
            allocate(VR(sum(Ns-1),sum(Ns-1)),WORK(2*sum(Ns-1)),RWORK(2*sum(Ns-1)))
            call zgeev('N','V',sum(Ns-1),L,sum(Ns-1),kr,VL,1,VR,sum(Ns-1),WORK,2*sum(Ns-1),RWORK,INFO)  	
            v2 = - matmul(L21,VR)
            deallocate(U, L11, L12, L21, L22, L, IPIV, WORK, RWORK)            
                        
            kr = sqrt(kr)

            INFO = 0
            do i = 1, Layers
                eigvector(INFO+i:        INFO+Ns(i)+i-2,:) = VR(INFO-i+2:INFO+Ns(i)-i, :)
                eigvector(INFO+Ns(i)+i-1:INFO+Ns(i)+i  ,:) = v2(   2*i-1:2*i, :)
                INFO = INFO + Ns(i)
            end do
            
            !L and VL store the sorted eigenvectors and eigenvalus, respectively.
            allocate(L(sum(Ns+1),sum(Ns-1)))
            do i = 1, sum(Ns-1)
                j       = maxloc(real(kr))
                VL(i)   = kr(j(1))
                L(:,i)  = eigvector(:,j(1))
                kr(j(1))= -9999.0_rkind
            end do
            
            kr        = VL
            eigvector = L		
          
            deallocate(VR,v2,L,VL)                             
        endif            

    end subroutine

    subroutine NumOfModes(Layers,freq,cpmax,kr,eigvector,nmodes)
        implicit none
        integer,       intent(in)   :: Layers
        integer,       intent(out)  :: nmodes	
        real(rkind),   intent(in)   :: freq
        real(rkind),   intent(in)   :: cpmax
        complex(rkind),intent(inout):: kr(:)
        complex(rkind),intent(inout):: eigvector(:, :)
        real(rkind),   allocatable  :: cp(:)
        integer                     :: i
        
        allocate(cp(size(kr)))
        cp = 2 * pi * freq / real(kr)
        
        where(cp > cpmax) kr = -9999.0_rkind 
        
        nmodes = 0
        do i=1, size(kr)
            if (kr(i) /= -9999.0_rkind) then
               nmodes               = nmodes + 1
               kr       (   nmodes) = kr(i)
               eigvector(:, nmodes) = eigvector(:, i)
            end if
        end do
        deallocate(cp)

    end subroutine

    subroutine GenerateModes(Layers,Ns,nmodes,eigvector,dz,z,dep,psi)

        implicit none
        integer,                         intent(in)  :: Layers 
        integer,                         intent(in)  :: Ns(Layers)        
        integer,                         intent(in)  :: nmodes
        real(rkind),    dimension(LL,CC),intent(in)  :: dep        
        complex(rkind), dimension(:, :), intent(in)  :: eigvector          
        real(rkind),                     intent(in)  :: dz
        real(rkind),    allocatable,     intent(out) :: z(:)          
        complex(rkind), allocatable,     intent(out) :: psi(:, :)        
        real(rkind),    allocatable, dimension(:)    :: xt
        real(rkind),    allocatable, dimension(:)    :: zt
        real(rkind),    allocatable, dimension(:)    :: za       
        complex(rkind), allocatable, dimension(:, :) :: psit
        integer                                      :: i, j, n, m
        
        allocate(z(ceiling(dep(Layers,Ns(Layers)+1)/dz+1)))
        do j = 1, size(z)
            z(j) = (j - 1) * dz
        end do
        allocate(za(size(z)), psi(size(z),nmodes))
        
        n = 1
        m = 0
        do i = 1, Layers
            allocate(zt(ceiling((dep(i,Ns(i)+1)-dep(i,1))/dz+1)))
            allocate(xt(size(zt)), psit(size(zt),nmodes))
            do j = 1, size(zt)
                zt(j) = dep(i, 1) + (j - 1) * dz
                xt(j) = -2.0_rkind / (dep(i,Ns(i)+1)-dep(i,1)) * zt(j) + &
                                     (dep(i,Ns(i)+1)+dep(i,1)) / (dep(i,Ns(i)+1)-dep(i,1))
            end do
            psit = InvChebMatrix(eigvector(n:n+Ns(i),1:nmodes), xt)
        
            if(i < Layers .and. zt(size(zt)) >= dep(i+1,1)) then
                za (m+1:m+size(zt)-1)  =  zt(1:size(zt)-1)
                psi(m+1:m+size(zt)-1,:)=  psit(1:size(zt)-1,:)
                m = m+size(zt)-1
            else                        
                za (m+1:m+size(zt))   = zt
                psi(m+1:m+size(zt),:) = psit
                m = m+size(zt)
            endif
            n = n + Ns(i) + 1
            deallocate(zt, xt, psit)
        end do
      
        do i = 1, nmodes
            call Interpolation_SSP(za,psi(:,i),z,psi(:,i))
        end do       

        deallocate(za)

    end subroutine

    subroutine Normalization(Layers,Ns,nmodes,psi,eigvector,rho,rhoh,kr,kh,Lowerboundary,dep)

        implicit none
        integer,                         intent(in)   :: Layers 
        integer,        dimension(:),    intent(in)   :: Ns
        integer,                         intent(in)   :: nmodes
        complex(rkind), dimension(:, :), intent(inout):: psi
        complex(rkind), dimension(:, :), intent(in)   :: eigvector
        real(rkind),    dimension(:, :), intent(in)   :: rho
        real(rkind),                     intent(in)   :: rhoh        
        complex(rkind),                  intent(in)   :: kr(nmodes)
        complex(rkind),                  intent(in)   :: kh        
        character(len=1),                intent(in)   :: Lowerboundary
        real(rkind),    dimension(LL,CC),intent(in)   :: dep        
        real(rkind)   , allocatable, dimension(:)     :: x
        real(rkind)   , allocatable, dimension(:)     :: P
        complex(rkind), allocatable, dimension(:,:)   :: f        
        real(rkind)   , allocatable, dimension(:,:)   :: Cov
        complex(rkind), allocatable, dimension(:)     :: norm        
        integer                                       :: n, i, j        
           
        allocate(norm(nmodes))
        norm = 0.0_rkind
        n = 1
        do i = 1, Layers
            allocate(x(Ns(i)+1), Cov(Ns(i)+1,Ns(i)+1), P(Ns(i)+1), f(Ns(i)+1,nmodes))
            do j = 1, Ns(i) + 1
                x(j) = cos((j-1) * pi / Ns(i))        
            end do
            Cov = 0.0_rkind
            Cov = Convolution(Cheb(1.0_rkind / rho(i,1:Ns(i)+1), x))            
            P   = 0.0_rkind
            do j = 0, Ns(i)-1, 2
                P(j+1) = - 2.0_rkind / (j * j - 1.0_rkind)
            end do
            f = 0.0_rkind
            do j = 1, nmodes
                f(:,j) = matmul(Convolution(eigvector(n:n+Ns(i), j)), eigvector(n:n+Ns(i), j))
            end do
            
            norm = norm + matmul(matmul(P,Cov), f) * (dep(i, Ns(i)+1) - dep(i,1)) * 0.5_rkind 
        
            n = n + Ns(i) + 1
            deallocate(x, Cov, P, f)
        end do        

        if(Lowerboundary == 'A') then
            do j = 1, nmodes
                norm(j) = norm(j) + 0.5_rkind / rhoh * psi(size(psi, 1), j) ** 2 / sqrt(kr(j) ** 2 - kh ** 2)            
            end do
        end if
        
        do j = 1, nmodes           
            psi(:, j) = psi(:, j) / sqrt(norm(j))                
        end do      

    end subroutine

    subroutine SynthesizeSoundField(nmodes,nr,r,kr,rhozs,zs,dz,psi,tl)

        implicit none
        integer,                  intent(in)    :: nmodes
        integer,                  intent(in)    :: nr
        real(rkind),              intent(in)    :: r(nr)
        real(rkind),              intent(in)    :: rhozs
        real(rkind),              intent(in)    :: zs
        real(rkind),              intent(in)    :: dz
        complex(rkind),           intent(in)    :: kr(nmodes)        
        real(rkind), allocatable, intent(inout) :: tl(:, :)
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
            
        allocate(p(size(psi,1), nr), tl(size(psi,1), nr))
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

    subroutine Msg
        write(*,*)
        write(*,*)'****************************************************************'
        write(*,*)' NM-CT is a program of range-independent acoustic propagation.'
        write(*,*)
        write(*,*)' Contact: Houwang Tu, Ph.D Candidate'
        write(*,*)'          National University of Defense Technology'
        write(*,*)'          Changsha, 410073, Hunan province, China'
        write(*,*)
        write(*,*)' Email  : tuhouwang@nudt.edu.cn'
        write(*,*)' Date   : 06/30/2022'
        write(*,*)'****************************************************************'
        write(*,*)
        return
    end subroutine
    
end module
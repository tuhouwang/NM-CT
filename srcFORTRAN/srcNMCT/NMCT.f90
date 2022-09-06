! Ocean acoustic normal modes.

! Copyright (C) 2020 Houwang Tu
! -------------------------------------------------------------------------
! This program is free software: you can redistribute it and/or modify it |
! under the terms of the GNU General Public License as published by the   |
! Free Software Foundation, either version 3 of the License, or (at your  |
! option) any later version.                                              |
!                                                                         |
! This code is distributed in the hope that it will be useful, but without|
! any warranty; without even the implied warranty of merchantability or   |
! fitness for a particular purpose. See the GNU General Public License for|
! more details.                                                           |
!                                                                         |
! You should have received a copy of the GNU General Public License along |
! with this program. If not, see <http://www.gnu.org/licenses/>.          |
!                                                                         |
! Originally developed as part of the author's article (H.Tu, Y.Wang, Q.  |
! Lan et al., A Chebyshev-Tau spectral method for normal modes of         |
! underwater sound propagation with a layered marine environment, Journal |
! of Sound and Vibration, https://doi.org/10.1016/j.jsv.2020.115784) under|
! the supervision of Prof. Yongxian Wang, National University of Defense  |
! Technology, China.                                                      |
!																		  |
! This Matlab/Scilab style code computes the layered and range-independent|
! modal acoustic field using the Chebyshev-Tau spectral method based on   |
! the normal modes.                                                       |
! -------------------------------------------------------------------------
! 11/08/2021                                                              |
!                                                                         |
! Roberto Sabatini et al proposed a new technique that can accurately     |
! solve acoustic semi-infinite space conditions based on the Chebyshev    |
! Collocation method in their paper "A multi-domain collocation method for| 
! the accurate computation of normal modes in open oceanic and atmospheric| 
! waveguides", Acta Acustica United with Acustica,105, 464–474 (2019)    |
! doi:10.3813/AAA.919328.                                                 |
!                                                                         |
! After carefully studying the technique, I integrated it into the NM-CT  |
! program and achieved good results. The latest NM-CT can accurately      |
! handle the perfectly free, perfectly rigid and infinite half-space      |
! conditions of the lower boundary! Simply specify the type of lower      |
! boundary condition in the input file.                                   |
!--------------------------------------------------------------------------
! 06/26/2022                                                              |
!                                                                         |
! New version of the program can calculate underwater acoustic propagation|
! in arbitrary horizontally stratified media. This improvement absorbs the|
! idea of domain-decomposition in the author's article (H.Tu, Y. Wang,    |
! Q. Lan et al., Applying a Legendre collocation method based on domain   |
! decomposition to calculate underwater sound propagation in a            |
! horizontally stratified environment, Journal of Sound and Vibration,    |
! https://doi.org/10.1016/j.jsv.2021.116364), and extends the range of    |
! solvable problems to the media with any number of layers.               |
!--------------------------------------------------------------------------

program NMCT
    use param_mod
    use nmct_mod
    implicit none
    
    external zgeev
    external zggev    
    external zgesv

    !---------------------------------------------------------------------------------
    ! Declare the variable needed later.
    character(len=MAX_FILENAME_LEN)              :: casename
    character(len=MAX_FILENAME_LEN)              :: data_file = "input.txt"
    character(len=MAX_FILENAME_LEN)              :: filename  = "tl.bin"
    character(len=1)                             :: Lowerboundary  
    integer(rkind)                               :: nr
    integer(rkind)                               :: nmodes
    real(rkind)                                  :: cpmax
    real(rkind)                                  :: dr
    real(rkind)                                  :: zs
    real(rkind)                                  :: zr
    real(rkind)                                  :: rmax
    real(rkind)                                  :: freq
    real(rkind)                                  :: dz
    real(rkind)                                  :: tlmin
    real(rkind)                                  :: tlmax
    real(rkind)                                  :: rhozs    
    real(rkind)                                  :: ch, rhoh, alphah
    complex(rkind)                               :: kh    
    integer(rkind)                               :: Layers    
    integer(rkind),allocatable, dimension(:)     :: Ns    
    real(rkind),   allocatable, dimension(:)     :: hinterface
    real(rkind),   allocatable, dimension(:)     :: r   
    real(rkind),                dimension(LL,CC) :: dep
    real(rkind),                dimension(LL,CC) :: c        
    real(rkind),                dimension(LL,CC) :: rho    
    real(rkind),                dimension(LL,CC) :: alpha
    complex(rkind),             dimension(LL,CC) :: ki
    complex(rkind),allocatable, dimension(:)     :: kr   
    complex(rkind),allocatable, dimension(:, :)  :: eigvector
    complex(rkind),allocatable, dimension(:, :)  :: psi
    real(rkind),   allocatable, dimension(:, :)  :: tl
    real(rkind),   allocatable, dimension(:)     :: z
    
    call Msg

    call ReadEnvParameter(casename,Layers,Ns,cpmax,freq,zs,zr,rmax,dr,hinterface,dz,&
            tlmin,tlmax,dep,c,rho,alpha,ch,rhoh,alphah,Lowerboundary,data_file)
                   
    call Initialization(Layers,Ns,freq,rmax,dr,zs,nr,r,rhozs,hinterface,ki,&
                        dep,c,rho,alpha,ch,alphah,kh)
                
    call EigenValueVector(Layers,Ns,dep,ki,rho,kh,rhoh,Lowerboundary,kr,eigvector)

    call NumofModes(Layers,freq,cpmax,kr,eigvector,nmodes)
           
    call GenerateModes(Layers,Ns,nmodes,eigvector,dz,z,dep,psi)
    
    call Normalization(Layers,Ns,nmodes,psi,eigvector,rho,rhoh,kr,kh,Lowerboundary,dep)
  
    call SynthesizeSoundField(nmodes,nr,r,kr,rhozs,zs,dz,psi,tl)        
        
    call SaveSoundField(filename,casename,tlmin,tlmax,r,z,tl)
    
    deallocate(Ns,hinterface,r,kr,eigvector,psi,z,tl)

end program

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

program NMCT
    use param_mod
    use nmct_mod
    implicit none
    
    external zgeev
	external zgesv
    !---------------------------------------------------------------------------------
    ! Declare the variable needed later.
    character(len=MAX_FILENAME_LEN)             :: casename
    character(len=MAX_FILENAME_LEN)             :: data_file = "input.txt"
    character(len=MAX_FILENAME_LEN)             :: filename  = "tl.bin"
    integer                                     :: Nw, Nb
    integer                                     :: Lowerboundary
    integer                                     :: nr
    integer                                     :: nmodes
    real(rkind)                                 :: cpmax
    real(rkind)                                 :: dr
    real(rkind)                                 :: zs
    real(rkind)                                 :: zr
    real(rkind)                                 :: rmax
    real(rkind)                                 :: freq
    real(rkind)                                 :: hinterface
    real(rkind)                                 :: Hb
    real(rkind)                                 :: dz
    real(rkind)                                 :: tlmin
    real(rkind)                                 :: tlmax
    real(rkind)                                 :: rhozs
    real(rkind),   allocatable, dimension(:)    :: rhow, alphaw, cw
    real(rkind),   allocatable, dimension(:)    :: rhob, alphab, cb
    real(rkind),   allocatable, dimension(:)    :: r
    complex(rkind),allocatable, dimension(:)    :: kw, kb
    complex(rkind),allocatable, dimension(:)    :: kr
    complex(rkind),allocatable, dimension(:, :) :: eigvectorw
    complex(rkind),allocatable, dimension(:, :) :: eigvectorb
    complex(rkind),allocatable, dimension(:, :) :: psi
    real(rkind),   allocatable, dimension(:)    :: z
    real(rkind),   allocatable, dimension(:, :) :: tl

    call ReadEnvParameter(casename,Nw,Nb,cpmax,freq,zs,zr,rmax,dr,hinterface,Hb,&
                    dz,Lowerboundary,tlmin,tlmax,rhow,rhob,alphaw,alphab,cw,cb,data_file)

    call Initialization(Nw,Nb,freq,rmax,dr,zs,rhow,rhob,cw,cb,alphaw,&
                    alphab,hinterface,Hb,nr,r,rhozs,kw,kb)

    allocate(kr(Nw+Nb-2),eigvectorw(Nw+1,Nw+Nb-2),eigvectorb(Nb+1,Nw+Nb-2))
    call EigenValueVector(Nw,Nb,hinterface,Hb,kw,kb,rhow,rhob,&
                    Lowerboundary,kr,eigvectorw,eigvectorb)

    call NumofModes(freq,kr,nmodes,cpmax)
    
    call GenerateModes(nmodes,dz,hinterface,Hb,eigvectorw,eigvectorb,psi,z)
    
    call Normalization(eigvectorw,eigvectorb,rhow,rhob,hinterface,Hb,Nw,Nb,nmodes,psi)
    
    call SynthesizeSoundField(nmodes,nr,r,kr,rhozs,zs,dz,psi,tl)        
        
    call SaveSoundField(filename,tlmin,tlmax,r,z,tl)
    
    deallocate(rhow,rhob,alphaw,alphab,cw,cb,r,kw,kb,kr,eigvectorw,eigvectorb,psi,z,tl)    

end program NMCT

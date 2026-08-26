       implicit none
       integer mxwork
       complex*16,allocatable:: v(:,:,:),coeff(:),wk(:)
       real,allocatable:: vp(:,:,:)
       real dx,dy,dz
       real dkx,dky,dkz
       integer ip,nxe,nye,nze,nkx,nky,nkz,np
       character*80  file1,file2
       character*80  file3,file4
       character*80  file5,file6
       mxwork=1000000
! 
        open(11,file="3d_full_fft.dat")
        read(11,*) nxe,nye,nze
        read(11,*) nkx,nky,nkz
        read(11,*) dx,dy,dz
        read(11,11) file1 
        read(11,11) file2 
        read(11,11) file3 
        read(11,11) file4 
        read(11,11) file5 
        read(11,11) file6 
        close(11) 
11      format(80a) 
       np=(2*nkx+1)*(2*nky+1)*(2*nkz+1)
       allocate(v(nxe,nye,nze),vp(nxe,nye,nze)) 
       allocate(coeff(np),wk(0:np))
        call fft_order(vp,v,wk,coeff,file1,nxe,nye,nze,np,nkx,nky,nkz,dx,dy,dz,dkx,dky,dkz) 
        open(11,file=trim(file2))
        write(11,*) nxe,nye,nze
        write(11,*) nkx,nky,nkz
        write(11,*) dx,dy,dz
        write(11,*) dkx,dky,dkz
        write(11,*) np
        do ip=1,np
        write(11,111) ip,coeff(ip)
        enddo
        close(11)
        call fft_order(vp,v,wk,coeff,file3,nxe,nye,nze,np,nkx,nky,nkz,dx,dy,dz,dkx,dky,dkz) 
        open(11,file=trim(file4))
        write(11,*) nxe,nye,nze
        write(11,*) nkx,nky,nkz
        write(11,*) dx,dy,dz
        write(11,*) dkx,dky,dkz
        write(11,*) np
        do ip=1,np
        write(11,111) ip,coeff(ip)
        enddo
        close(11)
        call fft_order(vp,v,wk,coeff,file5,nxe,nye,nze,np,nkx,nky,nkz,dx,dy,dz,dkx,dky,dkz) 
        open(11,file=trim(file6))
        write(11,*) nxe,nye,nze
        write(11,*) nkx,nky,nkz
        write(11,*) dx,dy,dz
        write(11,*) dkx,dky,dkz
        write(11,*) np
        do ip=1,np
        write(11,111) ip,coeff(ip)
        enddo
        close(11)
111     format(i10,1x,2(1x,e14.7))
      stop
      end
!
        subroutine fft_order(vp,v,wk,coeff,fname,nxe,nye,nze,nkxz,nkx,nky,nkz,dx,dy,dz,dkx,dky,dkz) 
        integer ix,iy,iz,nxe,nye,nze,nhx,nhy,nhz,nkx,nky,nkz,ikx,iky,ikz,ikxz,nkxz,irec
        real vp(nxe,nye,nze)
        real xmax,ymax,zmax,dkx,dky,dkz,dx,dy,dz,pi 
        real*8  s
        complex*16 v(nxe,nye,nze),coeff(*),wk(0:nkxz-1)
        character*80 fname
        pi=3.141592654
       xmax=(nxe-1)*dx
       ymax=(nye-1)*dy
       zmax=(nze-1)*dz
       dkx=1./xmax
       dky=1./ymax
       dkz=1./zmax
        v=0.0 
        open(11,file=trim(fname),recl=4*nze,access="direct",form="unformatted")
        irec=0
        do iy=1,nye
        do ix=1,nxe
        irec=irec+1
        read(11,rec=irec) (vp(ix,iy,iz),iz=1,nze)
        do iz=1,nze
        v(ix,iy,iz)=dcmplx(vp(ix,iy,iz),0.d0)
        enddo
        enddo
        enddo
        close(11) 
        print *," i am here"
         s=0.0
         do ix=1,nxe
         do iy=1,nye
         do iz=1,nze
         s=s+vp(ix,iy,iz)
         enddo
         enddo
         enddo
         s=s/nxe/nye/nze
         print *,"s=",s
!
! do fft in z direction
!
       print *,"z fft"
       do ix=1,nxe
       do iy=1,nye
       do iz=1,nze
       wk(iz-1)=v(ix,iy,iz)
       enddo
       call afft1d(wk,nze,1)
       do iz=1,nze
       v(ix,iy,iz)=wk(iz-1)
       enddo
       call reorder(v(ix,iy,:),nze)
       enddo
       enddo
!
       print *,"y fft"
       do ix=1,nxe
       do iz=1,nze
       do iy=1,nye
       wk(iy-1)=v(ix,iy,iz)
       enddo
       call afft1d(wk,nye,1)
       do iy=1,nye
       v(ix,iy,iz)=wk(iy-1)
       enddo
       call reorder(v(ix,:,iz),nye)
       enddo
       enddo
!
       print *,"x fft"
       do iy=1,nxe
       do iz=1,nze
       do ix=1,nye
       wk(ix-1)=v(ix,iy,iz)
       enddo
       call afft1d(wk,nxe,1)
       do ix=1,nye
       v(ix,iy,iz)=wk(ix-1)
       enddo
       call reorder(v(:,iy,iz),nxe)
       enddo
       enddo
!       v=v*dx*dy*dz
777    format(3i10,1x,2f10.3)
!
      if(mod(nxe,2).eq.0)  nhx=nxe/2
      if(mod(nye,2).eq.0)  nhy=nye/2
      if(mod(nze,2).eq.0)  nhz=nze/2
      if(mod(nxe,2).eq.1)  nhx=(nxe+1)/2
      if(mod(nye,2).eq.1)  nhy=(nye+1)/2
      if(mod(nze,2).eq.1)  nhz=(nze+1)/2
     
      ikxz=0
      do  ikz=-nkz,nkz
      do  iky=-nky,nky
      do  ikx=-nkx,nkx
      ikxz=ikxz+1 
      coeff(ikxz)=v(nhx+ikx,nhy+iky,nhz+ikz)
      enddo
      enddo
      enddo
      return
      end
      subroutine reorder(cv,nx)
      implicit none
      integer,intent(in) :: nx
      complex*16,intent(inout) :: cv(nx)
      complex*16 :: tmp(nx)
      integer n1,n2
      if(mod(nx,2).eq.0) then !! even nx reorder
        n1=nx/2+1
        n2=nx/2-1
        tmp(1:n2)=cv(n1+1:nx)
        tmp(n2+1:nx)=cv(1:n1)
        cv(1:nx)=tmp(1:nx)
      else !! odd nx reorder
        n1=(nx+1)/2
        tmp(1:n1-1)=cv(n1+1:nx)
        tmp(n1:nx)=cv(1:n1)
        cv(1:nx)=tmp(1:nx)
      endif
      return
      end subroutine
subroutine afft1d(a, n, isign)
  complex(16), intent(inout) :: a(0:n-1)
  integer,     intent(in)    :: n, isign

  integer,  parameter :: dp = kind(1.0d0)
  real(8), parameter :: pi = 3.141592653589793238462643383279502884_dp
  integer                  :: m, j, k
  real(8)                 :: s, ang
  complex(16), allocatable :: u(:), v(:), chirp(:)

  if (n <= 1) return

  !--- if N is already a power of two, nothing clever is needed -----------
  m=n
  do while (mod(m,2) == 0)
     m = m/2
  end do
  if (m == 1) then
     call fft1d(a, n, isign)
     return
  end if

  !--- smallest power of two that is at least 2N-1 ------------------------
  m = 1
  do while (m < 2*n - 1)
     m = 2*m
  end do

  allocate(u(0:m-1), v(0:m-1), chirp(0:n-1))

  !  the chirp.  The inverse carries the opposite sign, and its 1/N is
  !  applied at the very end.
  s = -1.0_dp
  if (isign > 0) s = 1.0_dp

  do j = 0, n-1
     ang      = s * pi * real(j, dp) * real(j, dp) / real(n, dp)
     chirp(j) = cmplx(cos(ang), sin(ang), dp)
  end do

  u = (0.0_dp, 0.0_dp)
  v = (0.0_dp, 0.0_dp)

  do j = 0, n-1
     u(j) = a(j) * chirp(j)
     v(j) = conjg(chirp(j))
  end do
  do j = 1, n-1
     v(m-j) = conjg(chirp(j))
  end do

  !  circular convolution through the power-of-two transform
  call fft1d(u, m, -1)
  call fft1d(v, m, -1)
  u = u * v
  call fft1d(u, m, +1)          ! this already divides by M

  do k = 0, n-1
     a(k) = chirp(k) * u(k)
  end do

  !--- NORMALISE ---------------------------------------------------------
  if (isign > 0) a = a / real(n, dp)

  deallocate(u, v, chirp)

end subroutine afft1d
subroutine fft1d(a, n, isign)
  complex*16, intent(inout) :: a(0:n-1)
  integer,     intent(in)    :: n, isign
  integer,  parameter :: dp = kind(1.0d0)
  real*8, parameter :: pi = 3.141592653589793238462643383279502884
  integer     :: i, j, m, mmax, istep, k
  real(dp)    :: theta, s
  complex(dp) :: w, temp

  if (n <= 1) return

  !--- N must be a power of two -----------------------------------------
  m = n
  do while (mod(m,2) == 0)
     m = m/2
  end do
  if (m /= 1) then
     write(*,*) 'fft1d: N is not a power of two.  N = ', n
     write(*,*) '       use fft1d_any instead.'
     stop
  end if

  !--- bit reversal ------------------------------------------------------
  j = 0
  do i = 0, n-2
     if (i < j) then
        temp = a(i);  a(i) = a(j);  a(j) = temp
     end if
     m = n/2
     do while (m >= 1 .and. j >= m)
        j = j - m
        m = m/2
     end do
     j = j + m
  end do

  !--- butterflies -------------------------------------------------------
  s = real(isign, dp)
  mmax = 1
  do while (mmax < n)
     istep = 2*mmax
     do k = 0, mmax-1
        theta = s * pi * real(k, dp) / real(mmax, dp)
        w = cmplx(cos(theta), sin(theta), dp)
        do i = k, n-1, istep
           j    = i + mmax
           temp = w * a(j)
           a(j) = a(i) - temp
           a(i) = a(i) + temp
        end do
     end do
     mmax = istep
  end do

  !--- NORMALISE ---------------------------------------------------------
  if (isign > 0) a = a / real(n, dp)

end subroutine fft1d

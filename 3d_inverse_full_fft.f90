       implicit none
       complex*16,allocatable:: v(:,:,:),wk(:)
       complex*16,allocatable:: trace(:,:),coeff(:) 
       complex dummy
       real*8 s
       real dx,dy,dz,dkx,dky,dkz,xmax,ymax,zmax,pi,fmax,df
       integer ip,ix,iy,iz,ikx,iky,ikz,nxe,nye,nze,nkx,nky,nkz,np,ikxz,nhx,nhy,nhz,irec
       integer ifreq,nfreq 
       character*80  file1,file2
       character*80  file3,file4
       character*80  file5,file6
       character*5 ifreq5,ishot5
! 
        call system("rm -rf snap3d")
        call system("rm -rf real_snap3d")
        call system("mkdir  snap3d")
        call system("mkdir  real_snap3d")
        open(11,file="vp_coeff.asc")
        read(11,*) nxe,nye,nze
        read(11,*) nkx,nky,nkz
        read(11,*) dkx,dky,dkz
        read(11,*) dx,dy,dz
        read(11,*) np
        close(11) 
11   format(a80)
       open(11,file="conect.dat")
       read(11,*) fmax,nfreq,df
       close(11)
       allocate(v(nxe,nye,nze),trace(nxe,nfreq)) 
       allocate(coeff(np),wk(0:np-1))
      if(mod(nxe,2).eq.0)  nhx=nxe/2
      if(mod(nye,2).eq.0)  nhy=nye/2
      if(mod(nze,2).eq.0)  nhz=nze/2
      if(mod(nxe,2).eq.1)  nhx=(nxe+1)/2
      if(mod(nye,2).eq.1)  nhy=(nye+1)/2
      if(mod(nze,2).eq.1)  nhz=(nze+1)/2
           do ifreq=1,nfreq
           write(*,*) ifreq
           write(ifreq5,"(i5.5)") ifreq
           write(ishot5,"(i5.5)") 1
           open(11,file="./true_wavefield/true."//ifreq5//"."//ishot5,form="formatted")
           do ip=1,np
           read(11,*) ix,coeff(ip)
           enddo
           close(11)
55   format(i10,4(1x,e14.7))
      v=0.0
      ikxz=0
      do ikz=nkz,-nkz,-1
      do iky=nky,-nky,-1
      do ikx=nkx,-nkx,-1
      ikxz=ikxz+1 
!      print *,ikx,iky,ikz,ikxz,coeff(ikxz)  
      v(nhx+ikx,nhy+iky,nhz+ikz)=coeff(ikxz)
      enddo
      enddo
      enddo
!
     print *,"i am here"
     call fftorder3d(v,nxe,nye,nze)
     print *,"i am here"
     v=dconjg(v)
     print *,"i am here"
!
! do fft in z direction
!
       print *,"z fft"
       do ix=1,nxe
       do iy=1,nye
       do iz=1,nze
       wk(iz-1)=v(ix,iy,iz)
       enddo
       call afft1d(wk,nze,-1)
       do iz=1,nze
       v(ix,iy,iz)=wk(iz-1)
       enddo
       enddo
       enddo
!
       print *,"y fft"
       do ix=1,nxe
       do iz=1,nze
       do iy=1,nye
       wk(iy-1)=v(ix,iy,iz)
       enddo
       call afft1d(wk,nye,-1)
       do iy=1,nye
       v(ix,iy,iz)=wk(iy-1)
       enddo
       enddo
       enddo
!
       print *,"x fft"
       do iy=1,nxe
       do iz=1,nze
       do ix=1,nxe
       wk(ix-1)=v(ix,iy,iz)
       enddo
       call afft1d(wk,nxe,-1)
       do ix=1,nxe
       v(ix,iy,iz)=wk(ix-1)
       enddo
       enddo
       enddo
!
       open(11,file="./snap3d/freq_snap."//ifreq5//"."//ishot5,recl=8*nze,access="direct",form="unformatted")
       irec=0
       do iy=1,nye
       do ix=1,nxe
       irec=irec+1
       write(11,rec=irec) (cmplx(v(ix,iy,iz)),iz=1,nze)
       enddo
       enddo
       close(11)  
 
       open(11,file="./real_snap3d/freq_snap."//ifreq5//"."//ishot5,recl=4*nze,access="direct",form="unformatted")
       irec=0
       do iy=1,nye
       do ix=1,nxe
       irec=irec+1
       write(11,rec=irec) (real(cmplx(v(ix,iy,iz))),iz=1,nze)
       enddo
       enddo
       close(11)  
       do ix=1,nxe
       trace(ix,ifreq)=v(ix,nye/2,15)
       enddo
!
       enddo ! ifreq 
!
       open(11,file="freq_sponge_trace.asc",form="formatted")
       do ix=1,nxe
       do ifreq=1,nfreq
       write(11,*)ix,ifreq+1,trace(ix,ifreq) 
       enddo
       enddo
       close(11)
       stop
       end
      subroutine fftorder(nx,cv)
      implicit none
      integer,intent(in) :: nx
      complex*16,intent(inout) :: cv(1)
      complex*16 :: tmp(10000)
      integer n1,n2
      tmp=0.0 
      if(mod(nx,2).eq.0) then !! even nx reorder
        n1=nx/2+1
        n2=nx/2-1
        tmp(1:n2)=cv(n1+1:nx)
        tmp(n2+1:nx)=cv(1:n1)
        cv(1:nx)=tmp(1:nx)
      else !! odd nx reorder
        n1=(nx+1)/2
        tmp(1:n1)=cv(n1:nx)
        tmp(n1+1:nx)=cv(1:n1-1)
        cv(1:nx)=tmp(1:nx)
      endif
      return
      end subroutine

      subroutine fftorder3d(green,nx,ny,nz)
      implicit none
      integer,intent(in) :: nx,ny,nz
      complex*16,intent(inout) :: green(nx,ny,nz)
      integer ix,iy,iz
      do ix=1,nx
      do iy=1,ny
        call fftorder(nz,green(ix,iy,:))
      enddo
      enddo
      do ix=1,nx
      do iz=1,nz
        call fftorder(ny,green(ix,:,iz))
      enddo
      enddo
      do iy=1,ny
      do iz=1,nz
        call fftorder(nx,green(:,iy,iz))
      enddo
      enddo
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

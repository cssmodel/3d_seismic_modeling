!
!      freq to time inverse fourier transform
!
        complex, allocatable ::snap(:,:,:,:)
        complex*16, allocatable ::temp(:),wk(:)
        character*6 ifreq6,idamp6
        character*5 ishot5,ifreq5
        character*250 fname
        character*10 ifreq10,idamp10
       open(11,file="conect.dat")  
       read(11,*)  fmax,nfreq,df
       close(11)
        dtnew=0.05
        nxe=256
        nye=256
        nze=64
        tmax=10.0
        dt=1./(2.*fmax)
        fmaxnew=1./(2.*dtnew)
        nnfreq=tmax/dtnew
        nnfreq=alog(float(nnfreq))/alog(2.0)+1
        nnfreq=2**nnfreq 
        print *,"nnfreq=",nnfreq
        allocate(snap(nxe,nye,nze,nnfreq),temp(0:nnfreq-1),wk(nnfreq))
        snap=0.0
        do ishot=1,1
        do ifreq=1,nfreq
        print *,"ifreq=",ifreq  
        write(ifreq5,"(i5.5)") ifreq
        write(ishot5,"(i5.5)") ishot
       open(11,file="./snap3d/freq_snap."//ifreq5//"."//ishot5,recl=8*nze,access="direct",form="unformatted")
       do ix=1,nxe
       do iy=1,nye
       irec=(ix-1)*nye+iy
       read(11,rec=irec) (snap(ix,iy,iz,ifreq),iz=1,nze)
       enddo
       enddo
       close(11)  
       enddo
       enddo
!
        do ix=1,nxe
        do iy=1,nye
        do iz=50,50
        print *,ix,iy,iz
        temp=dcmplx(0.0,0.0)
        temp(0)=0.0
        do kk=1,nfreq
        temp(kk)=dcmplx(snap(ix,iy,iz,kk))
        temp(nnfreq-kk)=conjg(temp(kk))
        enddo
!       temp=conjg(temp)
        call afft1d(temp,nnfreq,1)
        do kk=1,nnfreq
        t=(kk-1)*dtnew
        snap(ix,iy,iz,kk)=temp(kk-1)*exp(0.2*t)
        enddo
        enddo
        enddo
        enddo

!
        k=50
        open(14,file="time_trace_x_y_slice_z_18_sponge_20.bin",recl=4*nxe*nye,access="direct")
        do if=1,nnfreq
        write(14,rec=if) (((real(snap(i,j,k,if))),j=1,nye),i=1,nxe) 
        enddo
        close(14)
        open(14,file="time_trace.bin",recl=4*nnfreq,access="direct")
        do ix=1,nxe 
        write(14,rec=ix) (real(snap(ix,nye/2,k,ifreq)),ifreq=1,nnfreq) 
        enddo
        close(14)
        print *,"nnfreq=",nnfreq
        stop
        end
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

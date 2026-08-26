       implicit none 
       include 'mpif.h'
       complex,allocatable :: rhocoeff(:),vpcoeff(:),fxyzcoeff(:)
       complex*16,allocatable ::dampcoeff(:),source(:),wk(:)
       real*4,allocatable  :: wavlet(:)
       complex,allocatable :: crhs(:)
       complex,allocatable :: aa(:,:),cc(:,:),impedance(:,:),tmp_array(:,:)
       complex,allocatable :: kxx(:),kyy(:),kzz(:),kxx_kyy_kzz(:) 
       complex,allocatable :: kxkz(:),kzkx(:),kxky(:),kykx(:),kykz(:),kzky(:) 
       integer ipiv(1000000),info
       integer i,irec,nxe,nye,nze,nkx,nky,nkz,ip,jp,jloc,nh,np,npp,ierr,npe,ipe
       integer ifreq,nfreq
       integer is
       real*4 fmax,tmax,pi,dx,dy,dz,dkx,dky,dkz,df,dxx,dt,t1,t2,t3
       complex ci,omega
       character*5 ifreq5,ishot5
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,ipe,ierr)
      call mpi_comm_size(mpi_comm_world,npe,ierr)
       call system("rm -rf true_wavefield")
       call system("rm -rf snap3d")
       call system("rm -rf real_snap3d")
       call system("mkdir true_wavefield")
       call system("mkdir  snap3d")
       call system("mkdir real_snap3d")
       ci=cmplx(0.,1.0)
       pi=3.141592654
       open(11,file="vp_coeff.asc")
       read(11,*) nxe,nye,nze
       read(11,*) nkx,nky,nkz
       read(11,*) dx,dy,dz
       read(11,*) dkx,dky,dkz
       read(11,*) np 
       npp=2*np-1
allocate(aa(np,np),cc(np,np),impedance(np,np),tmp_array(np,npp),crhs(np),wk(np))  
allocate(dampcoeff(np),vpcoeff(np),rhocoeff(np),fxyzcoeff(np),kxx(np),kyy(np),kzz(np),kxx_kyy_kzz(np))
allocate(kxkz(np),kzkx(np),kxky(np),kykx(np),kykz(np),kzky(np)) 
       do i=1,np
       read(11,3333) irec,vpcoeff(irec)
       enddo
       close(11)
       open(11,file="damp_coeff.asc")
       read(11,*) nxe,nye,nze
       read(11,*) nkx,nky,nkz
       read(11,*) dx,dy,dz
       read(11,*) dkx,dky,dkz
       read(11,*) np 
       do i=1,np
       read(11,3333) irec,dampcoeff(irec)
       enddo
       close(11)
       open(11,file="fxyz_coeff.asc")
       read(11,*) nxe,nye,nze
       read(11,*) nkx,nky,nkz
       read(11,*) dx,dy,dz
       read(11,*) dkx,dky,dkz
       read(11,*) np 
       do i=1,np
       read(11,3333) irec,fxyzcoeff(irec)
       enddo
       close(11)
3333    format(i10,1x,2(1x,e14.7))
       npp=2*np-1
       dxx=(nxe*dx)/(2*nkx+1)/1.5
       fmax=1.5/(2.*dxx)  
       dt=1./(2.*fmax)
       tmax=10.0
       df=1./tmax
       nfreq=fmax/df+1
       nfreq=alog(float(nfreq))/alog(2.0)+1
       nfreq=2**nfreq
       open(11,file="conect.dat")
       write(11,*) fmax,nfreq,df
       close(11)
        cc=0.0
        nh=np/2
        do ip=1,np
        do jp=1,np
        jloc=ip+jp-1
        tmp_array(ip,jloc)=dampcoeff(jp)
        enddo
        enddo
         nh=np/2
 
        do ip=1,np
        do jp=nh+1,npp-nh+1    
        cc(ip,jp-nh)=tmp_array(ip,jp)
        enddo
        enddo
! 
!
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
allocate(wavlet(2*nfreq),source(0:2*nfreq-1))
       call fdgaus(wavlet,fmax,dt,2*nfreq)
       source=dcmplx(wavlet,0.0)
       call afft1d(source,2*nfreq,1)
       cc=cc/(dt/10.)
       call forward_3d_transform(np,nkx,nky,nkz,kxx,kyy,kzz,kxky,kxkz,kykz,kxx_kyy_kzz,dkx,dky,dkz)
       aa=0.0
       call genConvMat(np,vpcoeff,kxx_kyy_kzz,aa)
do 1000 ifreq=1+ipe,nfreq,npe
       print *,"ifreq=",ifreq
       write(ifreq5,"(i5.5)") ifreq
       omega=cmplx(2.*pi*ifreq*df,0.1)
       impedance=0.0
       impedance=aa+ci*omega*cc
       do 1100 ip=1,np
       impedance(ip,ip)=impedance(ip,ip)+omega**2
1100   continue
       is=0
       do ip=1,np
       do jp=1,np
       if(cabs(impedance(ip,jp)).eq.0.0) then
       is=is+1
       endif
       enddo     
       enddo     
!
call cpu_time(t1)
call cgetrf(np,np,impedance,np,ipiv,info)       
call cpu_time(t2)
write(ishot5,"(i5.5)") 1
crhs=0.0
crhs=fxyzcoeff*source(ifreq) 
call cgetrs('N',np,1,impedance,np,ipiv,crhs,np,info)
call cpu_time(t3)
print *,ifreq,"fact time=",t2-t1,"sol time=",t3-t2
open(11,file="./true_wavefield/true."//ifreq5//"."//ishot5,form="formatted")
do 1200 ip=1,np
write(11,*) ip,crhs(ip)
1200 continue
close(11)
1000 continue
!
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call mpi_finalize(ierr)
     stop
     end

      subroutine forward_3d_transform(np,nkx,nky,nkz,kxx,kyy,kzz,kxky,kxkz,kykz,kxx_kyy_kzz,dkx,dky,dkz)
      integer ikx,iky,ikz,nkx,nky,nkz,ikxz,np
      complex kxx(np),kyy(np),kzz(np),kxx_kyy_kzz(np)
      complex kxky(np),kxkz(np),kykz(np)
      real*4 pi,dkx,dky,dkz
      real*4 t1,t2
      complex ci
      pi=3.141592654
      ci=cmplx(0.0,1.0)
      call cpu_time(t1) 
           ikxz=0
           do ikz=-nkz,nkz
           do iky=-nky,nky
           do ikx=-nkx,nkx
           ikxz=ikxz+1
           kxx(ikxz)=cmplx(-(2.*pi*ikx*dkx)**2,0.0)
           kyy(ikxz)=cmplx(-(2.*pi*iky*dky)**2,0.0)
           kzz(ikxz)=cmplx(-(2.*pi*ikz*dkz)**2,0.0)
           kxky(ikxz)=cmplx(-(2.*pi*ikx*dkx*2.*pi*iky*dky),0.0)
           kxkz(ikxz)=cmplx(-(2.*pi*ikx*dkx*2.*pi*ikz*dkz),0.0)
           kykz(ikxz)=cmplx(-(2.*pi*iky*dky*2.*pi*ikz*dkz),0.0)
           kxx_kyy_kzz(ikxz)=kxx(ikxz)+kyy(ikxz)+kzz(ikxz)
           enddo
           enddo
           enddo
!
      call cpu_time(t2) 
      print *,"cpu_time=",t2-t1 ,"in forward transform"
      close(10)
!
      return
      end
          subroutine genConvMat(nx,vec1,vec2,mat)
          implicit none
          integer,intent(in) :: nx
          complex,intent(in) :: vec1(nx),vec2(nx)
          complex,intent(out):: mat(nx,nx)
          integer ilow,ihigh,i1,i2,icol,iz,irow,iv,ik
          mat(:,:)=cmplx(0.,0.)
          ilow=nx/2
          ihigh=nx/2+nx-1
       if(mod(nx,2).eq.1)then
          ilow=ilow+1
          ihigh=ihigh+1
       endif

          irow=0
          do iz=ilow,ihigh
            irow=irow+1

            i1=iz-(nx-1)
            i2=iz
            if(i1<1) i1=1
            if(i2>nx) i2=nx

            iv=i2+1
            do icol=i1,i2
              iv=iv-1
              ik=icol

              mat(irow,icol)=vec1(iv)*vec2(ik)
            enddo
          enddo
          end subroutine
      subroutine gauss(w,cutoff,dt,nt)
      integer i,icut,nt
      real w(nt)
      real a,cutoff,pi,amp,arg,dt,t,t0
      pi=3.141592654
      a=pi*(5.*cutoff/8.)**2
      amp=sqrt(a/pi)
      do 10 i=1,nt
      t=(i-1)*dt
      arg=-a*t**2
      if(arg.lt.-32.) arg=-32.
10    w(i)=amp*exp(arg)
      do 20 i=1,nt
      if(w(i).lt.0.001*w(1)) then
      icut=i
      t0=(icut-1)*dt
      go to 30
      endif
20    continue
30    continue
      do 40 i=1,nt
      t=(i-1)*dt
      t=t0-t
      arg=-a*t**2
      if(arg.lt.-32.) arg=-32.
40    w(i)=amp*exp(arg)
      return
      end
      subroutine fdgaus(w,cutoff,dt,nt)
      integer i,icut,nt
      real w(nt)
      real a,cutoff,pi,amp,arg,dt,t,t0,smax
      pi=4*atan(1.)
      a=pi*(5.*cutoff/8.)**2
      amp=sqrt(a/pi)
      do 10 i=1,nt
      t=(i-1)*dt
      arg=-a*t**2
      if(arg.lt.-32.) arg=-32.
10    w(i)=amp*exp(arg)
      do 20 i=1,nt
      if(w(i).lt.0.001*w(1)) then
      icut=i
      t0=(icut-1)*dt
      go to 30
      endif
20    continue
30    continue
      do 40 i=1,nt
      t=(i-1)*dt
      t=t-t0
      arg=-a*t**2
      if(arg.lt.-32.) arg=-32.
40    w(i)=-2.*sqrt(a)*a*t*exp(arg)/sqrt(pi)
      smax=0.
      do i=1,nt
      if(abs(w(i)).gt.smax) smax=abs(w(i)) 
      enddo
      do i=1,nt
      w(i)=w(i)/smax
      enddo
      return
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

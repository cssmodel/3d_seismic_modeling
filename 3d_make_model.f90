        parameter(nxe=256,nye=256,nze=64)
        real*4 vp(nxe,nye,nze),sponge(nxe,nye,nze),rho(nxe,nye,nze),dampxyz(nxe,nye,nze)
        real*4 fxyz(nxe,nye,nze)
!
    dt=1.00
    dampmax=0.98
    rho=1.0
        npml=10
    call assign_sponge(sponge,nxe,nye,nze,npml,dampmax,rho,dt,dampxyz)
        open(11,file="sponge.bin",recl=4*nze,access="direct",form="unformatted")
        irec=0
        do iy=1,nye
        do ix=1,nxe
        irec=irec+1
        write(11,rec=irec) (sponge(ix,iy,iz),iz=1,nze)
        enddo 
        enddo 
        close(11)
! 
        fxyz=0.0
        fxyz(128,128,20)=1.0 
        open(11,file="fxyz.bin",recl=4*nze,access="direct",form="unformatted")
        irec=0
        do iy=1,nye
        do ix=1,nxe
        irec=irec+1
        write(11,rec=irec) (fxyz(ix,iy,iz),iz=1,nze)
        enddo 
        enddo 
        close(11)
        stop
        end
        
      
    subroutine assign_sponge(sponge,nx,ny,nz,npml,dampmax,rho,dt,dampxyz)
    implicit none
    integer i,j,k 
    integer nx,ny,nz,npml
    real dampmax,rho(nx,ny,nz),dt
    real sponge(nx,ny,nz)
    real dampx(2*nx),dampy(2*ny),dampz(2*nz),dampxyz(nx,ny,nz)
    dampx=1.0
    dampy=1.0
    dampz=1.0
   
    print *," i am here"
    j=0
    do i=nx-npml+1,nx
        j=j+1
        dampx(i)= 1.0 - float(j)/npml * (1.0-dampmax)
        dampx(npml-j+1)=dampx(i)
    enddo
    j=0
    do i=ny-npml+1,ny
        j=j+1
        dampy(i)= 1.0 - float(j)/npml * (1.0-dampmax)
        dampy(npml-j+1)=dampy(i)
    enddo
    j=0
    do i=nz-npml+1,nz
        j=j+1
        dampz(i)= 1.0 - float(j)/npml * (1.0-dampmax)
        dampz(npml-j+1)=dampz(i)
    enddo
    open(22,file='dampx.out')
    do i=1,nx
        write(22,*) i, dampx(i)
    enddo
    close(22)
    open(22,file='dampy.out')
    do i=1,ny
        write(22,*) i, dampy(i)
    enddo
    close(22)
    open(22,file='dampz.out')
    do i=1,nz
        write(22,*) i, dampz(i)
    enddo
    close(22)
    do i=1,nx
      do j=1,ny
        do k=1,nz
            dampxyz(i,j,k)=dampx(i)*dampy(j)*dampz(k)
          enddo
        enddo
    enddo
    sponge = 2*(1-dampxyz)/dampxyz * rho**3/dt
    end subroutine
  

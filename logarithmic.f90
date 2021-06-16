program logarithmic

real::avmax,avmin
real::lambda,n
real::x,av,y,z,size,yzdim
integer::yzres,i,j,jj,k,tot,resolution
real::av0,cm,length
character(len=4)::choice
real,allocatable::r(:,:),rcut(:,:)
!************** V
!real::av_1,av_2,displacement
!integer::NN
!************** ^

cm=3.0856e18
av0=6.289e-22
write(6,*) 'give number density [cm^-3]'
read(5,*) n
write(6,*) 'give avmax'
read(5,*) avmax
length=avmax/cm/av0/n
write(6,*) 'length [pc] = ',length
write(6,*) 'give log10(avmin)'
read(5,*) avmin
10 continue
write(6,*) 'logarithmic [loga], uniform [unif], unilog [unlo] ?'
read(5,*) choice
write(6,*) 'give x-resolution'
read(5,*) resolution
if (choice.eq.'unif') write(6,*) 'displacement = ',avmax/real(resolution)
write(6,*) 'give yz-resolution'
read(5,*) yzres

tot=0
open(unit=1,file='outgrid.dat',status='replace')
write(1,'(4E15.7)') 0.,0.,0.,n
open(unit=2,file='av.dat',status='replace')
open(unit=3,file='cut.dat',status='replace')

if (yzres.eq.0) then
!one dimensional line
  if (choice.eq.'loga') then
!     lambda=(log10(avmax)+5.0)*real(resolution)
     lambda=(log10(avmax)-avmin)*real(resolution)
     do i=0,int(lambda)
!       av=10.0**(-5.+real(i)/real(resolution))
       av=10.0**(avmin+real(i)/real(resolution))
       x=av/cm/av0/n
       write(1,'(4E15.7)') x,0.,0.,n
       write(2,'(2E15.7)') x,av
       tot=tot+1
     enddo
  else if (choice.eq.'unif') then
     lambda=avmax*real(resolution)
     do i=0,int(lambda)
       av=real(i)/real(resolution)
       x=av/cm/av0/n
       tot=tot+1
       write(1,'(4E15.7)') x,0.,0.,n
       write(2,'(2E15.7)') x,av
     enddo
  else if (choice.eq.'unlo') then
     do i=int(-avmin),0,-1!0,3
       do j=1,resolution
        ! av=real(j)*(real(avmax)/10**(i))/real(resolution)
!        av=(real(j)/real(resolution))*real(avmax)/10**(i)
        av=(real(j)/real(resolution))*real(avmax)/10**(i)
        x=av/cm/av0/n
        tot=tot+1
        write(1,'(4E15.7)') x,0.,0.,n
        write(2,'(2E15.7)') x,av
       enddo
     enddo

   else
    write(6,*) 'choice not acceptable';write(6,*) ''
    goto 10
  endif

write(6,*) 'total x-points = ',tot+1
else
!three dimensional cube
  size=2.*avmax/cm/av0/n!/4.
  yzdim=size/real(yzres)
  if (choice.eq.'loga') then
       if (.not.allocated(r)) allocate(r(1:3,1:1000000))
    lambda=(log10(avmax)+3.0)*real(resolution)
    do i=0,int(lambda)
      av=10.0**(-3.+real(i)/real(resolution))
      x=av/cm/av0/n
      write(2,'(2E15.7)') x,av
      do j=0,yzres
         y=-size/2. + real(j)*yzdim
         do k=0,yzres
            z=-size/2.+real(k)*yzdim
            write(1,'(4E15.7)') x,y,z,n
            tot=tot+1
            if (tot.ge.1000000) stop 'increase r array'
            r(1,tot)=z;r(3,tot)=-x+length;r(2,tot)=y
         enddo
      enddo
    enddo
  else if (choice.eq.'unif') then
    lambda=avmax*real(resolution)
    do i=0,int(lambda)
      av=real(i)*real(resolution)
      x=av/cm/av0/n
      write(2,'(2E15.7)') x,av
      do j=0,yzres
        y=-size/2. + real(j)*yzdim
        do k=0,yzres
           z=-size/2.+real(k)*yzdim
            write(1,'(4E15.7)') x,y,z,n
            tot=tot+1
        enddo
      enddo
    enddo
  else if (choice.eq.'unlo') then
  size=2.*avmax/cm/av0/n
  yzdim=size/real(yzres)
     do i=int(-avmin),0,-1
!************************* V
!       av_1 = (1./real(resolution))*real(avmax)/10**(i)
!       av_2 = (2./real(resolution))*real(avmax)/10**(i)
!       displacement = 100.*(av_2-av_1)/cm/av0/n
!       write(6,*) 'Displacement at AV ~ 10**(',i,') = ',displacement
!       NN = int(size/real(displacement))
       if (.not.allocated(r)) allocate(r(1:3,1:1000000))
!       write(6,*) 'Particles to distribute = ',NN,NN**2
!************************* ^
       do j=1,resolution
        av=(real(j)/real(resolution))*real(avmax)/10**(i)
        x=av/cm/av0/n
        write(2,'(2E15.7)') x,av
!*************************V
!        do jj=0,NN
!          y=-size/2. + real(jj)*displacement*100.
!          do k=0,NN
!            z=-size/2. + real(k)*displacement*100.
!            write(1,'(4E15.7)') x,y,z,n
!            tot=tot+1
!            if (tot.ge.1000000) stop 'increase r array'
!            r(1,tot)=z;r(3,tot)=-x+length;r(2,tot)=y
!          enddo ! k=0,NN
!        enddo ! jj=0,NN
!************************* -----
        do jj=0,yzres
          y=-size/2. + real(jj)*yzdim
          do k=0,yzres
             z=-size/2.+real(k)*yzdim
              write(1,'(4E15.7)') x,y,z,n
              tot=tot+1
              if (tot.ge.1000000) stop 'increase r array'

           r(1,tot)=z;r(3,tot)=-x+length;r(2,tot)=y
          enddo
        enddo
!************************* ^
       enddo
     enddo
 endif

     k=0
     allocate(rcut(1:3,1:1000000))
     do i=1,tot
      if((r(3,i)).ge.abs(r(1,i)).and.r(3,i).ge.abs(r(2,i))) then
!      write(3,*) r(1:3,i)
       k=k+1
       rcut(1,k)=r(1,i);rcut(2,k)=r(2,i);rcut(3,k)=r(3,i)
      endif 
     enddo

     open(unit=4,file='box.dat',status='replace')
     j=0
     j=j+1
     write(6,*) 'Side = ',j
     do i=1,k
       write(4,*) rcut(1,i),rcut(2,i),rcut(3,i),n
     enddo
     j=j+1
     write(6,*) 'Side = ',j
     do i=1,k
       write(4,*) -rcut(1,i),-rcut(2,i),-rcut(3,i),n
     enddo
     j=j+1
     write(6,*) 'Side = ',j
     do i=1,k
       write(4,*) rcut(1,i),rcut(3,i),rcut(2,i),n
     enddo
     j=j+1
     write(6,*) 'Side = ',j
     do i=1,k
       write(4,*) rcut(1,i),-rcut(3,i),-rcut(2,i),n
     enddo
     j=j+1
     write(6,*) 'Side = ',j
     do i=1,k
       write(4,*) rcut(3,i),rcut(2,i),rcut(1,i),n
     enddo
     j=j+1
     write(6,*) 'Side = ',j
     do i=1,k
       write(4,*) -rcut(3,i),-rcut(2,i),-rcut(1,i),n
     enddo

write(6,*) 'total points = ',tot
write(6,*) 'file [box.dat] written!'
endif

write(6,*) 'file [outgrid.dat] written!'
end program

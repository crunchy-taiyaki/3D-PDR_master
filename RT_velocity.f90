program one
implicit none
character(len=50)::filepdr,fileline,filetau,filepop,prefix,spec(1:17),fileout
integer::itot,i,id,vv,j
real::dummy
double precision,allocatable::x(:),av(:),jnu(:,:),tau(:,:)
double precision,allocatable::N(:),Tex(:,:),Tgas(:),rho(:),Bnu(:,:)
double precision,allocatable::pop(:,:,:),abun(:,:)
double precision,allocatable::rx(:),rav(:),rjnu(:,:),rtau(:,:)
double precision,allocatable::rTgas(:),rrho(:)
double precision,allocatable::rpop(:,:,:),rabun(:,:),tr_incr(:,:)
double precision::c,kb,pc2cm,hp,pi,T,mhp,freq0(1:17),g(1:17,1:2),freq(1:17)
double precision::t_r(1:17),mh(1:17),Ntot,Ntgas
double precision::uv,dtau(1:17),t_a_i(1:17),t_a_ip1(1:17),Nwrite(1:17)
double precision::A(1:17),phi(1:17),sigma(1:17),tau_test(1:17),vturb,frac(1:17)
!double precision,allocatable::tau_test(:,:)

write(6,*) 'prefix?'
read(5,*) prefix
!prefix='Tt2'
write(6,*) '!!! vturb=1km/s !!!'
vturb=1*1d5 !cm/s

call constants
call readfile


allocate(Tex(1:17,1:itot))
allocate(N(1:33))
allocate(Bnu(1:17,1:itot))
allocate(tr_incr(1:itot,1:17))
do i=1,itot
  do j=1,17
    !if (g(j,2)*pop(j,1,i)/pop(j,2,i)/g(j,1).eq.1) then 
!    if (abs(g(j,2)*pop(j,1,i)/pop(j,2,i)/g(j,1)-1).lt.1e-2.or.pop(j,2,i).eq.0) then  !condition to handle explosions
!       Tex(j,i)=0
!       Bnu(j,i)=0
!       write(6,*) 'warning'
!    else 
       Tex(j,i)=(hp*freq0(j)/kb)/log(g(j,2)*pop(j,1,i)/pop(j,2,i)/g(j,1))
       Bnu(j,i)=(2.*hp*freq0(j)**3/c**2)/(exp(hp*freq0(j)/kb/Tex(j,i))-1)
!    endif
  enddo
enddo
!write(6,*) spec(7),freq0(7),g(7,2),pop(7,1,itot-2),pop(7,2,itot-2),g(7,1),Tex(7,itot-2)

do vv=0,0!-100,100
tau_test=0
freq=freq0/(1.+real(vv)*1e5/c/10.)
do i=1,itot-1
  sigma=(freq0/c)*sqrt(kb*Tgas(i)/mh+vturb**2/2.)
  phi=1./sigma/sqrt(2.*pi)
  !phi=phi*exp(-((1+0/c)*freq-freq0)**2/sigma**2/2.)
  frac=0.5*((pop(:,1,i)+pop(:,1,i+1))*g(:,2)/g(:,1)-(pop(:,2,i)+pop(:,2,i+1)))
  tau_test=tau_test+phi*(A*c**2/8./pi/freq0**2)*frac*abs(x(i+1)-x(i))*pc2cm
  !tau_test=tau_test+(0.1d0*(freq/1000d9)**2)*rho(i)*1.402*mhp*abs(x(i+1)-x(i))*pc2cm
  tau(:,i)=tau_test
enddo

open(unit=121,file='CII.out.dat',status='replace')
open(unit=122,file='CI.out.dat',status='replace')
open(unit=123,file='OI.out.dat',status='replace')
open(unit=124,file='CO.out.dat',status='replace')
fileout='NH2_dW.'//trim(adjustl(prefix))//'.dat'
open(unit=111,file=fileout,status='replace')
N=0;Ntot=0;Ntgas=0
t_r=0
do i=1,itot-2
  dtau=tau(:,i+1)-tau(:,i)
  do j=1,17
    if (dtau(j).gt.1d10) then
      t_r(j)=Bnu(j,i)
    else if (dtau(j).gt.1d-6) then
      t_a_i(j)=Bnu(j,i)*((1-exp(-dtau(j)))/dtau(j)-exp(-dtau(j)))
      t_a_ip1(j)=Bnu(j,i+1)*(1.-(1.-exp(-dtau(j)))/dtau(j))
      t_r(j)=t_r(j)*exp(-dtau(j))+t_a_i(j)+t_a_ip1(j)
    else
      t_r(j)=t_r(j)*(1-dtau(j))+(Bnu(j,i)+Bnu(j,i+1))*dtau(j)/2.
    endif
    if (t_r(j).gt.tr_incr(i-1,j)) then
      tr_incr(i,j)=t_r(j)
    else
      tr_incr(i,j)=tr_incr(i-1,j)
    endif
  enddo
enddo

do i=1,itot-2
  N=N+0.5*(rho(i)*abun(:,i)+rho(i+1)*abun(:,i+1))*abs(x(i+1)-x(i))*pc2cm
  Ntot=Ntot+0.5*(rho(i)+rho(i+1))*abs(x(i+1)-x(i))*pc2cm
  Ntgas=Ntgas+0.5*(rho(i)*tgas(i)+rho(i+1)*tgas(i+1))*abs(x(i+1)-x(i))*pc2cm
  !write(121,'(100ES11.3)') x(i),N(11),Tex(1,i),tau(1,i),tr_incr(i,1)*c**2/2./kb/freq0(1)**2
  !write(122,'(100ES11.3)') x(i),N(25),Tex(2,i),tau(2,i),tr_incr(i,2)*c**2/2./kb/freq0(2)**2
  !write(123,'(100ES11.3)') x(i),N(30),Tex(5,i),tau(5,i),tr_incr(i,5)*c**2/2./kb/freq0(5)**2
  !write(124,'(100ES11.3)') Ntot,N(28),Tex(8,i),tau(8,i),tr_incr(i,8)*c**2/2./kb/freq0(8)**2,rho(i)
  !write(111,'(7ES11.3)') Ntot*6.3e-22,N(31),t_r(1)*c**2/2./kb/freq0(1)**2,t_r(2)*c**2/2./kb/freq0(2)**2,&
!  write(111,'(7ES11.3)') x(i),N(25),t_r(1)*c**2/2./kb/freq0(1)**2,t_r(2)*c**2/2./kb/freq0(2)**2,&
!          &t_r(5)*c**2/2./kb/freq0(5)**2,t_r(8)*c**2/2./kb/freq0(8)**2,tgas(i)
enddo
!write(6,*)spec(7),N(30),rho(itot-2),abun(30,itot-2),x(itot-2),pc2cm
write(11,*) real(vv)/10., t_r(8)*c**2/2./kb/freq0(8)**2
enddo
write(6,*) '<Tgas>=',NTgas/Ntot
write(6,*) 'Linewidth (CO) = ',2.*sqrt(2.*log(2.))*sqrt(kb*(NTgas/Ntot)/mh(8)+vturb**2/2.)/1d5,' [km/s]'
Nwrite(1)=N(11);Nwrite(2:4)=N(25);Nwrite(5:7)=N(30);Nwrite(8:17)=N(28)
write(6,*) 'Species || Column density || Tex || tau || tau_test || Tr '
do i=1,17
  write(6,'(A10,2X,5ES11.3)') spec(i),Nwrite(i),Tex(i,itot-2),tau(i,itot-2),tau_test(i),tr_incr(itot-2,i)*c**2/2./kb/freq0(i)**2
enddo
t_r=t_r*c**2/2./kb/freq0**2
close(1);open(unit=1,file='CO_sled.dat',status='replace')
do i=1,10
  write(1,*) i,t_r(i+7)/t_r(8)
enddo

contains
subroutine readfile
!write(6,*) 'give prefix'
!read(5,*) prefix
filepdr='./ALL_TESTS/old/'//trim(adjustl(prefix))//'/'//"pdr.fin"
fileline='./ALL_TESTS/old/'//trim(adjustl(prefix))//'/'//"line.fin"
filetau='./ALL_TESTS/old/'//trim(adjustl(prefix))//'/'//"opdp.fin"
filepop='./ALL_TESTS/old/'//trim(adjustl(prefix))//'/'//"spop.fin"
open(unit=1,file=fileline,status='old')
open(unit=2,file=filepdr,status='old')
open(unit=3,file=filepop,status='old')
open(unit=4,file=filetau,status='old')
itot=0
do 
  read(1,*,end=100) dummy
  itot=itot+1
enddo
100 continue
rewind(1)
allocate(x(1:itot),av(1:itot),jnu(1:17,1:itot),abun(1:33,1:itot))
allocate(Tgas(1:itot),rho(1:itot),pop(1:17,1:2,1:itot),tau(1:17,1:itot))
allocate(rx(1:itot),rav(1:itot),rjnu(1:17,1:itot),rabun(1:33,1:itot))
allocate(rTgas(1:itot),rrho(1:itot),rpop(1:17,1:2,1:itot),rtau(1:17,1:itot))
do i=1,itot
  read(1,*) id,rx(i),rav(i),rjnu(1:17,i)
  read(2,*) id,dummy,dummy,rTgas(i),dummy,dummy,rrho(i),uv,rabun(1:33,i)
  read(3,*) id,dummy,rpop(1,1,i),rpop(1,2,i),dummy,dummy,dummy,&
       &rpop(2,1,i),rpop(2,2,i),rpop(3,2,i),dummy,dummy,&
       &rpop(5,1,i),rpop(5,2,i),rpop(6,2,i),dummy,dummy,&
       &rpop(8,1,i),rpop(8,2,i),rpop(9,2,i),rpop(10,2,i),rpop(11,2,i),&
       &rpop(12,2,i),rpop(13,2,i),rpop(14,2,i),rpop(15,2,i),rpop(16,2,i),rpop(17,2,i)
  rpop(3,1,i)=rpop(2,1,i);  rpop(4,2,i)=rpop(3,2,i);  rpop(4,1,i)=rpop(2,2,i)
  rpop(6,1,i)=rpop(5,1,i);  rpop(7,2,i)=rpop(6,2,i);  rpop(7,1,i)=rpop(5,2,i)
  rpop(9,1,i)=rpop(8,2,i);  rpop(10,1,i)=rpop(9,2,i);  rpop(11,1,i)=rpop(10,2,i)
  rpop(12,1,i)=rpop(11,2,i);  rpop(13,1,i)=rpop(12,2,i);  rpop(14,1,i)=rpop(13,2,i)
  rpop(15,1,i)=rpop(14,2,i);  rpop(16,1,i)=rpop(15,2,i);  rpop(17,1,i)=rpop(16,2,i)
  read(4,*) id,dummy,rtau(1:17,i)
enddo

x=rx;av=rav;jnu=rjnu;abun=rabun;Tgas=rTgas;rho=rrho;pop=rpop;tau=rtau

!do i=1,itot
!  j=itot-i
!  x(i)=rx(j)
!  av(i)=rav(j)
!  jnu(:,i)=rjnu(:,j)
!  abun(:,i)=rabun(:,j)
!  Tgas(i)=rTgas(j)
!  rho(i)=rrho(j)
!  pop(:,:,i)=rpop(:,:,j)
!  tau(:,i)=rtau(:,j)
!end do

write(6,*) 'Density=',rho(1)
write(6,*) 'Temperature=',Tgas(1)

return
end subroutine

subroutine constants

c=2.9979246d10 !cm/s
kb=1.380650d-16 !erg / K  *or*  g cm^2 / K s^2
hp=6.6260696e-27 !erg s  *or*  g cm^2 / s
pc2cm=3.085677d+18
pi=3.1415927
mhp=1.6726218d-24

freq0(1)=1900.5369d9     ;g(1,1)=2.0  ;g(1,2)=4.0   ; spec(1)="CII 158um"
freq0(2)=492.16065d9     ;g(2,1)=1.0  ;g(2,2)=3.0   ; spec(2)="CI (1-0)"
freq0(3)=1301.50262d9    ;g(3,1)=1.0  ;g(3,2)=5.0   ; spec(3)="CI (2-0)"
freq0(4)=809.34197d9     ;g(4,1)=3.0  ;g(4,2)=5.0   ; spec(4)="CI (2-1)"
freq0(5)=4744.77749d9    ;g(5,1)=5.0  ;g(5,2)=3.0   ; spec(5)="OI  1-0 "
freq0(6)=6804.84658d9    ;g(6,1)=5.0  ;g(6,2)=1.0   ; spec(6)="OI  2-0 "
freq0(7)=2060.06909d9    ;g(7,1)=3.0  ;g(7,2)=1.0   ; spec(7)="OI  2-1 "
freq0(8)=115.2712018d9   ;g(8,1)=1.0  ;g(8,2)=3.0   ; spec(8)="CO (1-0)"
freq0(9)=230.538d9       ;g(9,1)=3.0  ;g(9,2)=5.0   ; spec(9)="CO (2-1)"
freq0(10)=345.7959899d9  ;g(10,1)=5.0 ;g(10,2)=7.0  ; spec(10)="CO (3-2)"
freq0(11)=461.040768d9   ;g(11,1)=7.0 ;g(11,2)=9.0  ; spec(11)="CO (4-3)"
freq0(12)=576.2679305d9  ;g(12,1)=9.0 ;g(12,2)=11.0 ; spec(12)="CO (5-4)"
freq0(13)=691.4730763d9  ;g(13,1)=11.0;g(13,2)=13.0 ; spec(13)="CO (6-5)"
freq0(14)=806.6518060d9  ;g(14,1)=13.0;g(14,2)=15.0 ; spec(14)="CO (7-6)"
freq0(15)=921.7997000d9  ;g(15,1)=15.0;g(15,2)=17.0 ; spec(15)="CO (8-7)"
freq0(16)=1036.9123930d9 ;g(16,1)=17.0;g(16,2)=19.0 ; spec(16)="CO (9-8)"
freq0(17)=1151.9854520d9 ;g(17,1)=19.0;g(17,2)=21.0 ; spec(17)="CO (10-9)"

!Einstein A coefficients
A(1)=2.321d-06 ;  mh(1)=12.*mhp  
A(2)=7.880d-08 ;  mh(2)=12.*mhp  
A(3)=1.810d-14 ;  mh(3)=12.*mhp  
A(4)=2.650d-07 ;  mh(4)=12.*mhp  
A(5)=8.910d-05 ;  mh(5)=16.*mhp  
A(6)=1.340d-10 ;  mh(6)=16.*mhp  
A(7)=1.750d-05 ;  mh(7)=16.*mhp  
A(8)=7.203d-08 ;  mh(8)=28.*mhp  
A(9)=6.910d-07 ;  mh(9)=28.*mhp  
A(10)=2.497d-06;  mh(10)=28.*mhp 
A(11)=6.126d-06;  mh(11)=28.*mhp 
A(12)=1.221d-05;  mh(12)=28.*mhp 
A(13)=2.137d-05;  mh(13)=28.*mhp 
A(14)=3.422d-05;  mh(14)=28.*mhp 
A(15)=5.134d-05;  mh(15)=28.*mhp 
A(16)=7.330d-05;  mh(16)=28.*mhp 
A(17)=1.006d-04;  mh(17)=28.*mhp 

return
end subroutine

end program

!     Make a pdb with an array of dppc molecules
      Program GetLipBox
      IMPLICIT NONE

      integer,parameter :: nat = 130 ! number of atoms in the lipid
      integer,parameter :: nx = 4 ! number of lipids in x dimension
      integer,parameter :: ny = 4 
      integer,parameter :: nz = 1
      real,parameter :: dx = 10.0 ! distance between lipids in x-dim
      real,parameter :: dy = 14.0
      real,parameter :: dz = 5.0
      real,parameter :: PI = 3.14159265

      ! VARIABLES
      CHARACTER*23 s1(nat),s2(nat)
      CHARACTER*78 lpck(1000),lp2(1000)
      REAL*8 x(nat,3),xn(nat,3),p(1000,3),pn(1000,3)
      INTEGER r(nat)
      INTEGER i,j,k,l,m,cnt
      ! - - - - - - 

      open(file='dppc1.pdb',unit=1)
      open(file='dppc-mono.pdb',unit=2)

      read(1,*) ! first line in pdb is text only

      ! read in the coordinates of 1 lipid
      do i=1,nat
      read(1,'(a23,i3,f12.3,2f8.3,a24)') s1(i),r(i),
     &    (x(i,k),k=1,3), s2(i)
      enddo

      cnt=0

      ! make an aray (monolayer) of lipids
      do 40 m=1,nx ! x
      do 30 l=1,ny ! y
      do 20 j=1,nz ! z
      cnt=cnt+1
      do 10 i=1,nat
      xn(i,1) = x(i,1) + (m-1)*dx 
      xn(i,2) = x(i,2) + (l-1)*dy
      xn(i,3) = x(i,3) + (j-1)*dz
      write(2,'(a23,i3,f12.3,2f8.3,a24)')s1(i),cnt,(xn(i,k),k=1,3),
     &     s2(i)
10    continue
20    continue
30    continue
40    continue


      close(unit=1)
      close(unit=2)

      END
    

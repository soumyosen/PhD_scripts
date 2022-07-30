!     Make a pdb with an array of dppc molecules
      Program GetLipBox
      IMPLICIT NONE

      integer,parameter :: nat = 130 ! number of atoms in the lipid
      real,parameter :: PI = 3.14159265

      ! VARIABLES
      CHARACTER*23 s1(nat),s2(nat)
      CHARACTER*78 lpck(1000),lp2(1000)
      REAL*8 x(nat,3),xn(nat,3)
      INTEGER r(nat)
      INTEGER i,j,k,l,m,cnt
      ! - - - - - - 

      open(file='dppc1.pdb',unit=1)
      open(file='dppc-neg.pdb',unit=2)

      read(1,*) ! first line in pdb is text only

      ! read in the coordinates of 1 lipid, write out neg coords
      do i=1,nat
      read(1,'(a23,i3,f12.3,2f8.3,a24)') s1(i),r(i),
     &    (x(i,k),k=1,3), s2(i)

      xn(i,1)=x(i,1)
      xn(i,2)=x(i,2)
      xn(i,3)=-1.*x(i,3)
      write(2,'(a23,i3,f12.3,2f8.3,a24)')s1(i),1,(xn(i,k),k=1,3),
     &     s2(i)
      enddo

      close(unit=1)
      close(unit=2)

      END
    

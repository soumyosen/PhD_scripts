      PROGRAM propagate
      IMPLICIT NONE

      character*30 blah, ch(1000000,1)
      character*22 bck
      character*30 frt
      integer i, j, ml, kk,s,t,nx,ny,nz,r
      real*8 crd(1000000,3),c2(1000000,3),a,b,d

      open(file='PCL-PEG-NH2-1V.pdb',unit=1)
      open(file='PCL-PEG-NH2-30.pdb',unit=2)

c **********declare variables**************************************
      ml=1632
      a=80.
      b=5.
      d=5.   ! d=4*RADIUS+2*1.41 for nano tube
      nx=29
      ny=0
      nz=0
c read the header from the gjf
c      do 10 i=1,8
c       read(1,*)
c10    continue

c read and write the coordinates
       do 20 i=1,ml
         read(1,100) ch(i,1),(crd(i,kk),kk=1,3)
        write(2,100) ch(i,1),(crd(i,kk),kk=1,3)
20     continue
       
       do 40 s=1,nx
           do 30 i=1,ml
           
        crd(i+s*ml,1)=crd(i,1)+a*s
        crd(i+s*ml,2)=crd(i,2)
        crd(i+s*ml,3)=crd(i,3)
      
        ch(i+s*ml,1)=ch(i,1)
        write(2,100) ch(i,1),(crd(i+s*ml,kk),kk=1,3)
          
30     continue
40     continue     

       ml=(1+nx)*ml 
       
       do 50 t=1,ny
       do 60 i=1,ml

        crd(i+t*ml,1)=crd(i,1)
        crd(i+t*ml,2)=crd(i,2)+t*b
        crd(i+t*ml,3)=crd(i,3)

        ch(i+t*ml,1)=ch(i,1)
        write(2,100) ch(i,1),(crd(i+t*ml,kk),kk=1,3)
      
60     continue
50     continue

       ml=(1+ny)*ml 
       do 70 r=1,nz
       do 80 i=1,ml

        crd(i+r*ml,1)=crd(i,1)
        crd(i+r*ml,2)=crd(i,2)
        crd(i+r*ml,3)=crd(i,3)+r*d

c        ch(i+r*ml,1)=ch(i,1)
        write(2,100) ch(i,1),(crd(i+r*ml,kk),kk=1,3)

80     continue
70     continue
 
c******* formats and close out of open files*****************

       close(unit=1)
       close(unit=2)

100    format(a30 f8.3,f8.3,f8.3,a22,2f6.2)

       end


      module m1
      integer,parameter::BIG=80000
      real::dedge,edge,dBond,rPen,PenSpac,PenDist,PenDist72,PenDistE
      real::midSpac,midE,E,s2
      integer::magic,num,a=0,b=0
      real::o(BIG,3)=0,o2(BIG,3)=0
      integer::b1=0,b2=0,b3=0,b4=0,b5=0
      integer::ni,lay
      real::PI=3.1415926535,mass
      real:: dist,dist2,dist3
      integer::p(BIG,4)=0!1=shell,2=layer,3=position,4=leaf
      !========================
      real::Bond,dB,o3(BIG,3)=0,ang,rv,d
      real::c1,c2,c3,c4,c5
      integer::shell(BIG),shellr(999,2),atom=0,shellc=0,keep
      character*4 typ(BIG)
      end module m1
      
      program nanoparticle
      use m1
      implicit none
      
      open(file='core.pdb',unit=1)
      open(file='bond.0',unit=2)
      open(file='atom.0',unit=3)
      open(file='ang.0',unit=4)
      open(file='core.coor',unit=5)

c      Bond=2.74;dB=0.3 !Au-Au bond distance (A), acceptable range 
c      rv=Bond/(2*sin(31.7*PI/180))
c
c      !----------------------------------------
c      !center (0.95 smaller atomic radius)
c      !----------------------------------------
c      atom=atom+1;shell(atom)=0;typ(atom)="Au0 "
c      !----------------------------------------
c      !1st shell (icosahedron)
c      !----------------------------------------
c      shellc=shellc+1
c      ! mark shell range
c      shellr(1,1)=2;shellr(1,2)=13
c      do b4=1,12
c      atom=atom+1;shell(atom)=shellc;typ(atom)="Au3 "
c      if (b4.eq.1) ang=0
c      if (b4.eq.7) ang=36
c      c1=cos(ang*PI/180)*cos(26.56505*PI/180)
c      c2=sin(ang*PI/180)*cos(26.56505*PI/180)
c      c3=sin(26.56505*PI/180)
c      c4=1;if (b4.gt.6) c4=-1
c      if ((b4.eq.6).OR.(b4.eq.12)) then
c      o3(atom,1)=0;o3(atom,2)=0;o3(atom,3)=c4*rv
c      else
c      o3(atom,1)=rv*c1;o3(atom,2)=rv*c2;o3(atom,3)=c4*rv*c3
c      endif
c      ang=ang+72
c      enddo
c      !----------------------------------------
c      !2nd shell 
c      !tetrahedral packing = add 1 to each tri
c      !octahedral packing=add tri each tri
c      !----------------------------------------
c      shellc=shellc+1
c      do b1=shellr(shellc-1,1),shellr(shellc-1,2)
c      if (shell(b1).ne.shellc-1) cycle
c      do b2=b1+1,shellr(shellc-1,2)
c      if (shell(b2).ne.shellc-1) cycle
c      call Finddist(b1,b2);if ((d.lt.Bond-dB).OR.(d.gt.Bond+dB)) cycle
c      do b3=b2+1,shellr(shellc-1,2)
c      if (shell(b3).ne.shellc-1) cycle
c      call Finddist(b1,b3);if ((d.lt.Bond-dB).OR.(d.gt.Bond+dB)) cycle
c      call Finddist(b2,b3);if ((d.lt.Bond-dB).OR.(d.gt.Bond+dB)) cycle
c      ! for each triangle do
c      atom=atom+1;shell(atom)=shellc;typ(atom)="Au4 "
c      ! mark shell range
c      if (shell(atom-1).eq.shellc-1) shellr(shellc,1)=atom
c      shellr(shellc,2)=atom
c      ! find coor of new atom
c      c1=(o3(b1,1)+o3(b2,1)+o3(b3,1))/3 !tri center x
c      c2=(o3(b1,2)+o3(b2,2)+o3(b3,2))/3 !tri center y
c      c3=(o3(b1,3)+o3(b2,3)+o3(b3,3))/3 !tri center z
c      c4=sqrt(c1*c1+c2*c2+c3*c3) !length of (c1,c2,c3)
c      c5=sqrt(Bond**2*2/3) !height of pyramid
c      o3(atom,1)=c1*(1+c5/c4)
c      o3(atom,2)=c2*(1+c5/c4)
c      o3(atom,3)=c3*(1+c5/c4)
c      enddo;enddo;enddo
c
c      do b1=1,atom
c      write(5,'(i5,a10,a4,3f8.3)') shell(b1)," shell Au ",typ(b1),
c     &o3(b1,1),o3(b1,2),o3(b1,3)
c      enddo



      !==========================
      ! Icosahedral core
      !==========================
      dedge=22 ! edge to edge diameter
      edge=dedge/(2*0.80901699) !Length of edge (A)
      dBond=2.74 !Au-Au bond distance (A)
      num=nint(edge/dBond)+1 !number of atoms on one edge + 1
      lay=num-1
      magic=(10*lay**3+15*lay**2+11*lay+3)/3 !magic number of atom

      !central atom at (0,0,0)
c     atom=atom+1
      if (atom.ne.0) typ(atom)="Au1 "

      !b1 define building from which shell to which shell
      !do b1=num-1-1,num-1-1 !num-1 is the outer-most shell
      do b1=num-1-1,num-1-1 !num-1 is the outer-most shell
      rPen=dBond*b1/2 !height from origin to base of Pentagon
      PenSpac=rPen/b1 !vertical spacing between Pentagon
      !top and bottom vertice
      do b4=0,1
      atom=atom+1
      typ(atom)="Au3 "
      p(atom,1)=b1;p(atom,2)=1;p(atom,3)=b3;p(atom,4)=b4
      o(atom,1)=dBond*b1; o(atom,3)=180*b4
      enddo
      !pentagon region
      do b2=2,b1+1 !b2th layer from the top vertice
      PenDist72=PenSpac*(b2-1)*tan(60*PI/180) !center to edge
      PenDistE=2*PenDist72*sin(36*PI/180)/(b2-1) !edge spacing
      do b4=1,5 !5-fold symetry to build pentagon
      do b3=0,b2-2 !b3 is with atom from the edge of one leaf
      PenDist=sqrt((PenDistE*b3)**2+PenDist72**2-
     &2*(PenDistE*b3)*PenDist72*cos(54*PI/180))
      ang=acos(-1*((PenDistE*b3)**2-PenDist72**2-PenDist**2)/
     &(2*PenDist72*PenDist)) !in rad
      atom=atom+1 !top pentagon
      call FindType
      p(atom,1)=b1;p(atom,2)=b2;p(atom,3)=b3;p(atom,4)=b4
      o(atom,1)=sqrt((2*rPen-PenSpac*(b2-1))**2+PenDist**2)
      o(atom,2)=b4*72+ang*180/PI
      o(atom,3)=atan(PenDist/(2*rPen-PenSpac*(b2-1)))*180/PI
      atom=atom+1 !bottom pentagon
      call FindType
      p(atom,1)=b1;p(atom,2)=b2;p(atom,3)=b3;p(atom,4)=b4
      o(atom,1)=o(atom-1,1)
      o(atom,2)=o(atom-1,2)+36+72
      o(atom,3)=180-o(atom-1,3)
      !in between the two pentagons
      if (b2.ge.(b1+1)) cycle
      midSpac=2*rPen/b1
      atom=atom+1
      call FindType
      p(atom,1)=b1;p(atom,2)=b2;p(atom,3)=b3;p(atom,4)=b4
      o(atom,1)=o(atom-2,1)
      o(atom,3)=acos((rPen-midSpac*(b2-1))/o(atom,1))*180/PI
      s2=sin(o(atom,3)*PI/180)*o(atom,1)
      midE=PenDistE*(b2-1)/2
      E=abs(b3*PenDistE-midE)
      o(atom,2)=72*b4-asin(E/s2)*180/PI
      if (midE.gt.b3*PenDistE) o(atom,2)=72*b4+asin(E/s2)*180/PI
      atom=atom+1
      call FindType
      p(atom,1)=b1;p(atom,2)=b2;p(atom,3)=b3;p(atom,4)=b4
      o(atom,1)=o(atom-2,1)
      o(atom,2)=o(atom-1,2)+36+72
      o(atom,3)=180-o(atom-1,3)
      enddo; enddo; enddo; enddo
      !spherical(degree) to cartesian
      do b1=1,atom
      o2(b1,1)=o(b1,1)*sin(o(b1,3)*PI/180)*cos(o(b1,2)*PI/180)
      o2(b1,2)=o(b1,1)*sin(o(b1,3)*PI/180)*sin(o(b1,2)*PI/180)
      o2(b1,3)=o(b1,1)*cos(o(b1,3)*PI/180)
      enddo
      ! generate pdb
      do b1=1,atom
      ni=ni+1
      write(1,'(a4,i7,a19,3f8.3,a22)') 'ATOM',ni,
     &'  Au  NP1 G   1    ',o2(b1,1),o2(b1,2),
     &o2(b1,3),'  0.00  0.00      NANO'
      enddo
      ! Write atom.0
      mass=196.9665*magic/atom
      do b1=1,atom
      write(3,'(i8,a21,a4,a1,a11,f12.4,a12)') b1," NANO 1    NP   Au   "
     &,typ(b1)," ","  0.000000 ",mass,"           0"
      do b2=1,atom
      dist=sqrt((o2(b1,1)-o2(b2,1))**2+(o2(b1,2)-o2(b2,2))**2+
     &(o2(b1,3)-o2(b2,3))**2)
      ! Write bond.0
c     if (b2.gt.b1) then
      if ((dist.gt.dBond+1.0).OR.(dist.lt.dBond-1.0)) cycle
      if ((b2.gt.b1).AND.(typ(b1).ne."Au2").AND.
     &(typ(b2).ne."Au2")) then
      write(2,'(2i8)') b1,b2
      endif
      !Find angle
      do b3=b1+1,atom
      if (b3.eq.b2) cycle
      !same shell, linear?
      if ((p(b1,1).ne.p(b2,1)).OR.(p(b2,1).ne.p(b3,1))) cycle
      dist2=sqrt((o2(b3,1)-o2(b2,1))**2+(o2(b3,2)-o2(b2,2))**2+
     &(o2(b3,3)-o2(b2,3))**2)
      if ((dist2.gt.dBond+1.0).OR.(dist2.lt.dBond-1.0)) cycle
      dist3=sqrt((o2(b1,1)-o2(b3,1))**2+(o2(b1,2)-o2(b3,2))**2+
     &(o2(b1,3)-o2(b3,3))**2)
      ang=acos((dist3**2-dist**2-dist2**2)/(-2*dist*dist2))*180/PI
      if (ang.le.160) cycle
      !side edge?
      if ((p(b2,2).ge.1).AND.(p(b2,2).lt.p(b2,1)+1)) then
      if ((p(b2,3).eq.0).AND.(p(b1,3).eq.0).AND.(p(b3,3).eq.0)) then
      write(4,'(3i8)') b1,b2,b3
      endif; endif
      !bottom edge?
      if (p(b2,2).eq.p(b2,1)+1) then
      if (p(b2,3).gt.0) then
      if ((p(b1,2).eq.p(b2,2)).AND.(p(b3,2).eq.p(b2,2))) then
      write(4,'(3i8)') b1,b2,b3
      endif; endif; endif
      !face
      if ((p(b2,3).ge.1).AND.(p(b2,2).ne.p(b2,1)+1)) then
      write(4,'(3i8)') b1,b2,b3
      endif
      enddo; enddo; enddo


      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)      
      close(unit=5)      
      END


      !==========================
      ! Subroutines
      !==========================
      subroutine FindType()
      use m1
      implicit none
      typ(atom)="Au2 "
      if ((b3.eq.0).OR.(b2.eq.b1+1)) typ(atom)="Au1 "
      if (b2.eq.2) typ(atom)="Au3 "
      if ((b3.le.0).AND.(b2.eq.b1)) typ(atom)="Au3 "
      if ((b3.le.1).AND.(b2.eq.b1+1)) typ(atom)="Au3 "
      if ((b3.eq.b2-2).AND.(b2.eq.b1+1)) typ(atom)="Au3 "
      end subroutine

      real function distanceformula(x1,y1,z1,x2,y2,z2)
      real x1,y1,z1,x2,y2,z2
      distanceformula=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      return
      end function

      subroutine Finddist(atom1,atom2)
      use m1
      implicit none
      integer atom1,atom2
      real distanceformula
      d=0
      d=distanceformula(o3(atom1,1),o3(atom1,2),o3(atom1,3),
     &o3(atom2,1),o3(atom2,2),o3(atom2,3))
      end subroutine


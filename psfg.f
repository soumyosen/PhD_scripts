      MODULE MM       !-------- global variables and params
      !-----------constants-----------
      !--------- immutable constansts ----------
        INTEGER,PARAMETER :: BIG = 800000
        INTEGER,PARAMETER :: TYPSZ = 300
        REAL(8),PARAMETER :: PI = 3.141592
        REAL(8),PARAMETER :: d2r = PI/180.
        REAL(8),PARAMETER :: kc = 14.52488 !14.3988  [eV*A/e^2]
        REAL(8),PARAMETER :: kcdeb = 4.8032 ! [Debye/{A e}]
        REAL(8),PARAMETER :: kf = 1./kcdeb ! [{A e}/Debye] if E = V/A
        REAL(8),PARAMETER :: kdip = 0.6295 ! [eV A^3/D^2]
        REAL(8),PARAMETER :: kbev=0.00008617
        REAL(8),PARAMETER :: s2=1.414214     !sqrt(2.)
        REAL(8),PARAMETER :: s3=1.732050808  !sqrt(3.)
        REAL(8),PARAMETER :: cb=1.41  ! c-c bond
      !------- common for subroutines
        REAL(8) rc(3),rs(3)
        REAL(8) ts(BIG,3)  ! --- temporary (sub)structure
        INTEGER tst(BIG,2),tsz ! --- atom type/visibility
        INTEGER tbo(BIG,3),ntbo ! --- bonds
        INTEGER ta(BIG,4),nta  ! --- angles
        INTEGER td(BIG,5),ntd  ! --- dihedrals
        INTEGER cn ! --- final water atoms 
        INTEGER Nm ! --- number of atoms in molecule
      !-----------global structure-----------
        REAL(8) s(BIG,3)
        INTEGER st(BIG,5),ssz ! --- atom type/visibility/resname/resid/segid
        INTEGER bo(BIG,3),bsz
        INTEGER an(BIG,4),asz
        INTEGER di(BIG,5),dsz
        !--- excluded bonds
        CHARACTER*5 exb(TYPSZ,2) !--- types
        INTEGER nexb
        REAL(8) exbb(TYPSZ) !---- tolerance distance

        INTEGER ntyp
        REAL(8) masstb(TYPSZ),chartb(TYPSZ)
        CHARACTER*5 segtb(TYPSZ),resntb(TYPSZ),typtb(TYPSZ)
        CHARACTER*5 aname(BIG),rname(BIG)
      END MODULE MM

      !==================
      PROGRAM MAIN
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k,it,kt,ct,ct2
      ! - - - begin - - -

      CALL readfiles

      tsz = ssz!1988!2040!960
      do i=1,tsz
        do k=1,3
          ts(i,k) = s(i,k)
        enddo
c        write(*,*) (ts(i,k),k=1,3)
      enddo

      !CALL testpsf
      CALL psfgen
      
      ssz = tsz
      bsz = ntbo
      asz = nta
      dsz = ntd

      do i=1,ssz
        st(i,2) = 1
      enddo
      do i=1,bsz
        do k=1,3
          bo(i,k)=tbo(i,k)
        enddo
      enddo
      do i=1,asz
        do k=1,4
          an(i,k)=ta(i,k)
        enddo
      enddo
      do i=1,dsz
        do k=1,5
          di(i,k)=td(i,k)
        enddo
      enddo


c      CALL waterbox

      CALL savepsf
      CALL savepdb
c      CALL modpdb

      END PROGRAM
      !=================

      ! --------------------------------------------------
      SUBROUTINE readfiles() !read psf and pdb files
      ! --------------------------------------------------
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k, it,kt,ct,ct2
      INTEGER na
      REAL(8) th,si,co
      CHARACTER*4 buf4
      CHARACTER*5 buf5
      CHARACTER*6 buf6
      CHARACTER*64 buf64
     
      open(file='atompsf_sopc.cfg',unit=1)

      read(1,'(i5)') ntyp ! number of recognized atom types

      write(*,*) " [readfiles] recognized atom types:",ntyp
      do i=1,ntyp
        read(1,100)
     &        j,segtb(i),k,resntb(i),buf5,typtb(i),chartb(i),masstb(i)
c        write(*,100)
c     &        j,segtb(i),k,resntb(i),buf5,typtb(i),chartb(i),masstb(i)
      enddo

                    
      read(1,'(i5)') nexb

c      write(*,'(i5)') nexb
      write(*,*) " [readfiles] bond exclusion between atom types:",nexb
      do i=1,nexb
        read(1,110) exb(i,1),exb(i,2),exbb(i),th 
      enddo
      close(unit=1)

      open(file='chol1.pdb',unit=1)

c      open(file='polar5.pdb',unit=1)
      read(1,'(a6,i8)') buf6,ssz  ! REMARK number of atoms
      Nm = ssz 

      do i=1,ssz
c        read(1,101) buf4,j,aname(i),rname(i),
c     &   k,(s(i,k),k=1,3),th,th,buf5,buf4
        read(1,102) buf4,j,aname(i),rname(i),
     &   k,(s(i,k),k=1,3)

        if(rname(i) == ' BNT ')then
          st(i,1) = 1
        endif

        if(rname(i) == ' NBT ')then
          if (aname(i) == '   B ') then
          st(i,1) = 160
          endif
          if (aname(i) == '   N ') then
          st(i,1) = 161
          endif
        endif

        if(rname(i) ==   ' GLY ')then
          if(aname(i) == ' C   ')then
            st(i,1) = 2
          elseif(aname(i) == ' O   ')then
            st(i,1) = 3
          elseif(aname(i) == ' N   ')then
            st(i,1) = 47
          elseif((aname(i) == ' HA1 ').or.
     &           (aname(i) == ' HA2 '))then
            st(i,1) = 48
          elseif(aname(i) == ' HN  ')then
            st(i,1) = 49
          elseif(aname(i) == ' CA  ')then
            st(i,1) = 50
          endif
        endif
        if(rname(i) ==   ' ASP ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 4
          elseif((aname(i) == ' OD2 ').or.
     &          (aname(i) == ' OD1 '))then
            st(i,1) = 5
          elseif(aname(i) == ' CG  ')then
            st(i,1) = 6
          elseif((aname(i) == ' HB1 ').or.
     &          (aname(i) == ' HB2 '))then
            st(i,1) = 7
          elseif(aname(i) == ' HA  ')then
            st(i,1) = 8
          elseif(aname(i) == ' CB  ')then
            st(i,1) = 9
          elseif(aname(i) == ' HN  ')then
            st(i,1) = 10
          elseif(aname(i) == ' O   ')then
            st(i,1) = 11
          elseif(aname(i) == ' CA  ')then
            st(i,1) = 12
          elseif(aname(i) == ' C   ')then
            st(i,1) = 13
          endif
        elseif(rname(i) ==   ' SER ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 14
          elseif((aname(i) == ' HB2 ').or.
     &          (aname(i) == ' HB1 '))then
            st(i,1) = 15
          elseif(aname(i) == ' OG  ')then
            st(i,1) = 16
          elseif(aname(i) == ' HG  ')then
            st(i,1) = 17
          elseif((aname(i) == ' HB1 ').or.
     &          (aname(i) == ' HB2 ').or.
     &          (aname(i) == ' HB3 '))then
            st(i,1) = 18
          elseif(aname(i) == ' HA  ')then
            st(i,1) = 19
          elseif(aname(i) == ' CB  ')then
            st(i,1) = 20
          elseif(aname(i) == ' HN  ')then
            st(i,1) = 21
          elseif(aname(i) == ' O   ')then
            st(i,1) = 22
          elseif(aname(i) == ' CA  ')then
            st(i,1) = 23
          elseif(aname(i) == ' C   ')then
            st(i,1) = 24
          endif
        elseif(rname(i) ==   ' LYS ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 25
          elseif((aname(i) == ' HZ1 ').or.
     &          (aname(i) == ' HZ2 ').or.
     &          (aname(i) == ' HZ3 '))then
            st(i,1) = 26
          elseif((aname(i) == ' HE1 ').or.
     &          (aname(i) == ' HE2 ').or.
     &          (aname(i) == ' HE3 '))then
            st(i,1) = 27
          elseif(aname(i) == ' NZ  ')then
            st(i,1) = 28
          elseif(aname(i) == ' HA  ')then
            st(i,1) = 29
          elseif(aname(i) == '  CE ')then
            st(i,1) = 30
          elseif(aname(i) == ' HN  ')then
            st(i,1) = 31
          elseif(aname(i) == ' O   ')then
            st(i,1) = 32
          elseif(aname(i) == ' CA  ')then
            st(i,1) = 33
          elseif(aname(i) == ' C   ')then
            st(i,1) = 34
          endif
        elseif(rname(i) ==   ' VAL ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 35
          elseif((aname(i) == ' HG1 ').or.
     &          (aname(i) == ' HG2 ').or.
     &          (aname(i) == ' HB  ').or.
     &          (aname(i) == ' HG3 '))then
            st(i,1) = 36
          elseif((aname(i) == ' CG1 ').or.
     &          (aname(i) == ' CG2 ').or.
     &          (aname(i) == ' CG3 '))then
            st(i,1) = 37
          elseif(aname(i) == ' HA  ')then
            st(i,1) = 38
          elseif(aname(i) == ' CB  ')then
            st(i,1) = 39
          elseif(aname(i) == ' HN  ')then
            st(i,1) = 40
          elseif(aname(i) == ' O   ')then
            st(i,1) = 41
          elseif(aname(i) == ' CA  ')then
            st(i,1) = 42
          elseif(aname(i) == ' C   ')then
            st(i,1) = 43
          endif
        elseif(rname(i) ==   ' IP3 ')then
          if(aname(i) == '   O ')then
            st(i,1) = 44
          elseif((aname(i) == '  H1 '))then
            st(i,1) = 45
          elseif(aname(i) ==  '  H2 ')then
            st(i,1) = 46
          endif
        elseif(rname(i) ==   ' PHE ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 51
          elseif((aname(i) == '  HE '))then
            st(i,1) = 52
          elseif((aname(i) == '  HD '))then
            st(i,1) = 53
          elseif((aname(i) == ' HB  '))then
            st(i,1) = 54
          elseif((aname(i) == ' CZ  '))then
            st(i,1) = 55
          elseif((aname(i) == ' CE1 '))then
            st(i,1) = 56
          elseif((aname(i) == ' CE2 '))then
            st(i,1) = 57
          elseif((aname(i) == ' CD1 '))then
            st(i,1) = 58
          elseif((aname(i) == ' CD2 '))then
            st(i,1) = 59
          elseif((aname(i) == '  CG '))then
            st(i,1) = 60
          elseif((aname(i) == ' HZ  '))then
            st(i,1) = 61
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 62
          elseif((aname(i) == ' CB  '))then
            st(i,1) = 63
          elseif((aname(i) == ' HN  '))then
            st(i,1) = 64
          elseif(aname(i) == ' O   ')then
            st(i,1) = 65
          elseif(aname(i) == ' CA  ')then
            st(i,1) = 66
          elseif(aname(i) == ' C   ')then
            st(i,1) = 67
          elseif((aname(i) == '  GR '))then
            st(i,1) = 162
          endif
        elseif(rname(i) ==   ' PHO ')then
          if(aname(i) == ' PH  ')then
            st(i,1) = 157
          endif
        elseif(rname(i) ==   ' CYS ')then
          if(aname(i) == ' SG  ')then
            st(i,1) = 158
          endif
        elseif(rname(i) ==   ' ARG ')then
          if(aname(i) == ' NE  ')then
            st(i,1) = 159
          endif
        elseif(rname(i) ==   ' MTL ')then
          if(aname(i) == ' CT3 ')then
            st(i,1) = 68
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 69
          elseif((aname(i) == ' OH1 '))then
            st(i,1) = 70
          elseif((aname(i) == ' H   '))then
            st(i,1) = 71
          endif
        elseif(rname(i) ==   ' ETL ')then
          if(aname(i) == ' C2  ')then
            st(i,1) = 72
          elseif((aname(i) == ' HB  '))then
            st(i,1) = 73
          elseif((aname(i) == ' OH1 '))then
            st(i,1) = 74
          elseif((aname(i) == ' H   '))then
            st(i,1) = 75
          elseif(aname(i) == ' C3  ')then
            st(i,1) = 76
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 77
          endif
        elseif(rname(i) ==   ' SUG ')then
          if(aname(i) == ' C5  ')then
            st(i,1) = 78
          elseif((aname(i) == ' H5  '))then
            st(i,1) = 79
          elseif((aname(i) == ' O5  '))then
            st(i,1) = 80
          elseif((aname(i) == ' C1  '))then
            st(i,1) = 81
          elseif((aname(i) == ' H1  '))then
            st(i,1) = 82
          elseif((aname(i) == ' O1  '))then
            st(i,1) = 83
          elseif((aname(i) == ' HO1 '))then
            st(i,1) = 84
          elseif((aname(i) == ' C2  '))then
            st(i,1) = 85
          elseif((aname(i) == ' H2  '))then
            st(i,1) = 86
          elseif((aname(i) == ' O2  '))then
            st(i,1) = 87
          elseif((aname(i) == ' HO2 '))then
            st(i,1) = 88
          elseif((aname(i) == ' C3  '))then
            st(i,1) = 89
          elseif((aname(i) == ' H3  '))then
            st(i,1) = 90
          elseif((aname(i) == ' O3  '))then
            st(i,1) = 91
          elseif((aname(i) == ' HO3 '))then
            st(i,1) = 92
          elseif((aname(i) == ' C4  '))then
            st(i,1) = 93
          elseif((aname(i) == ' H4  '))then
            st(i,1) = 94
          elseif((aname(i) == ' O4  '))then
            st(i,1) = 95
          elseif((aname(i) == ' HO4 '))then
            st(i,1) = 96
          elseif((aname(i) == ' C6  '))then
            st(i,1) = 97
          elseif((aname(i) == ' H61 '))then
            st(i,1) = 98
          elseif((aname(i) == ' O6  '))then
            st(i,1) = 99
          elseif((aname(i) == ' HO6 '))then
            st(i,1) = 100
          elseif((aname(i) == ' H62 '))then
            st(i,1) = 101
          endif
        elseif(rname(i) ==   ' TRP ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 102
          elseif((aname(i) == ' HZ2 '))then
            st(i,1) = 103
          elseif((aname(i) == ' HE3 '))then
            st(i,1) = 104
          elseif((aname(i) == ' HE1 '))then
            st(i,1) = 105
          elseif((aname(i) == ' HD1 '))then
            st(i,1) = 106
          elseif((aname(i) == ' HB3 '))then
            st(i,1) = 107
          elseif((aname(i) == ' HB2 '))then
            st(i,1) = 108
          elseif((aname(i) == ' CH2 '))then
            st(i,1) = 109
          elseif((aname(i) == ' CZ3 '))then
            st(i,1) = 110
          elseif((aname(i) == ' CZ2 '))then
            st(i,1) = 111
          elseif((aname(i) == ' NE1 '))then
            st(i,1) = 112
          elseif((aname(i) == ' CE3 '))then
            st(i,1) = 113
          elseif((aname(i) == ' CE2 '))then
            st(i,1) = 114
          elseif((aname(i) == ' CD2 '))then
            st(i,1) = 115
          elseif((aname(i) == ' CD1 '))then
            st(i,1) = 116
          elseif((aname(i) == ' CG  '))then
            st(i,1) = 117
          elseif((aname(i) == ' HH2 '))then
            st(i,1) = 118
          elseif((aname(i) == ' HZ3 '))then
            st(i,1) = 119
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 120
          elseif((aname(i) == ' CB  '))then
            st(i,1) = 121
          elseif((aname(i) == ' HN  '))then
            st(i,1) = 122
          elseif((aname(i) == ' O   '))then
            st(i,1) = 123
          elseif((aname(i) == ' CA  '))then
            st(i,1) = 124
          elseif((aname(i) == ' C   '))then
            st(i,1) = 125
          endif
        elseif(rname(i) ==   ' IPRL')then
          if(aname(i) == ' CO  ')then
            st(i,1) = 126
          elseif((aname(i) == ' OH1 '))then
            st(i,1) = 127
          elseif((aname(i) == ' H   '))then
            st(i,1) = 128
          elseif((aname(i) == ' HS  '))then
            st(i,1) = 129
          elseif((aname(i) == ' C3  '))then
            st(i,1) = 130
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 131
          endif
        elseif(rname(i) ==   ' PRL ')then
          if(aname(i) == ' CO  ')then
            st(i,1) = 132
          elseif((aname(i) == ' OH1 '))then
            st(i,1) = 133
          elseif((aname(i) == ' H   '))then
            st(i,1) = 134
          elseif((aname(i) == ' C1  '))then
            st(i,1) = 135
          elseif((aname(i) == ' HB  '))then
            st(i,1) = 136
          elseif((aname(i) == ' C3  '))then
            st(i,1) = 137
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 138
          endif
       
        elseif(rname(i) ==   ' ETH ')then
          if(aname(i) == '  C1 ')then
            st(i,1) = 163
          elseif((aname(i) == ' H11 '))then
            st(i,1) = 164
          elseif((aname(i) == '  C2 '))then
            st(i,1) = 165
          elseif((aname(i) == ' H22 '))then
            st(i,1) = 166
          elseif((aname(i) == '  OG '))then
            st(i,1) = 167
          elseif((aname(i) == ' HG1 '))then
            st(i,1) = 168
          endif

         elseif(rname(i) ==   ' TYR ')then
          if(aname(i) == ' N   ')then
            st(i,1) = 139
          elseif((aname(i) == ' HE  '))then
            st(i,1) = 140
          elseif((aname(i) == ' HD  '))then
            st(i,1) = 141
          elseif((aname(i) == ' HB  '))then
            st(i,1) = 142
          elseif((aname(i) == ' CZ  '))then
            st(i,1) = 143 
          elseif((aname(i) == ' CE2 '))then
            st(i,1) = 144
          elseif((aname(i) == ' CE1 '))then
            st(i,1) = 145
          elseif((aname(i) == ' CD2 '))then
            st(i,1) = 146
          elseif((aname(i) == ' CD1 '))then
            st(i,1) = 147
          elseif((aname(i) == ' CG  '))then
            st(i,1) = 148
!            write(*,*) typtb(st(i,1)),chartb(st(i,1)),masstb(st(i,1))
          elseif((aname(i) == ' HH  '))then
            st(i,1) = 149
          elseif((aname(i) == ' OH1 '))then
            st(i,1) = 150
          elseif((aname(i) == ' HA  '))then
            st(i,1) = 151
          elseif((aname(i) == ' CB  '))then
            st(i,1) = 152
          elseif((aname(i) == ' HN  '))then
            st(i,1) = 153
          elseif((aname(i) == ' O   '))then
            st(i,1) = 154
          elseif((aname(i) == ' CA  '))then
            st(i,1) = 155
          elseif((aname(i) == ' C   '))then
            st(i,1) = 156
          endif

         elseif(rname(i) ==  ' SOPC')then
          if(aname(i) == '  N  ')then
            st(i,1) = 110
          elseif((aname(i) == '  C12'))then
            st(i,1) = 111
          elseif((aname(i) == '  C13').or.
     &           (aname(i) == '  C14').or.
     &           (aname(i) == '  C15'))then
            st(i,1) = 112
          elseif((aname(i) == ' H12A'))then
            st(i,1) = 115
          elseif((aname(i) == '  C11'))then
            st(i,1) = 123
          elseif((aname(i) == '  HA '))then
            st(i,1) = 124
          elseif((aname(i) == '  O12'))then
            st(i,1) = 125
          elseif((aname(i) == '  O13'))then
            st(i,1) = 126
          elseif((aname(i) == '  P  '))then
            st(i,1) = 127
          elseif((aname(i) == '  C1 '))then
            st(i,1) = 128
!            write(*,*) typtb(st(i,1)),chartb(st(i,1)),masstb(st(i,1))
          elseif((aname(i) == '  C2 '))then
            st(i,1) = 129
          elseif((aname(i) == '  HS '))then
            st(i,1) = 130
          elseif((aname(i) == '  C3 '))then
            st(i,1) = 131
          elseif((aname(i) == '  O21'))then
            st(i,1) = 132
          elseif((aname(i) == '  O22'))then
            st(i,1) = 133
          elseif((aname(i) == '  C21'))then
            st(i,1) = 134
          elseif((aname(i) == '  C32'))then
            st(i,1) = 135
          elseif((aname(i) == '  C29'))then
            st(i,1) = 136
          elseif((aname(i) == '  H91'))then
            st(i,1) = 137
          elseif((aname(i) == ' C218'))then
            st(i,1) = 138
          elseif((aname(i) == '  C33'))then
            st(i,1) = 139
          elseif((aname(i) == ' H18R'))then
            st(i,1) = 140
          endif

            elseif(rname(i) ==  ' DCE ')then
          if(aname(i) == '  C  ')then
            st(i,1) = 141
          elseif((aname(i) == '  Cl '))then
            st(i,1) = 142
          elseif((aname(i) == '  H  '))then
            st(i,1) = 143
          endif
    
          elseif(rname(i) ==  ' CHL1')then
           if(aname(i) == '  C3 ')then
            st(i,1) = 144
          elseif((aname(i) == '  O3 '))then
            st(i,1) = 145
          elseif((aname(i) == '  H3 '))then
            st(i,1) = 146
          elseif((aname(i) == '  H3 '))then
            st(i,1) = 147
          elseif((aname(i) == '  C4 '))then
            st(i,1) = 148
          elseif((aname(i) == '  H4A'))then
            st(i,1) = 149
          elseif((aname(i) == '  H4B'))then
            st(i,1) = 150
          elseif((aname(i) == '  C5 '))then
            st(i,1) = 151
          elseif((aname(i) == '  C6 '))then
            st(i,1) = 152
          elseif((aname(i) == '  H6 '))then
            st(i,1) = 153
          elseif((aname(i) == '  C7 '))then
            st(i,1) = 154
          elseif((aname(i) == '  H7A'))then
            st(i,1) = 155
          elseif((aname(i) == '  H7B'))then
            st(i,1) = 156
          elseif((aname(i) == '  C8 '))then
            st(i,1) = 157
          elseif((aname(i) == '  H8 '))then
            st(i,1) = 158
          elseif((aname(i) == '  C14'))then
            st(i,1) = 159
          elseif((aname(i) == '  H14'))then
            st(i,1) = 160
          elseif((aname(i) == '  C15'))then
            st(i,1) = 161
          elseif((aname(i) == ' H15A'))then
            st(i,1) = 162
          elseif((aname(i) == ' H15B'))then
            st(i,1) = 163
          elseif((aname(i) == '  C16'))then
            st(i,1) = 164
          elseif((aname(i) == ' H16A'))then
            st(i,1) = 165
          elseif((aname(i) == ' H16B'))then
            st(i,1) = 166
          elseif((aname(i) == '  C17'))then
            st(i,1) = 167
          elseif((aname(i) == '  H17'))then
            st(i,1) = 168
          elseif((aname(i) == '  C13'))then
            st(i,1) = 169
          elseif((aname(i) == '  C18'))then
            st(i,1) = 170
          elseif((aname(i) == ' H18A'))then
            st(i,1) = 171
          elseif((aname(i) == ' H18B'))then
            st(i,1) = 172
          elseif((aname(i) == ' H18C'))then
            st(i,1) = 173
          elseif((aname(i) == '  C12'))then
            st(i,1) = 174
          elseif((aname(i) == ' H12A'))then
            st(i,1) = 175
          elseif((aname(i) == ' H12B'))then
            st(i,1) = 176
          elseif((aname(i) == '  C11'))then
            st(i,1) = 177
          elseif((aname(i) == ' H11A'))then
            st(i,1) = 178
          elseif((aname(i) == ' H11B'))then
            st(i,1) = 179
          elseif((aname(i) == '  C9 '))then
            st(i,1) = 180
          elseif((aname(i) == '  H9 '))then
            st(i,1) = 181
          elseif((aname(i) == '  C10'))then
            st(i,1) = 182
          elseif((aname(i) == '  C19'))then
            st(i,1) = 183
          elseif((aname(i) == ' H19A'))then
            st(i,1) = 184
          elseif((aname(i) == ' H19B'))then
            st(i,1) = 185
          elseif((aname(i) == ' H19C'))then
            st(i,1) = 186
          elseif((aname(i) == '  C1 '))then
            st(i,1) = 187
          elseif((aname(i) == '  H1A'))then
            st(i,1) = 188
          elseif((aname(i) == '  H1B'))then
            st(i,1) = 189
          elseif((aname(i) == '  C2 '))then
            st(i,1) = 190
          elseif((aname(i) == '  H2A'))then
            st(i,1) = 191
          elseif((aname(i) == '  H2B'))then
            st(i,1) = 192
          elseif((aname(i) == '  C20'))then
            st(i,1) = 193
          elseif((aname(i) == '  H20'))then
            st(i,1) = 194
          elseif((aname(i) == '  C21'))then
            st(i,1) = 195
          elseif((aname(i) == ' H21A'))then
            st(i,1) = 196
          elseif((aname(i) == ' H21B'))then
            st(i,1) = 197
          elseif((aname(i) == ' H21C'))then
            st(i,1) = 198
          elseif((aname(i) == '  C22'))then
            st(i,1) = 199
          elseif((aname(i) == ' H22A'))then
            st(i,1) = 200
          elseif((aname(i) == ' H22B'))then
            st(i,1) = 201
          elseif((aname(i) == '  C23'))then
            st(i,1) = 202
          elseif((aname(i) == ' H23A'))then
            st(i,1) = 203
          elseif((aname(i) == ' H23B'))then
            st(i,1) = 204
          elseif((aname(i) == '  C24'))then
            st(i,1) = 205
          elseif((aname(i) == ' H24A'))then
            st(i,1) = 206
          elseif((aname(i) == ' H24B'))then
            st(i,1) = 207
          elseif((aname(i) == '  C25'))then
            st(i,1) = 208
          elseif((aname(i) == '  H25'))then
            st(i,1) = 209
          elseif((aname(i) == '  C26'))then
            st(i,1) = 210
          elseif((aname(i) == ' H26A'))then
            st(i,1) = 211
          elseif((aname(i) == ' H26B'))then
            st(i,1) = 212
          elseif((aname(i) == ' H26C'))then
            st(i,1) = 213
          elseif((aname(i) == '  C27'))then
            st(i,1) = 214
          elseif((aname(i) == ' H27A'))then
            st(i,1) = 215
          elseif((aname(i) == ' H27B'))then
            st(i,1) = 216
          elseif((aname(i) == ' H27C'))then
            st(i,1) = 217
          endif          


           

         endif
      enddo
      close(unit=1)

c ----- PSF ATOM FORMAT:
c     atomID, segmentID,resID,resname,atomname,atomtype,charge,mass
100   format(i8,a5,i5,a5,a5,a5,f12.6,f14.4)
c ----- PDB FORMAT:
c     PDB: ATOM, atomID,atomname,resname,resID,coordX,Y,Z,fixed,something,
c          empty space, segmentID
101   format(a4,i7,a5,a5,i5,f12.3,2f8.3,f6.2,f6.2,a5,a4)
102   format(a4,i7,a5,a5,i7,f10.3,2f8.3)
110   format(2a5,2f6.2)

      RETURN
      END SUBROUTINE

      ! --------------------------------------------------
      SUBROUTINE psfgen() ! generate psf structure
      ! --------------------------------------------------
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k, it,kt,ct,ct2
      INTEGER na,excl
      INTEGER p1,p2,p3,p4
      INTEGER ign 
      INTEGER ok,a(4) 
      REAL(8) cut,bd2
      REAL(8) cut2l,cut2r,rtol,ltol
      CHARACTER*4 buf4
      CHARACTER*5 buf5,t1,t2
      CHARACTER*6 buf6

      cut = 1.5 
      ltol = 1.61
      rtol = 0.19
      cut2l = (cut-ltol)**2.
      cut2r = (cut+rtol)**2.

      ntbo = 0
      nta = 0
      ntd = 0
      !--- find all possible bonds
      excl = 0
      do i=1,tsz
        do j=1,i
          bd2 = (ts(i,1)-ts(j,1))*(ts(i,1)-ts(j,1))+
     &          (ts(i,2)-ts(j,2))*(ts(i,2)-ts(j,2))+
     &          (ts(i,3)-ts(j,3))*(ts(i,3)-ts(j,3))
c            write(*,*) bd2
          if((cut2l < bd2).and.(bd2 < cut2r))then
            t1 = typtb(st(i,1))
            t2 = typtb(st(j,1))
            ign = 0
            do it=1,nexb
              if(((t1 == exb(it,1)).and.
     &          (t2 == exb(it,2))).or.
     &          ((t2 == exb(it,1)).and.
     &          (t1 == exb(it,2))))then
                ign = 1
                !write(*,*) t1,t2,exb(it,1),exb(it,2),ign
                excl = excl + 1
                goto 1050! important, otherwise exclusion lost in the
                         ! next iteration
              else
                ign = 0
              endif
            enddo
1050        if(ign == 0)then
              ntbo = ntbo+1
              tbo(ntbo,1) = j
              tbo(ntbo,2) = i
              tbo(ntbo,3) = 1
            endif
          endif
        enddo
      enddo
      write(*,*) " [psfgen] the number of excluded bonds:", excl

      !--- find all possible angles
      do i=1,ntbo
        p1 = tbo(i,1)
        p2 = tbo(i,2)
        do j=1,i
          if(i /=j)then
            if((tbo(j,1) == p1))then
              nta = nta + 1
              ta(nta,1) = tbo(j,2)
              ta(nta,2) = p1 
              ta(nta,3) = tbo(i,2)
              ta(nta,4) = 1 
              !write(*,'(i2,3i8)') 1,(ta(nta,k),k=1,3)
            endif
            if((tbo(j,2) == p1))then
              nta = nta + 1
              ta(nta,1) = tbo(j,1)
              ta(nta,2) = p1 
              ta(nta,3) = tbo(i,2)
              ta(nta,4) = 1 
              !write(*,'(i2,3i8)') 2,(ta(nta,k),k=1,3)
            endif
            if((tbo(j,1) == p2))then
              nta = nta + 1
              ta(nta,1) = tbo(j,2)
              ta(nta,2) = p2 
              ta(nta,3) = tbo(i,1)
              ta(nta,4) = 1 
              !write(*,'(i2,3i8)') 3,(ta(nta,k),k=1,3)
            endif
            if((tbo(j,2) == p2))then
              nta = nta + 1
              ta(nta,1) = tbo(j,1)
              ta(nta,2) = p2 
              ta(nta,3) = tbo(i,1)
              ta(nta,4) = 1 
              !write(*,'(i2,3i8)') 4,(ta(nta,k),k=1,3)
            endif
          endif
        enddo
      enddo

      !--- find all possible dihedrals ----
      do i=1,ntbo
        p1 = tbo(i,1)
        p2 = tbo(i,2)
        do j=1,nta
          ! --- attachment to the angle j from the left side
          if((p1 == ta(j,1)).and.(p2 /= ta(j,2))
     &      .and.(p2 /= ta(j,3)))then
            ign = 0
            p3 = ta(j,2) !middle point of angle
            p4 = ta(j,3)
            do it=1,ntd
              if((p2 == td(it,1)).and.(p1 == td(it,2))
     &           .and.(p3 == td(it,3)).and.(p4 == td(it,4)))then
                ign = 1
                goto 1000
              endif
              if((p4 == td(it,1)).and.(p3 == td(it,2))
     &           .and.(p1 == td(it,3)).and.(p2 == td(it,4)))then
                ign = 1
                goto 1000
              endif
            enddo
1000        if(ign /= 1)then
              ntd = ntd + 1
              td(ntd,1) = p2
              td(ntd,2) = p1 
              td(ntd,3) = p3 
              td(ntd,4) = p4
              td(ntd,5) = 1
            endif
            !write(*,'(a2,i2,4i8)') 'p1',1,(td(ntd,k),k=1,4) !show the dihedral
          endif
          ! --- attachment to the angle j from the left side
          if((p2 == ta(j,1)).and.(p1 /= ta(j,2))
     &       .and.(p1 /= ta(j,3)))then
            ign = 0
            p3 = ta(j,2) !middle point of angle 
            p4 = ta(j,3)
            do it=1,ntd
              if((p1 == td(it,1)).and.(p2 == td(it,2))
     &           .and.(p3 == td(it,3)).and.(p4 == td(it,4)))then
                ign = 1
                goto 1010
              endif
              if((p4 == td(it,1)).and.(p3 == td(it,2))
     &           .and.(p2 == td(it,3)).and.(p1 == td(it,4)))then
                ign = 1
                goto 1010
              endif
            enddo
1010        if(ign /= 1)then
              ntd = ntd + 1
              td(ntd,1) = p1
              td(ntd,2) = p2 
              td(ntd,3) = p3 
              td(ntd,4) = p4
              td(ntd,5) = 1
            endif
            !write(*,'(a2,i2,4i8)') 'p2',2,(td(ntd,k),k=1,4) !show the dihedral
          endif
          ! --- attachment to the angle j from the right side
          if((p1 == ta(j,3)).and.(p2 /= ta(j,2))
     &       .and.(p2 /= ta(j,1)))then
            ign = 0
            p3 = ta(j,1)
            p4 = ta(j,2)  !middle point of angle
            do it=1,ntd
              if((p3 == td(it,1)).and.(p4 == td(it,2))
     &           .and.(p1 == td(it,3)).and.(p2 == td(it,4)))then
                ign = 1
                goto 1020
              endif
              if((p2 == td(it,1)).and.(p1 == td(it,2))
     &           .and.(p4 == td(it,3)).and.(p3 == td(it,4)))then
                ign = 1
                goto 1020
              endif
            enddo
1020        if(ign /= 1)then
              ntd = ntd + 1
              td(ntd,1) = p3
              td(ntd,2) = p4 
              td(ntd,3) = p1 
              td(ntd,4) = p2
              td(ntd,5) = 1
            endif
            !write(*,'(a2,i2,4i8)') 'p1',3,(td(ntd,k),k=1,4) !show the dihedral
          endif
          ! --- attachment to the angle j from the right side
          if((p2 == ta(j,3)).and.(p1 /= ta(j,2))
     &       .and.(p1 /= ta(j,1)))then
            ign = 0
            p3 = ta(j,1)
            p4 = ta(j,2) !middle point of angle
            do it=1,ntd
              if((p3 == td(it,1)).and.(p4 == td(it,2))
     &           .and.(p2 == td(it,3)).and.(p1 == td(it,4)))then
                ign = 1
                goto 1030
              endif
              if((p1 == td(it,1)).and.(p2 == td(it,2))
     &           .and.(p4 == td(it,3)).and.(p3 == td(it,4)))then
                ign = 1
                goto 1030
              endif
            enddo
1030        if(ign /= 1)then
              ntd = ntd + 1
              td(ntd,1) = p3
              td(ntd,2) = p4 
              td(ntd,3) = p2 
              td(ntd,4) = p1
              td(ntd,5) = 1
            endif
            !write(*,'(a2,i2,4i8)') 'p2',4,(td(ntd,k),k=1,4) !show the dihedral
          endif
        enddo
      enddo

      RETURN
      END SUBROUTINE

      ! --------------------------------------------------
      SUBROUTINE savepsf !write out PSF structure
      ! --------------------------------------------------
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k, it,kt,ct,ct2
      INTEGER nb,nbp,temp,t1
      INTEGER boa(BIG),ana(BIG)
      INTEGER exat,exbo,exan,exdi 
      REAL(8) th,totq
      REAL(8) chartx(BIG) 
      CHARACTER*15 buf15

      exat = 0
      exbo = 0
      exan = 0
      exdi = 0
      totq = 0

      open(file='chol1.psf',unit=1)

      write(1,*) 'PSF'
      write(1,*) '1 !NTITLE'
      write(1,*) ' REMARKS original generated
     &              structure x-plor psf file'
      write(1,*) ' '
      write(1,'(i12,a8)') (ssz-exat),' !NATOM'

c      open(file='char.out',unit=2)
c      do i=1,Nm
c      read(2,'(a15,f8.5)') buf15,chartx(i)
c      enddo
c      close(unit=2)


      it = 0
      do i=1,ssz
        if(st(i,2) /= 0)then
          it = it + 1
          t1 = st(i,1)
          write(1,100)
     &    it,segtb(1),1,rname(i),aname(i),typtb(t1),
     &    chartb(t1),masstb(t1)
c     &    0.000000,masstb(t1)

c          totq = totq + chartb(i)
          totq =0.
        endif
      enddo
      write(*,*) " [savepsf] total charge:",totq

      write(1,*)
      write(1,'(i8,a15)') (bsz-exbo),'   !NBOND: bonds'
      j = 0
      do i=1,bsz
        if(bo(i,3) /= 0) then
          j = j+1
          boa(j) = i
        endif
      enddo
      write(*,*) '    bonds to write:',j
      write(1,'(8i8)') ((bo(boa(i),k),k=1,2),i=1,j)

      write(*,*) '   angles to write:',asz
      write(1,*)
      write(1,'(i8,a16)') (asz-exan),'  !NTHETA: angles'
      write(1,'(9i8)') ((an(i,k),k=1,3),i=1,asz)
      
      temp = 0
      do i=1,dsz
        if(di(i,5) == 1) then
          temp = temp + 1
        endif
      enddo
      write(*,*) 'dihedrals to write:',temp!dsz
      write(1,*)
      write(1,'(i8,a18)') (dsz),'   !NPHI: dihedrals'
      write(1,'(8i8)') ((di(i,k),k=1,4),i=1,dsz)
      write(1,*)

      write(1,*) '           0 !NIMPHI: impropers'
      write(1,*) '    0 !NDON: donors'
      write(1,*) '    0 !NACC: acceptors'
      write(1,*) '    0 !NNB'

      close(unit=1)

630     format(i8,a5,i2,a7,a4,a5,f13.6,f14.4,a12)

c ----- PSF ATOM FORMAT:
c     atomID, segmentID,resID,resname,atomname,atomtype,charge,mass
100   format(i8,a5,i5,a5,a5,a5,f12.6,f14.4)
c ----- PDB FORMAT:
c     PDB: ATOM, atomID,atomname,resname,resID,coordX,Y,Z,fixed,something,
c          empty space, segmentID
101   format(a4,i7,a5,a5,i5,f12.3,2f8.3,f6.2,f6.2,a5,a4)
110   format(2a5,2f6.2)

      RETURN
      END SUBROUTINE

      ! --------------------------------------------------
      SUBROUTINE savepdb() !write out PDB structure
      ! --------------------------------------------------
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k, it,kt,ct,ct2
      INTEGER nb,nbp,temp,t1
      INTEGER boa(BIG),ana(BIG)
      INTEGER exat,exbo,exan,exdi 
      REAL(8) th,totq,fix,rborder

      open(file='chol1.pdb',unit=1)

      write(1,'(a6,i8)') 'REMARK',ssz
      write(*,*) " [savepdb] saving",ssz,"atoms" 

      it = 0
      do i=1,ssz
c        if(st(i,2) /= 0)then
          it = it + 1
          t1 = st(i,1)
          write(1,101)
     &    'ATOM',it,aname(i),rname(i),3,(s(i,k),k=1,3),0.00,0.00
c        endif
      enddo

      close(unit=1)

c      open(file='recf.pdb',unit=1)
c
c      write(1,'(a6,i8)') 'REMARK',ssz
c      write(*,*) " [savepdb] saving",ssz,"atoms" 
c
c      it = 0
c      rborder = 20.5!9.89
c      do i=1,ssz
cc        if(st(i,2) /= 0)then
c          it = it + 1
c          t1 = st(i,1)
c          fix = 0.00
c          if((s(i,3)>rborder).and.(st(i,1)==1))then
c            fix = 1.00
c          endif
c          if((i>960).and.(i<1081))then
c            fix = 1.00
c          endif
c          write(1,101)
c     &    'ATOM',it,aname(i),rname(i),3,(s(i,k),k=1,3),fix,0.00
cc        endif
c      enddo
c
c      close(unit=1)

630     format(i8,a5,i2,a7,a4,a5,f13.6,f14.4,a12)

c ----- PSF ATOM FORMAT:
c     atomID, segmentID,resID,resname,atomname,atomtype,charge,mass
100   format(i8,a5,i5,a5,a5,a5,f12.6,f14.4)
c ----- PDB FORMAT:
c     PDB: ATOM, atomID,atomname,resname,resID,coordX,Y,Z,fixed,something,
c          empty space, segmentID
101   format(a4,i7,a5,a5,i5,f12.3,2f8.3,f6.2,f6.2,a5,a4)
110   format(2a5,2f6.2)

      RETURN
      END SUBROUTINE
      ! --------------------------------------------------
      SUBROUTINE waterbox() ! create waterbox 
      ! --------------------------------------------------
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k,ix,iy,iz,it
      INTEGER ci,ign
      INTEGER nw(3)
      REAL(8) wv(3,3),loc,loc2,r2
      REAL(8) latv(3),sp(3),v(3)
      REAL(8) bl,ba
      
      bl = 0.957
      ba = 104.5*d2r
      loc = 2.5
      loc2 = loc*loc

      nw(1) = 10 ! number of water in each dimension
      nw(2) = 10
      nw(3) = 10
      write(*,*) " [waterbox]",nw(1)*nw(2)*nw(3),"waters to add"

      latv(1) = 2.8 ! lattice vector
      latv(2) = 2.8
      latv(3) = 2.8
      
      do k=1,3
        sp(k) = -nw(k)*latv(k)/2. !starting point 
      enddo

      wv(1,1) = 0. 
      wv(1,2) = 0.
      wv(1,3) = 0.

      wv(2,1) = -bl*cos(ba/2.)
      wv(2,2) = +bl*sin(ba/2.)
      wv(2,3) = 0. 

      wv(3,1) = -bl*cos(ba/2.)
      wv(3,2) = -bl*sin(ba/2.)
      wv(3,3) = 0.

      cn = 0
      do ix=1,nw(1)
        do iy=1,nw(2)
          do iz=1,nw(3)

            v(1) = sp(1) + latv(1)*(ix-1)
            v(2) = sp(2) + latv(2)*(iy-1)
            v(3) = sp(3) + latv(3)*(iz-1)

            ign = 0
            do i=1,ssz
              r2 = 0.
              do k=1,3
                r2 = r2 + (v(k)-s(i,k))**2.
              enddo
              if(r2 < loc2)then
                ign = 1
                goto 1060
              endif
            enddo

1060        if(ign == 0)then
              ci = ssz+cn
              do it=1,3
                do k=1,3
                  s(ci+it,k) = v(k) + wv(it,k) 
                enddo
                !write(*,*) ci,ci+it,(s(ci+it,k),k=1,3)
                st(ci+it,1) = 43+it
                st(ci+it,2) = 1 
                rname(ci+it) =' IP3 '
              enddo
              aname(ci+1) = '  O  '
              aname(ci+2) = ' H1  '
              aname(ci+3) = ' H2  '
              cn = cn + 3
              bo(bsz+1,1) = ci+1
              bo(bsz+1,2) = ci+2
              bo(bsz+1,3) = 1
              bo(bsz+2,1) = ci+1
              bo(bsz+2,2) = ci+3
              bo(bsz+2,3) = 1
              bsz = bsz + 2
              an(asz+1,1) = ci+2
              an(asz+1,2) = ci+1
              an(asz+1,3) = ci+3
              asz = asz + 1
            endif

          enddo
        enddo
      enddo

      write(*,*) " [waterbox] added",cn," water atoms"

      ssz = ssz + cn

      RETURN
      END SUBROUTINE


      ! --------------------------------------------------
      SUBROUTINE testpsf !test psfgen
      ! --------------------------------------------------
      USE MM
      IMPLICIT NONE
      ! - - - local declarations - - -
      INTEGER i,j,k, it,kt,ct,ct2
      INTEGER nb,nbp,test 
      INTEGER boa(BIG),ana(BIG)
      INTEGER exat,exbo,exan,exdi 

      test = 3

      if(test == 1) then
        tsz = 2
        do k=1,3
          ts(1,k) = 0.
          ts(2,k) = 0.
          ts(3,k) = 0.
        enddo

        ts(1,1) = 1.5
        ts(3,1) = -1.5
      elseif(test == 2)then
        tsz = 3
        do k=1,3
          ts(1,k) = 0.
          ts(2,k) = 0.
          ts(3,k) = 0.
        enddo

        ts(1,1) = -0.75
        ts(2,1) =  0.75
        ts(3,2) = -0.95
      elseif(test == 3)then
        tsz = 4
        do k=1,3
          ts(1,k) = 0.
          ts(2,k) = 0.
          ts(3,k) = 0.
          ts(4,k) = 0.
        enddo

        ts(1,1) = -1.0
        ts(2,1) =  1.0
        ts(3,2) =  1.0
        ts(4,2) = -1.0
      endif

      RETURN
      END SUBROUTINE



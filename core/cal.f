      program cal
      real::edge,area,volume

      edge=41.113468

      area=20*(sqrt(3.0)*edge**2/4)

      volume=(15+5*sqrt(5.0))*edge**3/12

      write(*,*) "edge length: ",edge,", surface area: ",area,
     &", volume: ",volume

      end

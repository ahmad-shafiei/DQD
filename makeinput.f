      program makeinput
c
      implicit none
c     iIsolem is the maximum number of islands;
c     iEnlevm is the maximum number of levels per island
      integer iIsolem,iEnlevm
      parameter (iIsolem=6,iEnlevm=30)
c     iPuntatori are the pointers to the top occupied level: the first 
c     index refers to the island, the second to the spin orientation;
c     iIsole is the number of islands;
c     iEnlev is the number of levels per island
      integer iPuntatori(iIsolem,2),iNn(iIsolem),i,iIsole,j,iEnlev,
     #nx,ny,l,imin
c     dEpsilon is the energy for each level
      double precision pi,m,hbar,wx,wy,dEpsilon(iEnlevm**2),emin,eauxi
c
      iIsole=2
      iEnlev=20
c
      open(unit=10,file='enlev.dat',status='unknown')
      write(10,*) iEnlev
      close(10)
c
      if (iIsole.gt.iIsolem) then
         write(*,*) 'iIsole>iIsolem'
         goto 10
      end if
      if (iEnlev.gt.iEnlevm) then
         write(*,*) 'iEnlev>iEnlevm'
         goto 10
      end if
c
C     The occupancy are initialized at arbitrary values
      do i=1,iIsole
         iPuntatori(i,1)=10
         iPuntatori(i,2)=10
c        puntatori al massimo livello occupato nell'isola i
c        dagli elettroni con spin up (i,1) e down (i,2)
         if (iPuntatori(i,1).le.0) then
            write(*,*) 'iPuntatori(',i,'1)<=0'
            goto 10
         end if
         if (iPuntatori(i,1).ge.iEnlev) then
            write(*,*) 'iPuntatori(',i,'1)>=',iEnlev
            goto 10
         end if
         if (iPuntatori(i,2).le.0) then
            write(*,*) 'iPuntatori(',i,'2)<=0'
            goto 10
         end if
         if (iPuntatori(i,2).ge.iEnlev) then
            write(*,*) 'iPuntatori(',i,'2)>=',iEnlev
            goto 10
         end if
         iNn(i)=-iPuntatori(i,1)-iPuntatori(i,2)
c        cariche iniziali sulle isole (con loro segno) espresse in
c        termini di cariche elementari ("-" perche' si tratta di elettroni)
      end do
c
      open(unit=11,file='puntatori.dat',status='unknown')
      do i=1,iIsole
         write(11,*) i,1,iPuntatori(i,1)
c        puntatore dello spin up nell'isola i
         write(11,*) i,0,iPuntatori(i,2)
c        puntatore dello spin down nell'isola i
      end do
      close(11)
c
      open(unit=14,file='caricheiniz.dat',status='unknown')
      do i=1,iIsole
         write(14,*) i,iNn(i)
      end do
      close(14)

      open(unit=12,file='tabellaoccupazioni.dat',status='unknown')
      do i=1,iIsole
         do j=1,iEnlev

            if (j.le.iPuntatori(i,1)) then
               write(12,*) i,j,1,1
c              isola, livello energia, spin up(1) o down(0),se c'e'(1) o no(0)
            else
               write(12,*) i,j,1,0
            end if

            if (j.le.iPuntatori(i,2)) then
               write(12,*) i,j,0,1
           else
               write(12,*) i,j,0,0
            end if

         end do
         write(12,*) ''
      end do
      close(12)

      pi=acos(-1.0d0)
      m=6.103291099d-32
c     massa dell'elettrone (in Kg)
      hbar=1.05457266d-34
c     hbar (in J)
      wx=20.0d-9
      wy=30.0d-9
c     dimensioni (in m) lungo x e y dei dot di forma rettangolare
c
      i=0
      open(unit=33,file='barab.dat',status='unknown')
      do nx=1,iEnlev
         do ny=1,iEnlev
            i=i+1
            dEpsilon(i)=((hbar*pi)**2/(2*m))*((nx/wx)**2+(ny/wy)**2)
            write(33,*) dEpsilon(i) 
         end do
      end do
      close(33)
c     ordinamento crescente delle energie
      i=0
 11   i=i+1
      imin=i
      emin=dEpsilon(i)
      do j=i,iEnlev**2
         if (dEpsilon(j).lt.emin) then
            imin=j
            emin=dEpsilon(j)
         end if
      end do
      eauxi=dEpsilon(i)
      dEpsilon(i)=dEpsilon(imin)
      dEpsilon(imin)=eauxi
      if (i.lt.iEnlev) goto 11
c     perche', pur avendo calcolato iEnlev**2 livelli, ci bastano
c     gli iEnlev livelli piu' bassi tra questi iEnlev**2
c
c
      open(unit=13,file='tabellaenergie.dat',status='unknown')
      do i=1,iIsole
         do j=1,iEnlev
            do l=1,2
               write(13,*) i,j,l,dEpsilon(j)
c               write(13,*) i,j,l,0.0d0
c              energie (in Joule)
            end do
         end do
         write(13,*) ''
      end do
      close(13)
c
      open(unit=14,file='eminmax.dat',status='unknown')
      write(14,*) 0.0d0
c     limite inferiore del range di energie considerato
      write(14,*) dEpsilon(iEnlev)+1.0d-19
c     limite superiore del range di energie considerato
      close(14)
c
 10   end


  



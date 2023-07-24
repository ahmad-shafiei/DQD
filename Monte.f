      program simulazione
*------------------------------------
*     Dichiarazioni
*------------------------------------
      implicit double precision(d)
      implicit integer(i)
c     dPotenziali is a temporary matrix for the storage of the 
c     elastance matrix 
      dimension dPotenziali(6,6)
c@@@@
c     dCapacita stores the capacitance matrix and, later, the 
c     elastance matrix
      dimension dCapacita(6,6),dCostante(6),dEsterno(10),iTunnel(6,2)
      dimension dCollegamenti(6,16),dBias(6),iNn(6)
      
c     dMassProb is the probability of a transition from a state
c     (i,j,l) to (i,j+1,l)
      dimension dMassProb(2,6),dMxtmp(2,6)
      
      dimension dQmedia(6),dImedia(6)
      dimension dQProvv(6,300000),dIProvv(6,300000)
c     dVoltaggio is the voltage that conected to a node. 
c     teh initial (1) or the final vlaue (2) of the voltage. 
c     the second index for the number of external contacts. 
      dimension dVoltaggio(2,10),dRampa(10)

      dimension dRate(6,2,2,30,30)
      dimension iTab(6,30,2)
      dimension dEpsilon(6,30,2)

c      istodep is the accumulation array for the occupancy histograms
c      the first index is the voltage step, the second is the island, 
c      the third is the energy level, and the fourth is the spin
      dimension istodep(3000,30,100,2)
c      deposon is the array for the occupancy histograms normalized
c      dividing by iMedie
      dimension deposon(3000,30,100,2)
      common /satrapo/isatrapon,istodep
      common /primo/ dCollegamenti
      common /secondo/ dEsterno
      common /terzo/dCapacita
      common /bisterzo/dCostante
      common /quarto/iIsole,iNesterni
      common /quinto/dBias
      common /sesto/ iNn,dBeta
      common /ottavo/ dRate

      common /decimo/ iSeme
      common /undicesimo/ iNtunnel
      common /dodicesimo/ iTunnel
      common /quindicesimo/ dQProvv,dIProvv,dRampa,iPassi
      common /sedicesimo/ dLarghezza
      common /diciassettesimo/dVoltaggio
      common /diciottesimo/ iMedie
      common /ventesimo/iEnlev
      common /ventunesimo/iTab

      common /ventitreesimo/dEpsilon
      common /ventiquattresimo/iCont3
      common /ventottesimo/dEmin,dEmax
      dimension dtensio(10,10000)
*------------------------------------
*     Configurazione
*------------------------------------
      namelist/data/iIsole,iNesterni,iPassi,iMedie,dLarghezza,
     $     dCollegamenti,iNtunnel,iTunnel,
     $     dCostante,dVoltaggio,iNn,dTemperatura,iSeme

      character*20 Fcircuito
      print*, "Nome del file di circuito (mx 20 caratteri)"
      read*,Fcircuito
      print*,Fcircuito	
      
      open(22,file=Fcircuito,status="old")
      read(22,data)
      close(22)
        
      print*,imedie,iseme      
      dBeta=1.855d3/dTemperatura
      print*, "beta",dBeta

c     iEnlev is the number of energy levels in each dot
      open(unit=51,file='enlev.dat',status='unknown')
      read(51,*) iEnlev
      close(51)
      if (iEnlev.gt.30) then
         write(*,*) 'iEnlev>30'
         goto 99
      end if

      open(unit=54,file='tabellaenergie.dat',status='unknown')
      open(unit=154,file='tabellaenergie2.dat',status='unknown')
      do i=1,iIsole
         do j=1,iEnlev
            do l=1,2
               read(54,*) iaux1,iaux2,iaux3,dEpsilon(i,j,l)
               dEpsilon(i,j,l)=dEpsilon(i,j,l)/2.56d-20
               write(154,*) iaux1,iaux2,iaux3,dEpsilon(i,j,l)
            end do
         end do
      end do
      close(54)
      close(154)
      open(unit=56,file='eminmax.dat',status='unknown')
      read(56,*) dEmin
      read(56,*) dEmax
      dEmax=dEmax/(2.56d-20)
      dEmin=dEmin/(2.56d-20)
      close(56)
*------------------------------------
*     Inizializzazioni
*------------------------------------
      do i=1,iNesterni
         if (iPassi.ne.1) then
            dRampa(i+iIsole)=((dVoltaggio(2,i+Iisole)/1.6d-1) 
     $        - (dVoltaggio(1,i+Iisole)/1.6d-1))/(iPassi-1)
         else
            dRampa(i+iIsole)=0.0d0
         end if
      enddo
      write(6,*) 'dVolt', dVoltaggio(1,3), dVoltaggio(2,3)
      write(6,*) 'dVolt', dVoltaggio(1,4), dVoltaggio(2,4)
      write(6,*) 'dVolt', dVoltaggio(1,5), dVoltaggio(2,5)  !###################
      write(6,*) 'dVolt', dVoltaggio(1,6), dVoltaggio(2,6)  !###################
      write(6,*) 'dRampa', dRampa(3), dRampa(4)  
      write(6,*) 'dRampa', dRampa(5),dRampa(6)              !###############
      
      do i=1,iPassi
        do j=1,iIsole
          do k=1,iEnlev
            do l=1,2
              istodep(i,j,k,l)=0
            end do
          end do
        end do
      end do 
*-----------------------------------------
*     Calcolo della matrice delle capacita
*-----------------------------------------
      call MatriceCapacita(dCapacita)
c     the vedo subroutine writes the capacitance matrix to the standard out
      call vedo(iIsole,dCapacita)
*-----------------------------------------
*     Calcolo della matrice dei potenziali
*-----------------------------------------
C     Calculation of the elastance matrix      
      call inversione(iIsole,dCapacita,6,dPotenziali)
      call vedo(iIsole,dPotenziali)
      do i=1,iIsole
         do j=1,iIsole
           dCapacita(i,j)=dPotenziali(i,j)
         enddo
      enddo	
      call vedo(iIsole,dCapacita)
*------------------------------------
*------------------------------------      
      isatrapon=0
      open(unit=77,file='volt.dat',status='unknown')
*------------------------------------
*     Singola Configurazione
*------------------------------------
      do iCont3=1,iPassi
         do ig=1,iNesterni
            dEsterno(ig) = dVoltaggio(1,ig+iIsole)/1.6d-1+
     $       dRampa(ig+iIsole)*(iCont3-1)
             dtensio(ig,iCont3)=dEsterno(ig)       
         enddo
         open(unit=55,file='caricheiniz.dat',status='unknown')
         do i=1,iIsole
            read(55,*) iaux1,iNn(i)
c           cariche iniziali sulle isole (con loro segno) espresse in
c           termini di cariche elementari
         end do
         close(55)
c        iTab is the occupancy array: the first index is the Island number
c        the second index is the energy level number, the third index is the
c        spin
         open(unit=53,file='tabellaoccupazioni.dat',status='unknown')
         do i=1,iIsole
            do j=1,iEnlev
               do l=1,2
                  read(53,*) iaux1,iaux2,iaux3,iTab(i,j,l)
               end do
            end do
         end do
         close(53)
         do i=1,iIsole
            do l=1,2
               do j=1,iEnlev-1
                  ipre=iTab(i,j,l)
                  ipost=iTab(i,j+1,l)
                  if (ipre.lt.ipost) then
c                    (i.e., ipre=0 and ipost=1)
                     write(*,*) 'error: the initial state is excited!'
                     goto 99
                  end if
               end do
            end do
         end do
C        q tilde is computed                 
         call bias(dBias)
         call montecarlo(*99,dQmedia,dImedia)         
         do i=1,iIsole
            dQProvv(i,iCont3) = dQmedia(i)
         enddo         
         do i=1,iNtunnel
            dIProvv(i,iCont3) = dImedia(i)
         enddo
      enddo
      open(unit=97,file='istogram.dat',status='unknown')
      write(6,*) 'iPassi= ',iPassi
      write(6,*) 'iIsole= ',iIsole
      write(6,*) 'iEnlev= ',iEnlev
      do i=1,iPassi
        do j=1,iIsole
          do k=1,IEnlev
            do l=1,2
              deposon(i,j,k,l)=float(istodep(i,j,k,l))/iMedie
            end do 
          end do
        end do
      end do
      do j=1,iIsole
        do i=1,iPassi
          write(97,*) ((deposon(i,j,k,l), l=1,2), k=1,iEnlev)
        end do 
      end do
      close(97)      
      call scrividati
      do ikon=1,iPassi
        write(77,*) ikon, dtensio(1,ikon), dtensio(2,ikon)
      end do
      close(77)
 99   end      
*------------------------------------
*     Subroutine per il calcolo
*     della matrice delle capacita
*------------------------------------      
c      subroutine MatriceCapacita(dCapacita,dCapacita1)
      subroutine MatriceCapacita(dCapacita)     
      implicit double precision(d)
      implicit integer(i)      
      dimension dCollegamenti(6,16)
      common /primo/dCollegamenti
      common /quarto/ iIsole,iNesterni
      dimension dCapacita(6,6)
      do i=1,iIsole
         dCapacita(i,i)=0d0
         do j=1,(iIsole+iNesterni)
            dCapacita(i,i)=dCapacita(i,i)+dCollegamenti(i,j)
         enddo
      enddo
      do i=1,iIsole
         do j=1,iIsole
            if(i.ne.j) then
               dCapacita(i,j)=-dCollegamenti(i,j)
            endif
         enddo
      enddo  
      end      
*-------------------------------------
* Visualizzazione di una matrice
*-------------------------------------
      subroutine vedo(iNn,dMatrice)      
      implicit double precision(d)
      implicit integer(i)      
      dimension dMatrice(6,6)
 10   format(6(1x,f10.4))
      do i=1,iNn
         write(*,10) (dMatrice(i,j),j=1,iNn)
      enddo
      end
*-------------------------------------
* Subroutine per le cariche di bias
*-------------------------------------
      subroutine bias(dBias)
      implicit double precision(d)
      implicit integer(i)      
      dimension dCollegamenti(6,16)
      dimension dEsterno(10)
      dimension dBias(6)
      common /primo/ dCollegamenti
      common /quarto/ iIsole,iNesterni
      common /secondo/ dEsterno      
      do i=1,iIsole
         dBias(i) = 0 
         do j=1,iNesterni
            dBias(i) = dBias(i)+dCollegamenti(i,iIsole+j)*dEsterno(j)
        enddo
      enddo      
      end      
*------------------------------------------
*     Probabilita' di transizione
*------------------------------------------
      subroutine probabilita
      
      implicit double precision(d)
      implicit integer(i)
      common /satrapo/isatrapon,istodep
      common /secondo/dEsterno
      common /terzo/dCapacita
      common /bisterzo/ dCostante
      common /quarto/ iIsole,iNesterni
      common /quinto/dBias
      common /sesto/iNn,dBeta
      common /ottavo/ dRate
      common /undicesimo/ iNtunnel
      common /dodicesimo/ iTunnel
      common /sedicesimo/ dLarghezza
      common /ventesimo/iEnlev
      common /ventunesimo/iTab
      common /ventitreesimo/dEpsilon
      common /ventottesimo/dEmin,dEmax
      dimension dCapacita(6,6),dCostante(6),iNn(6),dBias(6),iM(6)
      dimension dRate(6,2,2,30,30)
      dimension iTab(6,30,2)
      dimension dEpsilon(6,30,2)
      dimension iTunnel(6,2),dEsterno(10)
      dimension istodep(3000,30,100,2)
      call energia(dEnergia,iNn)
      do i=1,iNtunnel
         iK=iTunnel(i,1)
         do j=1,2
c        transizione forward (1) o backward (2) di un elettrone
            if (iTunnel(i,2).le.iIsole) then
c           secondo nodo isola
               iP=iTunnel(i,2)
               if (j.eq.1) then
c              transizione forward da isola a isola
                  do iii=1,iIsole
                     iM(iii) = iNn(iii)
                  enddo
                  iM(iK) = iNn(iK) + 1
                  iM(iP) = iNn(iP) - 1
                  call energia(dEnew,iM)
                  dDenergia=dEnergia-dEnew
                  do l=1,2
c                 elettrone con spin up (1) o down (2)
                     do ne1=1,iEnlev
c                    indice del livello di energia nel primo nodo
                        if (iTab(iK,ne1,l).eq.0) then
                           do ne2=1,iEnlev
c                          indice del livello di energia nel secondo nodo
                              dRate(i,j,l,ne1,ne2)=0.0d0
                           end do
                        else
                           dEin=dEpsilon(iK,ne1,l)
                           dEfin=dEin+dDenergia
                           if ((dEfin.lt.dEmin).or.(dEfin.gt.dEmax))
     #                        then
                              iEnear=0
c                             dEfin fuori del range considerato
                           else
                              dDiff=abs(dEfin-dEpsilon(iP,1,l))
                              iEnear=1
c                             indice del livello piu' vicino a dEfin
                              do ne2=2,iEnlev
c                             indice del livello di energia nel secondo nodo
                                 dDiffnew=abs(dEfin-dEpsilon(iP,ne2,l))
                                 if (dDiffnew.lt.dDiff) then
                                    dDiff=dDiffnew
                                    iEnear=ne2
                                 end if
                              end do
                           end if
                           do ne2=1,iEnlev
c                          indice del livello di energia nel secondo nodo
                              if ((ne2.ne.iEnear).or.
     #                            (iTab(iP,ne2,l).eq.1)) then
                                 dRate(i,j,l,ne1,ne2)=0.0d0
                              else
                                 dRate(i,j,l,ne1,ne2)=dCostante(i)
                              end if
                           end do
                        end if
                     end do
                  end do
               else
c              transizione backward da isola a isola

                  do iii=1,iIsole
                     iM(iii) = iNn(iii)
                  enddo
                  iM(iK) = iNn(iK) - 1
                  iM(iP) = iNn(iP) + 1
                  call energia(dEnew,iM)
                  dDenergia=dEnergia-dEnew
                  do l=1,2
c                 elettrone con spin up (1) o down (2)
                     do ne2=1,iEnlev
c                    indice del livello di energia nel secondo nodo
                        if (iTab(iP,ne2,l).eq.0) then
                           do ne1=1,iEnlev
c                          indice del livello di energia nel primo nodo
                              dRate(i,j,l,ne1,ne2)=0.0d0
                           end do
                        else
                           dEin=dEpsilon(iP,ne2,l)
                           dEfin=dEin+dDenergia
                           if ((dEfin.lt.dEmin).or.(dEfin.gt.dEmax))
     #                        then
                              iEnear=0
c                             dEfin fuori del range considerato
                           else
                              dDiff=abs(dEfin-dEpsilon(iK,1,l))
                              iEnear=1
c                             indice del livello piu' vicino a dEfin
                              do ne1=2,iEnlev
c                             indice del livello di energia nel primo nodo
                                 dDiffnew=abs(dEfin-dEpsilon(iK,ne1,l))
                                 if (dDiffnew.lt.dDiff) then
                                    dDiff=dDiffnew
                                    iEnear=ne1
                                 end if
                              end do
                           end if
                           do ne1=1,iEnlev
c                          indice del livello di energia nel primo nodo
                              if ((ne1.ne.iEnear).or.
     #                            (iTab(iK,ne1,l).eq.1)) then
                                 dRate(i,j,l,ne1,ne2)=0.0d0
                              else
                                 dRate(i,j,l,ne1,ne2)=dCostante(i)
                              end if
                           end do
                        end if
                     end do
                  end do
               end if
            else
c           secondo nodo lead
               iP=iTunnel(i,2)-iIsole            
               if (j.eq.1) then
c              transizione forward da isola a lead
                  do iii=1,iIsole
                     iM(iii) = iNn(iii)
                  enddo
                  iM(iK) = iNn(iK) + 1
                  call energia(dTmpEnew,iM)
                  dEnew = dTmpEnew - dEsterno(iP)
                  dDenergia=dEnergia-dEnew
                  do l=1,2
c                 elettrone con spin up (1) o down (2)
                     do ne1=1,iEnlev
c                    indice del livello di energia nel primo nodo (isola)
                        if (iTab(iK,ne1,l).eq.0) then
                           dRate(i,j,l,ne1,1)=0.0d0
                        else
                           dEin=dEpsilon(iK,ne1,l)
                           dEfin=dEin+dDenergia
                           dfun=1/(1+exp(dBeta*dEfin))
                           dRate(i,j,l,ne1,1)=dCostante(i)*(1-dfun)
                        end if
                     end do
                  end do
               else
c              transizione backward da lead a isola
                  do iii=1,iIsole
                     iM(iii) = iNn(iii)
                  enddo
                  iM(iK) = iNn(iK) - 1
                  call energia(dTmpEnew,iM)
                  dEnew =dTmpEnew + dEsterno(iP)
                  dDenergia=dEnergia-dEnew
                  do l=1,2
c                 elettrone con spin up (1) o down (2)
                     do ne1=1,iEnlev
c                    indice del livello di energia nel primo nodo (isola)
                        if (iTab(iK,ne1,l).eq.1) then
                           dRate(i,j,l,ne1,1)=0.0d0
                        else
                           dEfin=dEpsilon(iK,ne1,l)
                           dEin=dEfin-dDenergia
                           dfun=1/(1+exp(dBeta*dEin))
                           dRate(i,j,l,ne1,1)=dCostante(i)*dfun
                        end if
                     end do
                  end do
               end if
            end if
         end do
      end do
      end      
*-------------------------------
*     Calcolo secco dell'energia
*-------------------------------
      subroutine energia(dEnergia,iNn)      
      implicit double precision(d)
      implicit integer(i)
      dimension dCapacita(6,6),iNn(6),dBias(6)
      common /terzo/dCapacita
      common /quarto/ iIsole,iNesterni
      common /quinto/dBias
      dEnergia = 0d0
      do i=1,iIsole
         do j=1,iIsole
            dEnergia=dEnergia+ 
     $   5d-1*(iNn(i)+dBias(i))*dCapacita(i,j)*(iNn(j)+dBias(j))
         enddo
      enddo
      end
*---------------------------------------------
*     Controllo della probabilita massima
*     evento
*-------------------
      subroutine montecarlo(*,dQmedia,dImedia)      
      implicit double precision(d)
      implicit integer(i)
      common /satrapo/isatrapon,istodep
      common /quarto/iIsole,iNesterni
      common /sesto/ iNn,dBeta
      common /ottavo/ dRate
      common /decimo/ iSeme
      common /undicesimo/ iNtunnel
      common /dodicesimo/ iTunnel
      common /sedicesimo/ dLarghezza
      common /diciottesimo/ iMedie
      common /ventesimo/iEnlev
      common /ventunesimo/iTab
      common /ventitreesimo/dEpsilon
      common /ventiquattresimo/iCont3
      dimension iTunnel(6,2),iNn(6),iCorrente(6)
      dimension dTempi(300000),dRate(6,2,2,30,30),dQdelta(6,300000)
      dimension iTab(6,30,2)
      dimension dEpsilon(6,30,2)
      dimension dQmedia(6),dImedia(6),dQsomma(6)
      dimension iChOld(6)
       dimension istodep(3000,30,100,2)
      do i=1,iIsole
         dQsomma(i) =0d0
      enddo
      do i=1,iNtunnel
         iCorrente(i) = 0
      enddo      
      dTempoTotale=0d0      
      do iCont2=1,iMedie
         call probabilita
         dLarghezza = 0d0
         do i=1,iNtunnel
c        indice della giunzione tunnel
           if (iTunnel(i,2).le.iIsole) then
             ne2max=iEnlev
           else
             ne2max=1
c            se il secondo nodo e' un terminale esterno considero
c            un solo rate (nel terminale esterno non considero piu'
c            livelli di energia)
            end if
            do j=1,2
c           transizione forward (1) o backward (2) di un elettrone
              do l=1,2
c             elettrone con spin up (1) o down (2)
                do ne1=1,iEnlev
c               indice del livello energetico nel primo nodo
                  do ne2=1,ne2max
c                 indice del livello energetico nel secondo nodo
                    dLarghezza = dLarghezza + dRate(i,j,l,ne1,ne2)
                  end do
                end do
              end do
            enddo
          enddo
         d1Random = dRan3(iSeme)
         d2Random = dRan3(iSeme)*dLarghezza        
         do i=1,iIsole
          iChOld(i)=iNn(i)
         enddo
         if (dLarghezza.eq.0.0d0) then
            dLarghezza=1.0d-6
         end if
         dTempi(iCont2) = - log(d1Random)/dLarghezza
         dTempoTotale =dTempoTotale + dTempi(iCont2)
         dLivello = 0d0
         do i=1,iNtunnel
            do j=1,2
               do l=1,2
                  do ne1=1,iEnlev
                     do ne2=1,iEnlev
                        dLivello = dLivello + dRate(i,j,l,ne1,ne2)
                        if (d2Random.lt.dLivello) then
                           if (j.eq.1) then
                              iCorrente(i)=iCorrente(i)+1
                           else
                              iCorrente(i)=iCorrente(i)-1
                           endif
                           iK=iTunnel(i,1)               
                           if (iTunnel(i,2).le.iIsole) then
                              iP=iTunnel(i,2)
                              if (j.eq.1) then
                                 iNn(iK)=iNn(iK)+1
                                 iNn(iP)=iNn(iP)-1
                                 iTab(iK,ne1,l)=0
                                 iTab(iP,ne2,l)=1
                              else
                                 iNn(iK)=iNn(iK)-1
                                 iNn(iP)=iNn(iP)+1
                                 iTab(iK,ne1,l)=1
                                 iTab(iP,ne2,l)=0
                              endif
                           else
                              if (j.eq.1) then
                                 iNn(iK)=iNn(iK)+1
                                 iTab(iK,ne1,l)=0
                              else
                                 iNn(iK)=iNn(iK)-1
                                 iTab(iK,ne1,l)=1
                              endif
                           endif
c                          relaxation:
                           inum=0
c                          number of occupied energy levels for isle iK and spin l
                           do iene=1,iEnlev
                              inum=inum+iTab(iK,iene,l)
                           end do
                           do iene=1,inum
                              iTab(iK,iene,l)=1
                           end do
                           do iene=inum+1,iEnlev
                              iTab(iK,iene,l)=0
                           end do
                           if (iTunnel(i,2).le.iIsole) then
                              inum=0
c                             number of occupied energy levels for isle iP and spin l
                              do iene=1,iEnlev
                                 inum=inum+iTab(iP,iene,l)
                              end do
                              do iene=1,inum
                                 iTab(iP,iene,l)=1
                              end do
                              do iene=inum+1,iEnlev
                                 iTab(iP,iene,l)=0
                              end do
                           end if
                           goto 30
                        endif
                     end do
                  end do
               end do
            enddo
         enddo
 30      continue
         do i=1,iIsole
           do j=1,iEnlev
             do l=1,2
               istodep(iCont3,i,j,l)=istodep(iCont3,i,j,l)+
     #         iTab(i,j,l)
             end do
           end do
         end do
         do j=1,iIsole
            dQdelta(j,iCont2) = iChOld(j)*dTempi(iCont2) 
         enddo
      enddo
C     do over iMedie completed
 50   do j=1,iIsole
         do i=1,iMedie
            dQsomma(j)=dQsomma(j)+dQdelta(j,i)
         enddo
         dQmedia(j)=dQsomma(j)/dTempoTotale
      enddo
      do i=1,iNtunnel
         dImedia(i) = iCorrente(i)/dTempoTotale
      enddo
      end
*----------------------------------
*     Scrittura file di uscita
*----------------------------------
      subroutine scrividati      
      implicit double precision(d)
      implicit integer(i)
      dimension dQProvv(6,300000),dIProvv(6,300000)
      dimension dRampa(10),dVoltaggio(2,10)
      common /quarto/iIsole,iNesterni
      common /undicesimo/ iNtunnel
      common /quindicesimo/ dQProvv,dIProvv,dRampa,iPassi
      common /diciassettesimo/dVoltaggio
      open(1,file="correnti1.dat",status="unknown")
      open(2,file="cariche.dat",status="unknown")
      open(3,file="riferimento1.dat",status="unknown")
  10   format(3(e14.4))
      do i=1,iPassi
         write(2,*) (dQProvv(ik,i),ik=1,iIsole)
         write(1,*) (dIProvv(ik,i)*1.6e-6,ik=1,3)
         write(3,*) (dVoltaggio(1,ik+iIsole)+(i-1)*dRampa(iIsole+ik)
     #   *1.6e-1,ik=1,4)
      enddo
      close(1)
      close(2)
      close(3)
      end
*---------------------------------------------
*     subroutine di inversione che si trova
*     sul numerical recipes
*---------------------------------------------
      subroutine inversione(iN,dA,iNfisso,dY)      
      implicit double precision(d)
      implicit integer(i)       
      dimension dA(iNfisso,iNfisso),dY(iNfisso,iNfisso),iIndx(20)
      do i=1,iN
         do j=1,iN
            dY(i,j)=0d0
         enddo
         dY(i,i)=1d0
      enddo
      call ludcmp(dA,iN,6,iIndx,dD)
      do j=1,iN
         call lubksb(dA,iN,6,iIndx,dY(1,j))
      enddo
      end
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      implicit double precision(a-h,o-z)
      implicit integer(i-n)      
      DIMENSION A(NP,NP),INDX(N),B(N)      
      II=0
      DO 12 I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0)THEN
            DO 11 J=II,I-1
               SUM=SUM-A(I,J)*B(J)
 11         CONTINUE
         ELSE IF (SUM.NE.0.) THEN
            II=I
         ENDIF
         B(I)=SUM
 12   CONTINUE
      DO 14 I=N,1,-1
         SUM=B(I)
         DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
 13      CONTINUE
         B(I)=SUM/A(I,I)
 14   CONTINUE
      RETURN
      END      
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)      
      implicit double precision(a-h,o-z)
      implicit integer(i-n)      
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)      
      D=1.
      DO 12 I=1,N
         AAMAX=0.
         DO 11 J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
 11      CONTINUE
         IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
         VV(I)=1./AAMAX
 12   CONTINUE
      DO 19 J=1,N
         DO 14 I=1,J-1
            SUM=A(I,J)
            DO 13 K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
 13         CONTINUE
            A(I,J)=SUM
 14      CONTINUE
         AAMAX=0.
         DO 16 I=J,N
            SUM=A(I,J)
            DO 15 K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
 15         CONTINUE
            A(I,J)=SUM
            DUM=VV(I)*ABS(SUM)
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
 16      CONTINUE
         IF (J.NE.IMAX)THEN
            DO 17 K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
 17         CONTINUE
            D=-D
            VV(IMAX)=VV(J)
         ENDIF
         INDX(J)=IMAX
         IF(A(J,J).EQ.0.)A(J,J)=TINY
         IF(J.NE.N)THEN
            DUM=1./A(J,J)
            DO 18 I=J+1,N
               A(I,J)=A(I,J)*DUM
 18         CONTINUE
         ENDIF
 19   CONTINUE
      RETURN
      END      
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: dRan3.f  
C TYPE   : function
C PURPOSE: generate random numbers
C I/O    :
C VERSION: 27 MAI 94
C COMMENT: Initialize idum with negative integer
C=========+=========+=========+=========+=========+=========+=========+=$
      double precision function dRan3(idum)      
      implicit integer(i-n)
      implicit double precision(a-h,o-z)      
      Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)
      Dimension MA(55)
      save      
      if (idum.lt.0.or.iff.eq.0) then
         iff=1
         mj=mseed-iabs(idum)
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if (mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if (ma(i).lt.mz) ma(i)=ma(i)+mbig
 12         continue
 13      continue
         inext=0        
         inextp=31
         idum=1
      end if
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if (mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      dRan3=mj*fac
      end   
      

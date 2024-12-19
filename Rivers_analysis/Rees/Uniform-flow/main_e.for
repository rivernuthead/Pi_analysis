*######################################################################
* 	  calcolo delle caratteristiche di una sezione secondo Engelund
	  program carsez
	  use comuni
	  use files
	  implicit none
	  integer i,numpunti,j,quantifile,k,numport
	  real*8 port_iniz,profmax
	  character*16 nomesh

	  filectrl='ctrleng.txt'
	  filedati='dati.dat'

c	cost_ks=3.
c	range=0


	  CALL DATE_AND_TIME(DATE=c_data,TIME=c_tempo)

	  write(*,*)'Calcolo delle caratteristiche delle sezioni'
	  write(*,*)'File di controllo:',filectrl
	  open(unit=11,file=filectrl,status='unknown')
	  write(11,*)'data: ',c_data,' ',c_tempo

	  call lettura

      open(unit=14,file=toposez2,status='old')
	  read(14,*) quantifile

	  open(unit=15,file=results,status='unknown')
	  write(15,'(a5, a14,5(a20))') '#', 'slope','discharge','Qs   ',
     @'      meandepth', 'meanwidth', '   activewidth' 
     
c     @, 'active depth''TBI    ','
c     @ABI    ','ragg_idr  ','maxdepth  '


	  do i=1,quantifile


	  	  pend_fin=0.
	  read(14,*) toposez,pend,portfis,numsez,numpunti,deltax,
     @ deltay,partenza
c	write(*,*) toposez,pend,portfis,numsez,numpunti,deltax,deltay,partenza
	
		allocate(numpun(numsez))
c		allocate(numpun(numsez),pp(numsez))
		numpun=numpunti
		maxpun=numpunti
		call lettsez
c		write(*,*) 'letta sezione'
c		write(*,*) i
c		pause

		select case(modo)
			case(6)
		    pend=(pend-pend_fin)/1.1				! per NZ lab
c			if(pend_fin.ne.0.) pend=-pend_fin/1.1	! per lab tn

			case(7)
		    pend=(pend-pend_fin)/1.1

			case(8)
			pend=(pend)/1.1

			case(3)
			pend=(pend)/1.1

		end select
		
!		write(16,*) -9999
c		write(*,*) pend, portfis
c	pause

		port_iniz=portfis
c		call cls
		allocate(hu(numsez),hsh(numsez))
		allocate(jsl(numsez),tetasez(numsez),fimedio(numsez),ba(numsez)
     @		,batt(numsez))

		open(2,file=portate,status='old')
		read(2,*) numport

		do k=1,numport
			read(2,*) portfis
c	write(*,*) 'portata', k
		profmax=0.
		sotto=0.00001
c	write(*,*) portfis, pend
			do j=partenza,numsez
				sez=j

c	write(*,*) j,pend
c	pause
				call eng
				profmax=profmax+(hu(j)-talweg(j))
c	write(*,*) pend, hu(j), talweg(j)
c	pause
			enddo

			profmax=profmax/(numsez+1.-partenza)
c	write(*,*) 'profmax = ', profmax

c	write(*,*) portata, hu
c	pause

			hsh=hu
	
c	write(*,*) 'finito engelund entro in shields'
c	write(*,*) 'sezione', i

			call shields
            fim=fim*sqrt(9.81*delta*diam**3)
			write(15,'(i3, 11(f12.7,2x))') i,pend,portfis,fim,ymedio,
     @bmedio,battmedio2 !,yattmedio,ib,iba2,ragidr,profmax

c		write(*,'(4(f12.7,2x))') pend,portfis,fim,bmedio
			deallocate(largh_sd,ysin,ydes)

		enddo

c		write(15,*)
*		pause
		deallocate(hu,hsh,jsl,tetasez,fimedio,ba,batt)
		deallocate(hh,ll,ks,dx,numpun)
		deallocate(talweg,argini,argin,xx)
		deallocate(numero)
		close(2)
!	close(16)

	  enddo


	  write(11,*)
	  write(11,*)'Fine'
	  close(11)
	  close(15)
c	close(20)

	  stop
	  end
*######################################################################



*######################################################################
	  subroutine lettura
	  use comuni
	  use files
	  implicit none


*	file di dati
  	  open(unit=1,file=filedati,status='old')


*	file di dati sul tratto considerato
      read(1,*)
      read(1,*)toposez2,modo
      read(1,*)
      read(1,*)portate
	  read(1,*)
	  read(1,*)results
      read(1,*)
      read(1,*)cost_ks
 	  read(1,*)
	  read(1,*)delta,diam,formulatrasp
	  read(1,*)
	  read(1,*)sogliateta
	  read(1,*)
	  read(1,*)miny

	  close(1)

*	lettura file sezioni
	  write(11,*)'Dati topografia sezioni in: ',toposez2
	  write(*,*)'Dati topografia sezioni in: ',toposez2


	  return
	  end
*######################################################################



*######################################################################
	  subroutine lettsez
	  use comuni
	  use files
c	  use imsl
	  implicit none
	  real*8:: or
	  integer:: i,j,jj,nzone,tipo(7),oi,oii,n,r
	  real*8:: arg(4),cambi(7),coeff(2),somme(2),statp(10)
	  character*20 seztxt

*	devo aumentare maxpun perché nella lettura dell'Adige sfora in qualche punto
	  maxpun=maxpun+2

	  allocate(hh(maxpun,numsez),ll(maxpun,numsez),ks(maxpun,numsez)
     @,dx(numsez))
	  allocate(talweg(numsez),argini(numsez,2),xx(numsez),argin(numsez))
	  allocate(numero(numsez))
	  numero=-9

*	lettura file sezioni
	  open(unit=1,file=toposez,status='old')
	  write(11,*)
	  write(11,*)'Lettura topografia sezioni modo:',modo
	  select case(modo)
	  case(1)
		write(11,*)'Modo 1, tipo Po di Goro'
		do i=1,numsez
			read(1,*)
			read(1,*)
			read(1,*) dx(i),or,or,nzone
			read(1,*)
			read(1,*) (tipo(j),j=1,nzone)
			read(1,*) (cambi(j),j=1,nzone-1)
			cambi(nzone)=9999.d9
			jj=1
			numpun(i)=0
			do j=1,maxpun+1
				read(1,*) ll(j,i),hh(j,i)
				if (ll(j,i).eq.-10000.) goto 11
				if (tipo(jj).eq.2) ks(j,i)=ksg       !golena
				if (tipo(jj).eq.1) ks(j,i)=ksa       !alveo
				if (ll(j,i).ge.cambi(jj)) jj=jj+1
				numpun(i)=numpun(i)+1
			enddo
			print*,i,j
			stop 'ciclo punti'
11			continue
		enddo
*		quote iniziali del fondo e lunghezza totale
		lung=0.
		xx(1)=dx(1)
		do i=2,numsez
			lung=lung+dx(i)
			xx(i)=lung
		enddo

	  case(2)
		write(11,*)'Modo 2, tipo Adige originario'
		oii=-999999
		n=0
		read(1,*)
12		continue
		read(1,*,end=22)oi,(arg(jj),jj=1,3)
		if(oi.ne.oii) then
			oii=oi
			if(n.ne.0) numpun(n)=j
			n=n+1
			j=1
			xx(n)=arg(1)*1000.
		else
			j=j+1
		endif
		ll(j,n)=arg(2)
		hh(j,n)=arg(3)
		goto 12
22		continue
		numpun(numsez)=j
		or=xx(1)
		xx=xx-or
		lung=xx(numsez)-xx(1)
		dx(1)=0.
		do i=2,numsez
			dx(i)=xx(i)-xx(i-1)
		enddo
		ks=ksa	!@@@@@@@@@@@@@@@@@@@@@@@@@

	  case(3)
		write(11,*)'Modo 3, uscita magicsez'
		read(1,*)
		do n=1,numsez
			read(1,*)i,numpun(n),xx(n)
			do j=1,numpun(n)
				read(1,*)ll(j,n),hh(j,n),ks(j,n)
			enddo
		enddo
		lung=xx(numsez)-xx(1)
		dx(1)=0.
		do i=2,numsez
			dx(i)=xx(i)-xx(i-1)
		enddo

	  case(5)
c	A_tn da database
		read(1,*)
		do n=1,numsez
			read(1,*)numero(n),xx(n)
			xx(n)=xx(n)*1000.
		enddo

		close(1)
		write(11,*)'Lettura topografia da: A_tn_dat.dat'
		!punti
		open(unit=1,file='A_tn_dat.dat',status='old')
		read(1,*)
		n=1
		j=1
		read(1,*,end=925)oii,ll(j,n),hh(j,n)
915		continue
		read(1,*,end=925)oi,(arg(jj),jj=1,2)
		if(oi.ne.oii) then
			numpun(n)=j
			oii=oi
			n=n+1
			j=1
			ll(j,n)=arg(1)
			hh(j,n)=arg(2)
		else
			j=j+1
			ll(j,n)=arg(1)
			hh(j,n)=arg(2)
		endif
		goto 915
925		continue
		numpun(n)=j
		if(n.ne.numsez) stop 'errore lettura 5 (numero sezioni)'


		lung=xx(numsez)-xx(1)
		dx(1)=0.
		do i=2,numsez
			dx(i)=xx(i)-xx(i-1)
		enddo
		ks=ksa

c@@@@@@
	  write(*,*)
	  write(*,*)'Controllo sezioni con i file txt [1->sì]'
	  read(*,*)oi
	  if(oi.eq.1)then
		write(11,*)'Controllo sezioni con i file txt'
		n=0
		write(*,*)'Correggi: 1->sì, 2->sì tutto'
		do i=1,numsez
			!seztxt=numero(i)
			write(seztxt,'(i8,a)')numero(i),'.txt'
			oi=8-int(log10(1.*numero(i)))
			seztxt=seztxt(oi:12)
			write(seztxt,'(a,a)')'A_txt\',trim(seztxt)
			write(11,*)seztxt,i,numero(i)
			open(unit=8,file=trim(seztxt),status='old')
			do j=1,numpun(i)
				read(8,*)oi,oii,(arg(jj),jj=1,4)
	  if(oi.ne.numero(i)) then
		write(11,*)'Errore numero',numero(i),oi
		stop 'Controllo numero'
	  endif
	  if(oii.ne.j) then
		write(11,*)'Errore candela',j,oii
		stop 'Controllo candela'
	  endif
	  if(arg(1).ne.ll(j,i)) then
		write(11,'(a,i5,2e12.3)')'Errore ll',j,ll(j,i),arg(1)
	  write(*,'(a,3i5,2f12.3)')'Errore ll',numero(i),i,j,ll(j,i),arg(1)
		if(n.ne.2) then
			write(*,*)'Correggi?'
			read(*,*)n
			if(n.eq.1.or.n.eq.2)then
				ll(j,i)=arg(1)
			endif
		else
			ll(j,i)=arg(1)
		endif
		write(11,*)'ll nuovo:',ll(j,i)
	  endif
	  if(arg(4).ne.hh(j,i)) then
		write(11,'(a,i5,2e12.3)')'Errore hh',j,hh(j,i),arg(4)
	   write(*,'(a,3i5,2f12.3)')'Errore hh',numero(i),i,j,hh(j,i),arg(4)
		if(n.ne.2) then
			write(*,*)'Correggi?'
			read(*,*)n
			if(n.eq.1.or.n.eq.2)then
				hh(j,i)=arg(4)
			endif
		else
			hh(j,i)=arg(4)
		endif
		write(11,*)'hh nuovo:',hh(j,i)
	  endif
			enddo
87			continue
			read(8,*,end=88)oi,oii,(arg(jj),jj=1,4)
			write(11,'(a,i5,2f12.3)')'oltre',oii,arg(1),arg(4)
			if(oii.eq.numpun(i)+1) write(*,*)'Oltre sez',i,numero(i)
	  if(n.ne.2) then
	  write(*,'(a,3i5,2f12.3)')'Errore oltre',numero(i),i,oii,arg(1),arg(4)
		write(*,*)'Correggi?'
		read(*,*)n
		if(n.eq.1.or.n.eq.2)then
			ll(oii,i)=arg(1)
			hh(oii,i)=arg(4)
		endif
	  else
		ll(oii,i)=arg(1)
		hh(oii,i)=arg(4)
	  endif
			goto 87
88			continue
			if(oii.ne.numpun(i))then
	  write(11,*)'Correzione numero punti',numpun(i),oii
				if(oii.le.maxpun) then
					numpun(i)=oii
				else
					stop 'corr > maxpun'
				endif
			endif
			close(8)
		enddo
	  endif
c@@@@@@

C$$$$$$       case(6)
C$$$$$$ c   file fondo canaletta
C$$$$$$         read(1,*)
C$$$$$$         read(1,*)
C$$$$$$         dx=deltax
C$$$$$$         do i=1,numsez
C$$$$$$             xx(i)=deltax*(i-1)
C$$$$$$             read(1,*) (hh(n,i),n=2,numpun(i)+1)
C$$$$$$             hh(1,i)=hh(2,i)+30.
C$$$$$$             hh(numpun(i)+2,i)=hh(numpun(i)+1,i)+30.
C$$$$$$             do n=1,numpun(i)+2
C$$$$$$                 ll(n,i)=deltay*(n-1)
C$$$$$$                 hh(n,i)=hh(n,i)/1000. !.-xx(i)*ifm
C$$$$$$                 ks(n,i)=ksa
C$$$$$$             enddo
C$$$$$$         enddo
C$$$$$$         lung=xx(numsez)-xx(1)
C$$$$$$ 
C$$$$$$       r=5
C$$$$$$ 
C$$$$$$ c   do i=1,numsez
C$$$$$$ c       do n=2,int((numpun(i)+2)/r)-1
C$$$$$$ c           ll(n,i)=ll(r*n-r+1,i)
C$$$$$$ c           hh(n,i)=hh(r*n-r+1,i)
C$$$$$$ c       enddo
C$$$$$$ c       ll(int((numpun(i)+2)/r),i)=ll(int((numpun(i)+2)/r-1),i)+.05
C$$$$$$ c       hh(int((numpun(i)+2)/r),i)=hh(int((numpun(i)+2)/r-1),i)+0.3
C$$$$$$ c
C$$$$$$ c       do n=int((numpun(i)+2)/r)+1,numpun(i)+2
C$$$$$$ c           ll(n,i)=ll(n-1,i)+.05
C$$$$$$ c           hh(n,i)=hh(n-1,i)
C$$$$$$ c       enddo
C$$$$$$ c
C$$$$$$ c   enddo   
C$$$$$$ 
C$$$$$$ 
C$$$$$$ 
C$$$$$$ 
C$$$$$$       do i=1,numsez
C$$$$$$         argin(i)=0.
C$$$$$$         do j=50,numpun(i)-50
C$$$$$$             argin(i)=argin(i)+hh(j,i)
C$$$$$$         enddo
C$$$$$$         argin(i)=argin(i)/(numpun(i)-99.)
C$$$$$$       enddo
C$$$$$$ 
C$$$$$$ 
C$$$$$$ c   do i=1,numsez
C$$$$$$ c       argin(i)=0.
C$$$$$$ c       do j=int(20/r),int((numpun(i)+2)/r)-int(20/r)
C$$$$$$ c           argin(i)=argin(i)+hh(j,i)
C$$$$$$ c       enddo
C$$$$$$ c       argin(i)=argin(i)/(int((numpun(i)+2)/r)-int(40/r-1))
C$$$$$$ c   enddo
C$$$$$$ 
C$$$$$$ 
C$$$$$$ 
C$$$$$$       CALL drcurv(numsez-20,xx,argin,1,coeff,somme,statp)
C$$$$$$       pend_fin=coeff(2)
c		WRITE(*,*) pend_fin



C$$$$$$       case(7)
C$$$$$$ c   file fondo canaletta Ashmore
C$$$$$$         read(1,*)
C$$$$$$         read(1,*)
C$$$$$$         dx=deltax
C$$$$$$         do i=1,numsez
C$$$$$$ c           xx(i)=deltax*(i-1)
C$$$$$$             read(1,*) xx(i), (hh(n,i),n=2,numpun(i)+1)
C$$$$$$             hh(1,i)=hh(2,i)-100.
C$$$$$$             hh(numpun(i)+2,i)=hh(numpun(i)+1,i)-100.
C$$$$$$             do n=1,numpun(i)+2
C$$$$$$                 ll(n,i)=deltay*(n-1)
C$$$$$$                 hh(n,i)=-hh(n,i)/1000. !.-xx(i)*ifm
C$$$$$$                 ks(n,i)=ksa
C$$$$$$             enddo
C$$$$$$         enddo
C$$$$$$         lung=xx(numsez)-xx(1)
C$$$$$$ 
C$$$$$$     
C$$$$$$       do i=1,numsez
C$$$$$$         argin(i)=0.
C$$$$$$         do j=2,numpun(i)-2
C$$$$$$             argin(i)=argin(i)+hh(j,i)
C$$$$$$         enddo
C$$$$$$         argin(i)=argin(i)/(numpun(i)-3.)
C$$$$$$       enddo
C$$$$$$ 
C$$$$$$       CALL drcurv(numsez,xx,argin,1,coeff,somme,statp)
C$$$$$$       pend_fin=coeff(2)


	  case(8)
c	Sunwapta 2003
		read(1,*)
		read(1,*)
		dx=deltax
		do i=1,numsez
c			xx(i)=deltax*(i-1)
			read(1,*) xx(i), (hh(n,i),n=1,numpun(i))
			do n=1,numpun(i)
				ll(n,i)=deltay*(n-1)
				ks(n,i)=ksa
			enddo
		enddo
		lung=xx(numsez)-xx(1)

	
	  do i=1,numsez
		argin(i)=0.
		do j=2,numpun(i)-2
			argin(i)=argin(i)+hh(j,i)
		enddo
		argin(i)=argin(i)/(numpun(i)-3.)
	  enddo

	  case default
		stop 'modo di lettura non previsto'
	  end select
	  close(1)
*	fine lettura

	  do i=1,numsez
		talweg(i)=hh(1,i)
		argini(i,2)=hh(1,i)
		do j=2,numpun(i)
			if (hh(j,i).lt.talweg(i)) talweg(i)=hh(j,i)
			if (hh(j,i).gt.argini(i,2)) argini(i,2)=hh(j,i)
		enddo
		argini(i,1)=min(hh(1,i),hh(numpun(i),i))
		if(modo.eq.6) argini(i,1)=min(hh(1,i),hh(numpun(i)+2,i))
		if(modo.eq.7) argini(i,1)=min(hh(1,i),hh(numpun(i)+2,i))
	  enddo

      write(11,*)'Fase di lettura sezioni conclusa'
	

*	nome file riassunto lettura sezioni
	  riassez=toposez
	  oi=len_trim(riassez)-3
	  write(riassez(oi:),'(a,a)')'.sez'

c	write(*,*)'Riassunto dopo la lettura sezioni in: ',riassez
	  write(11,*)'Riassunto dopo la lettura sezioni in: ',riassez

c	open(8,file=riassez,status='unknown')
c	write(8,*)'data: ',c_data,' ',c_tempo

c      write(8,'(a30,a30)') 'Numero di sezioni','Numero massimo di punti'
c      write(8,'(i30,i30)') numsez,maxpun
c	write(8,'(a30)') 'Lunghezza totale'
c	write(8,'(f30.3)') lung
c     write(8,'(2a30)') 'Pendenza degli argini','Pendenza del fondo'
c	write(8,*) (argini(1,1)-argini(numsez,1))/lung
c     @,(talweg(1)-talweg(numsez))/lung
c      write(8,'(a)')'ks (a,g)'
c	write(8,*) ksa,ksg
c      write(8,*)
c	write(8,'(a4,5a15,a8,a10)')'sez','x','talweg','argine1','argine2','dx'
c     @,'npunti','numero'
c	do i=1,numsez
c		write(8,'(i4,5f15.3,i8,i10)')i,xx(i),talweg(i),argini(i,1),argini(i,2)
c     @	,dx(i),numpun(i),numero(i)
c	enddo
c	close(8)

	  return
	  end
*######################################################################


*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
	  subroutine eng
	  use comuni
	  use files
	  implicit none
	  character*12 fsez
	  real*8 quota,fr,hcrit,qtrad
  	  integer i,n


*	calcolo e stampa caratteristiche
c	write(*,*)'Inizio calcolo'

c	open(unit=9,file='engris.txt',status='unknown')
*	write(9,*)'data: ',c_data,' ',c_tempo

*	allocazione
	  allocate(area(1),beta(1),cqr(1),larg(1))

cc	  n=(argini(sez,1)-talweg(sez))/dz

*	calcolo portata con il metodo della bisezione

		quota=argini(sez,1)-sotto
		call engtutto(sez,quota,1)
		qalta=portata
		quotaalta=quota
		quota=talweg(sez)+sotto
		call engtutto(sez,quota,1)

		qbassa=portata
		quotabassa=quota
123	  continue
		quota=(quotaalta+quotabassa)*0.5
		call engtutto(sez,quota,1)
c		write(*,*) quota
	  if(abs((portata-portfis)/portfis).lt.0.001) goto 124
		if(portata.lt.portfis) then
			qbassa=portata
			quotabassa=quota
			goto 123
		else
			qalta=portata
			quotaalta=quota
			goto 123
		endif

124		continue

	  hu(sez)=quota

	  deallocate(area,beta,cqr,larg)

	  return
	  end

*###################################################################################################


*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
*     ENGELUND calcola C^2*Rh, AREA, ALFA, BETA, LARG      *
*     HORTON calcola Qhort, KSeq, cba, cbg                 *
*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
*     cqr = C^2 * R idraulico, cb = contorno bagnato
*     cba = cont.bagn. alveo, cbg = cont.bagn. golena
      subroutine engtutto(n,quota,nn)
      use comuni
	  implicit none
	  integer i,n,nn
	  real*8 dl,ya,yb,quota,dhh,cba,cbg,rh,ym

      area(nn)=0.
      larg(nn)=0.
      alfa=0.
      beta(nn)=0.
      cqr(nn)=0.
	  cb=0.
	  cba=0.
	  cbg=0.
      yb=max(quota-hh(1,n),0.d0)
	  do i=2,numpun(n)
	  ya=yb
	  yb=max(quota-hh(i,n),0.d0)
        if((ya+yb).gt.0.d0) then
		dhh=abs(hh(i,n)-hh(i-1,n))
		if(dhh.lt.0.0001) then
			dl=ll(i,n)-ll(i-1,n)
		else
			dl=abs(ya-yb)*(ll(i,n)-ll(i-1,n))/dhh
		endif
		ym=(ya+yb)/2.
          rh=ym/sqrt(1.+((ya-yb)/dl)**2.)
		larg(nn)=larg(nn)+dl
          area(nn)=area(nn)+ym*dl
		if(cost_ks.gt.0.1) then
			ks(i,n)=6.+2.5*log(rh/(cost_ks*diam))			!ks logaritmico
			ks(i,n)=ks(i,n)*(9.81**.5)/(rh**(1./6.))
		endif

		cqr(nn)=cqr(nn)+ks(i,n)*rh**(2./3.)*ym*dl
          alfa=alfa+ks(i,n)**3.*rh**2.*ym*dl
          beta(nn)=beta(nn)+ks(i,n)**2.*rh**(4./3.)*ym*dl
		cb=cb+(dl**2.+(ya-yb)**2.)**.5
C$$$$$$         if(ks(i,n).eq.ksa) then
C$$$$$$             cba=cba+(dl**2.+(ya-yb)**2.)**.5
C$$$$$$         else
C$$$$$$             cbg=cbg+(dl**2.+(ya-yb)**2.)**.5
C$$$$$$         endif
	  endif

      enddo
      portata=cqr(nn)*sqrt(pend)
	  beta(nn)=area(nn)*beta(nn)/cqr(nn)**2.
      alfa=area(nn)**2.*alfa/cqr(nn)**3.
      cqr(nn)=(cqr(nn)/area(nn))**2./gconst
ccc	  kseq=((cba+cbg)/(cba/ksa**1.5+cbg/ksg**1.5))**(2./3.)
ccc      qhort=area(nn)*kseq*sqrt(pend)*(area(nn)/(cba+cbg))**(2./3.)
      return
      end


*#####################################################################



	  subroutine shields

	  use comuni
	  use files
  	  implicit none

	  integer(2) color
	  integer i,n,sottoacqua,w
	  real*8 yunif,yunif2,ksln,ff,dl,ym,rh,cq,klog,www,fi_C
	  real*8 inizio(1000),fine(1000),soglialargh,q_can(100),sogliaport
	  real*8 ymediosez,ya,yb,dhh,omega,cbagn
	  real*8,allocatable:: fi_can(:),batt_c(:)
	  integer,allocatable:: canali(:),c_att(:),can_att(:)
	  logical acqua,sabbia

	  allocate(largh_sd(numsez),ysin(numsez),ydes(numsez),canali(numsez)
     @	,c_att(numsez),can_att(numsez),fi_can(100),batt_c(numsez))


c	do w=1,10
c		hsh=hsh+.003


	  soglialargh=80.*diam
	  sogliaport=portfis*0.005
c	sogliateta=0.047
	  tetac=0.047
	  ba=0.d0
	  batt=0.d0
	  batt_c=0.
	  fim=0.d0
	  ib=0
	  iba=0
	  iba2=0
	  bmedio=0.
	  battmedio=0.
	  battmedio2=0.
	  yattmedio=0.
	  fimedio=0.d0
	  ymedio=0.
	  sottoacqua=0
	  ragidr=0.
31	  continue
	  do i=partenza,numsez

c		if(modo.eq.3) pend=pp(i)/1.1
	
		canali(i)=0
		c_att(i)=0
		can_att(i)=0
		acqua=.false.
		sabbia=.false.
		cq=0.
		fi_c=0.
		sottoacqua=0
		omega=0.
		cbagn=0.
		ymediosez=0.
		yattmediosez=0.
c	write(*,*) '1',i
!		jsl(i)=pend !max((hsh(i-1)-hsh(i+1))/(xx(i+1)-xx(i-1)),1.e-5)
		ysin(i)=ll(numpun(i),i)
		ydes(i)=ll(1,i)
	    yb=max(hsh(i)-hh(1,i),0.d0)


		do n=2,numpun(i)
c	write(*,*) n,ya,yb,rh,dl,teta
c	pause
			ya=yb
			yb=max(hsh(i)-hh(n,i),0.d0)

			ym=(ya+yb)*.5
			if(ym.gt.(miny*diam)) then
				sottoacqua=sottoacqua+1
				ysin(i)=min(ysin(i),ll(n,i))
				ydes(i)=max(ydes(i),ll(n,i))

				dhh=abs(hh(n,i)-hh(n-1,i))
					if(dhh.lt.0.0001) then
				dl=ll(n,i)-ll(n-1,i)
				else
					dl=abs(ya-yb)*(ll(n,i)-ll(n-1,i))/dhh
				endif
				omega=omega+ym*dl
				cbagn=cbagn+sqrt(dl**2+dhh**2)

c				dl=ll(n,i)-ll(n-1,i)
				ymediosez=ymediosez+ym
				rh=ym*dl/sqrt(dl**2+(ya-yb)**2)
				teta=pend*rh/delta/diam
c	write(16,*) teta
				fi=0.d0
				call calcfi(n,i,ym,rh)
				fimedio(i)=fimedio(i)+fi*dl
				ba(i)=ba(i)+dl
				batt_c(i)=batt_c(i)+dl
				if(teta.gt.sogliateta) batt(i)=batt(i)+dl
				if(teta.gt.sogliateta) yattmediosez=yattmediosez+ym*dl

				if(acqua) then
					fine(canali(i))=ll(n,i)

					if(.not.sabbia) then
						if(teta.gt.sogliateta) then
							sabbia=.true.
							c_att(i)=c_att(i)+1
						endif
					endif
					if(cost_ks.eq.0) then
						klog=ks(n,i)
					else
	 					klog=6.+2.5*log(rh/(cost_ks*diam))	!ks logaritmico
						klog=klog*(9.81**.5)/(rh**(1./6.))
					endif
					cq=cq+klog*rh**(2./3.)*ym*dl
					fi_c=fi_c+fi*dl

				else 
					canali(i)=canali(i)+1
					can_att(i)=can_att(i)+1
					inizio(canali(i))=ll(n,i)
					acqua=.true.

					if(.not.sabbia) then
						if(teta.gt.sogliateta) then
							sabbia=.true.
							c_att(i)=c_att(i)+1
						endif
					endif

 					if(cost_ks.eq.0) then
						klog=ks(n,i)
					else
	 					klog=6.+2.5*log(rh/(cost_ks*diam))		!ks logaritmico
						klog=klog*(9.81**.5)/(rh**(1./6.))
					endif
					cq=cq+klog*rh**(2./3.)*ym*dl
					fi_c=fi_c+fi*dl

				endif
			
			else
				if(acqua) then
					fine(canali(i))=ll(n,i)

				    q_can(canali(i))=cq*sqrt(pend)
	  fi_can(canali(i))=fi_c/(fine(canali(i))-inizio(canali(i)))
					cq=0.
					fi_c=0.
*	write(*,*) canali(i), q_can(canali(i))
*	pause

					www=fi_can(canali(i))
					call calcteta(www) 
					if(tetamedio.lt.sogliateta) then
						can_att(i)=can_att(i)-1
		batt_c(i)=batt_c(i)-(fine(canali(i))-inizio(canali(i)))
					endif
		if((fine(canali(i))-inizio(canali(i))).lt.soglialargh) then
						canali(i)=canali(i)-1
						if (sabbia) c_att(i)=c_att(i)-1
					elseif(q_can(canali(i)).lt.sogliaport) then
						canali(i)=canali(i)-1
						if (sabbia) c_att(i)=c_att(i)-1
					endif
					acqua=.false.
					sabbia=.false.
				endif
			endif
			
		enddo
		ragidr=ragidr+omega/cbagn
c	write(*,*) ba(i)
		ymediosez=ymediosez/(sottoacqua*1.d0)
		largh_sd(i)=ydes(i)-ysin(i)
		fimedio(i)=max(fimedio(i)/ba(i),0.)
		tetasez(i)=(fimedio(i)/8.)**(2./3.)+tetac
		ff=fimedio(i)
		ib=ib+canali(i)
		iba=iba+can_att(i)
		iba2=iba2+c_att(i)
		fim=fim+fimedio(i)*ba(i)
		bmedio=bmedio+ba(i)
		battmedio=battmedio+batt_c(i)
		battmedio2=battmedio2+batt(i)
c		yattmedio=yattmedio+yattmediosez/batt(i)
		ymedio=ymedio+ymediosez
	  enddo

	  bmedio=bmedio/(numsez+1.-partenza)
	  fim=fim/(numsez+1.-partenza)
	  ymedio=ymedio/(numsez+1.-partenza)
	  ib=ib/(numsez+1.-partenza)
	  iba=iba/(numsez+1.-partenza)
	  iba2=iba2/(numsez+1.-partenza)
	  battmedio=battmedio/(numsez+1.-partenza)
	  battmedio2=battmedio2/(numsez+1.-partenza)
c	  yattmedio=yattmedio/(numsez+1.-partenza)
c	  ragidr=ragidr/(numsez+1.-partenza)



10	  format(i5,9f12.5,3i5)

c	deallocate(canali,c_att)

c	enddo

 	  return

	  end


******************************************************************************

	  subroutine calcfi(n,i,ym,rh)

	  use comuni
	  use files
	  implicit none

	  integer n,i
	  real*8 fi0,csi,ym,kslog,rh,q_str,om_p,om_c,Y_c,c_bag

	  if (formulatrasp.eq.1) then !MPM corrected by Wong
		if(teta.gt.tetac) then
			fi=4.93*(teta-tetac)**1.6
		else
			fi=0.
		endif
		return
	  endif
	
	  if (formulatrasp.eq.2) then !Parker 1990
		fi0=0.00218*teta**1.5
		csi=teta/0.0386
		fi=fi0*EXP(14.2*(csi-1)-9.28*(csi-1.)**2)
		if(csi.gt.1.5878) fi=fi0*5474.*(1.-0.853/csi)**4.5
		if(csi.lt.1.) fi=fi0*csi**14.2
		return
	  endif

	  if (formulatrasp.eq.3) then	!Bagnold 1980
		if(cost_ks.eq.0) then
			kslog=ks(n,i)
		else
	 		kslog=6.+2.5*log(rh/(cost_ks*diam))		!ks logaritmico
			kslog=kslog*(9.81**.5)/(rh**(1./6.))
		endif
		q_str=kslog*rh**(2./3.)*ym*sqrt(pend)
		
		om_p=1000*pend*q_str

		Y_c=0.04*1.65*diam/pend
		c_bag=5.75*LOG(12*Y_c/diam)/log(10.)
		om_c=((0.04*1.65)**1.5)*sqrt(9.81)*(diam**1.5)*1000*c_bag
		fi=0.
		if(om_p.gt.om_c) then
			fi=0.1*(((om_p-om_c)/(0.5))**1.5)*((ym/.1)**(-.6))*((
     &			diam/.0011)**(-.5))/sqrt(9.81*1.65*diam**3)/1650.
		endif
		return
	  endif





	  write(*,*) 'formula di trasporto errata'

	  stop

	  end


******************************************************************************

	  subroutine calcteta(ffi)

	  use comuni
	  use files
	  implicit none
	  real*8 ffi,rh,ym
	  integer n,i
	
	  rh=1.
	  ym=1.
	  n=1
	  i=numsez
	  teta=0.0001d0

5	  continue
	  call calcfi(n,i,ym,rh)

	  if(fi.lt.ffi) then
		teta=teta+0.0001d0
		goto 5
	  endif

	  tetamedio=teta
	  return
	  end





	
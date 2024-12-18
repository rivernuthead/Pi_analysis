*######################################################################
	  module comuni
	  implicit none

	  real*8,parameter:: gconst=9.81d0,ro=1000

*	dati sezioni
	  integer:: numsez,maxpun,modo,partenza,nclassi,range
	  integer,allocatable:: numpun(:)
	  real*8,allocatable:: hh(:,:),ll(:,:),ks(:,:),dx(:),argin(:)
	  real*8,allocatable:: xx(:),talweg(:),argini(:,:)
	  real*8 ifm,ifv,ksa,ksg,lung,cost_ks,deltax,deltay

*	caratteristiche sezioni
	  real*8 portata,pend,dz,pend_fin

	  real*8,allocatable::area(:),beta(:),cqr(:),larg(:),hu(:) !u(:),bu(:)

	  real*8 alfa,cb,port,portfis,sotto,ragidr
	  real*8 kseq,qhort,qalta,qbassa,quotaalta,quotabassa

	  integer sez

*	numero cicli, tolleranza
	  integer nccrit
	  real*8 tolcr

*	data e ora
	  character:: c_data*8,c_tempo*10
*	shileds
	  real*8,allocatable:: jsl(:)
	  real*8,allocatable:: fimedio(:),ba(:),hsh(:),tetasez(:),batt(:)
	  real*8,allocatable:: largh_sd(:),ysin(:),ydes(:),pp(:)
	  real*8 tetamedio,delta,diam,fim,tetac,bmedio,tetaunif,ib,teta,fi
	  real*8 iba2,battmedio,battmedio2,ymedio,yattmedio,yattmediosez
	  real*8 sogliateta,miny,iba
 	  integer formulatrasp



	  integer,allocatable:: numero(:)

	  save
	  end
*######################################################################

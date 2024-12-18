*######################################################################
	  module files
 	  implicit none

*	INPUT
	  character:: filectrl*12 !file di controllo
	  character:: filedati*120 !file dati iniziali

	  character:: datisez*24 !file con dati (if,ks...) e nomi altri files
	  character:: toposez*120 !file con la topografia delle sezioni
	  character:: toposez2*120 !file con la topografia delle sezioni per lab
	  character:: portate*120 !file con le portate
	  character:: results*120 !file con i risultati


*	OUTPUT
	  character:: riassez*120 !riassunto sezioni dopo la lettura
	  character:: filebeta*12 !beta interpolato
	  character:: filearea*12 !area interpolata
	  character:: filelarg*12 !larghezza interpolata
	  character:: fileresi*12 !resistenza interpolata

	  save
	  end
*######################################################################

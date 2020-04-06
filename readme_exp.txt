VERSIONE 8_03_2020

Funzionamento dettagliato del programma Experimentum_united.py
Librerie di Appoggio: lib_experimentum, Models, Alessandria

Fase 0) PREPARAZIONE AL LANCIO

Si utilizza un altro programma: pre_Experimentum_united.py:
	- Partendo dai dati di una certa acquisizione ("matrici" di spettri ordine 10^3 10^4)
	  Il programma grazie a find_peaks() di scipy trova tutti i picchi per ogni spettro 
	  con delle caratteristiche poco restrittive su altezza (circa 20 conteggi) mentre su width e distance si usano dei valori
	  che sono risultati essere caratteristici: rispettivamente 2.5 e 30 circa
	  --> VENGONO CLASSIFICATI GLI SPETTRI IN QUATTRO GRUPPI: SATURATED, FOUR, MORE THAN FUOR E LESS_THAN_FOUR
	- un ciclo successivo trova 4 picchi per ogni spettro, riempiendo quindi FOUR e svuotando gli altri
	  eventuali problemi individuano spettri che non vanno (magari è da ritoccare la width, o etc..)
	- sugli spettri individuati si effettua una statistica su altezza, distanza e width dei 4 picchi in gioco:
	  MI ASPETTO CHE LE DISTRIBUZIONI SIANO NORMALI
	  grazie anche a rappresentazione grafica, individuo spettri che non vanno (errore frequente è il take into account del 
    	  picco brillouin del vetro che quindi ovviamente presenta anomalia nelle distanze tra i picchi)
	- una volta che non ci sono più anomalie, il programma stampa i valori di parametri che sono fondamentali per 
 	  il funzinamento del programma Experimentum_united.py:
	  valore minimo dell'altezza brillouin
	  valore minimo altezza degli elastici
	  valore medio distanza tra picchi brillouin-elastici per il taglio
          valori minimo distanza e width dei picchi

	- bisogna trovare negli spettri i primi (in senso cronologico di acquisizione, tendenzialmente a serpentina da in alto a destra)
	  prototipi di normals, almost_height e brillouin_high* 
          questi andranno trattati con Experimentum_single.ipynb e andranno inseriti come input di Experimentum_united.py i risultati dei paramentri dei fit 
          di questi, una volta che sono buoni

Fase 1) ACQUISIZIONE E MODIFICA SPETTRI
	
	- Initialize_Matrix genera una tupla di tuple contenente oggetti Spectrum (def in lib_experimentum) 
	  per ora vuoti, con il nome (ii,jj) indici di riga, colonna che sarà la nomenclatura std per chiamare il singolo 
	  elemento
 	- genero una cartella nella directory dove si trovano i dati e ci salvo un file di log sul quale scriverò varie info
	- Get_Spectrum() prende gli spettri dai dati (solitamente formato matlab)
	- Get_Spectrum_Peaks(**syg_kwargs) prende tutti i picchi con altezza maggiore uguale del Brillouin più basso (vd sopra)
 	  --> sono sicuro che gli spettri con meno di 4 picchi sono del tipo 
	      brillouin_higher o brillouin_highest ossia hanno brillouin più alti di uno o di entrambi gli elastici
	- eliminati gli spettri saturati, faccio un ciclo in cui definisco definitivamente le sotto categorie della macro categ boni
	  che contiene gli spettri sui quali farò il fit
 
		normals = spettri con elastici più alti ma non troppo, i più comuni
		almost_height = spettri con elastici alti oltre una treshold che cominciano a mangiarsi i brillouin
		brillouin_higher = brillouin più alto di un elastico
		brillouin_highest = brillouin più alto di tutti gli elastici 

	- per tutti gli spettri sono quindi selezionati 4 picchi: 
	 	- per chi ne ha di più, prendo i 4 più alti
		- per chi ne ha di meno (i brillouin_high*) faccio un retake Get_Spectrum_Peaks(**syg_kwrgs_brillouin) dove quest
		  kwargs hanno l'altezza più bassa degli elastici registrata da pre_exp() (vd sopra)
		  --> SONO SICURO che sto prendendo solo i 4 boni
	
	- parallelamente agisco sul VIPA(solo su elemento (0,0) perchè per ora ssumiamo tutto uguale)
	  1. faccio FIT per conversione da pixel in GHz con Fit_Pix2GHz() basandomi sulla conoscenza teorica che i picchi elastici 
	  sono distanti il free spectral range
	  2. fitto la gaussiana che inviluppa lo spettro con 3 picchi elastici
	- converto tutti spettri in frequenze, li allineo, nel senso che a seconda dell'inviluppo gaussiano so se il secondo ordine di 
	  elastici deve stare a sx o dx del principale, poi taglio gli spettri prendendo la metà della distanza media tra brillouin ed elastici data dal 	   pre_exp()	 
	
Fase 2) FIT MARKOVIANO
	 
	procedo con il fit modello markoviano ( perchè so che viene bene e potrà poi essere il punto di partenza per quell'altro)
	
	- la cosa fondamentale sono le condizioni iniziali:
	  faccio il fit scorrendo come ordine di acquisizione e prendo come p0 il migliore (in termini di chiquadro) tra i vicini che hanno 
	  effettivamente dei parametri Markov (quindi tendenzialemente non tanti quanti sono i primi vicini, 8, ma solo quelli sopra) e anche il p0 
	  prototipo che ho definito all'inizio. 
 	  A SECONDA DELL'APPARTENENZA A UNA CATEGORIA, AGGIORNO IL P0 PROTOTIPO
 	  i bounds sono scelti in modo che la gaussiana possa muoversi poco

Fase 3) Fit totale




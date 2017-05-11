Nella cartella tomtool/PattRec ci sono tutte le classi necessarie per fare la Pattern Recognition.

Nella cartella tomtool/Monitor si trova la classe utilizzata per il calcolo dei TTrig.


1- Prima di compilare e far girare il programma "runPR" e' necessario settare le variabili 
   d'ambiente corrette, che si trovano in /usr/local/bin/root.sh:

      source /usr/local/bin/root.sh


2- Nella cartella tomtool/PattRec compilare il programma (nel caso si siano fatte modifiche
   al programma):

      gmake clean
      gmake


3- Lanciare il programma di ricostruzione "runPR": 
   selezionando in modo opportuno la flag "ttrig" in RawAnalyzer.C, e' possibile scegliere se
   calcolare i ttrig o se far andare la ricostruzione delle tracce. 
   E' possibile settare la flag "ttrig" da riga di comando, come segue:
   
   Per calcolare i ttrig:
      ./runPR NOME_FILE_SENZA_ESTENSIONE N.RUN N.RUN_TTRIG NUMERO_EVENTI RUNID_MIN 1=CALCOLA_TTRIG  >&  log_trig_Nrun.txt
 
   Per lanciare la ricostruzione delle tracce:   
      ./runPR NOME_FILE_SENZA_ESTENSIONE N.RUN N.RUN_TTRIG NUMERO_EVENTI RUNID_MIN 0=COMPUTE_PR     >&  log_Nrun.txt

   Prima di lanciare la ricostruzione delle tracce per un certo N.RUN e' necessario aver creato 
   il file con i ttrig relativi al run che si vuole ricostruire!
   I file con i ttrig si trovano in tomtool/PattRec/ttrig/.

   Esempio:
      Calcolo dei ttrig per il run 1253:
      ./runPR 18tubi_3300_2900_30mv_r1253 1253 1253 60000 1 1  >&  log_trig_1253.txt
      Ricostruzione run 1253:
      ./runPR 18tubi_3300_2900_30mv_r1253 1253 1253 600000 1 0 >&  log_1253.txt


NB: da riga di comando e' possibile decidere:
       - il numero di run da analizzare (N.RUN);
       - il run da utilizzare per la sottrazione dei ttrig (N.RUN_TTRIG);
       - il numero massimo di eventi da analizzare (NUMERO_EVENTI);
       - il run da cui partire per l'analisi (RUNID_MIN), 
            es: RUNID_MIN=2 significa che il primo run analizzato e' NOMEFILE_Nrun.i2;
       - se calcolare il ttrig o fare la ricostruzione (CALCOLA_TTRIG).

       Il programma analizza 40 file a partire da RUNID_MIN, pero' se il NUMERO_EVENTI impostato
       dall'utente e' minore del numero di eventi presenti in 40 file, il programma si ferma prima.


In /software/tomtool/PattRec/ c'e' uno script che puo' essere modificato ed usato per la ricostruzione:
go_runPR_new.csh



printf "Processing Trigger for test Run 2276 \n"
./runPR /home/sara/data/experiments/LEMMA/raw/MuBlast_phi_carotavuota_riempita_newrif_r2276  2276 2276 100000 0 1 >& log_trig_2276.txt
printf "Processing Tracks for test Run 2276 \n"
./runPR /home/sara/data/experiments/LEMMA/raw/MuBlast_phi_carotavuota_riempita_newrif_r2276  2276 2276 1000000 0 0 >& log_2276.txt
printf "THE END"

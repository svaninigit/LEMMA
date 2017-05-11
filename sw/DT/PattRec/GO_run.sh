#source /usr/local/bin/root32.csh
#printf "Processing Trigger for Run 2150\n"
#./runPR /mudata/LNL/muonegrafia/MuBlast_M204_r2150 2150 2150 1000000 0 1 >& log_trig_2150.txt
#printf "Processing Tracks for Run 2150\n"
#./runPR /mudata/LNL/muonegrafia/MuBlast_M204_r2150 2150 2150 5000000 0 0 >& log_2150.txt
#
#printf "Processing Trigger for Run 2151\n"
#./runPR /mudata/LNL/muonegrafia/MuBlast_M204_phi_r2151 2151 2151 1000000 0 1 >& log_trig_2151.txt
#printf "Processing Tracks for Run 2151\n"
#./runPR /mudata/LNL/muonegrafia/MuBlast_M2161_phi_r2151 2151 2151 5000000 0 0 >& log_2151.txt
printf "Processing Trigger for Run 2161\n"
./runPR /mudata/LNL/muonegrafia/MuBlast_phi_204_r2161 2161 2161 5000000 0 1 >& log_trig_2161.txt
printf "Processing Tracks for Run 2161\n"
./runPR /mudata/LNL/muonegrafia/MuBlast_phi_204_r2161 2161 2161 55000000 0 0 >& log_2161.txt
printf "THE END"

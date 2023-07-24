set echo
gfortran -w -O3 makeinput.f -o makeinput
./makeinput
gcc -w -O3 Interfaccia.c -o Interfaccia
gfortran -w -O3 Monte.f -o Monte
./Montecarlo
gfortran -w -O3 elabora.f -o elabora
./elabora


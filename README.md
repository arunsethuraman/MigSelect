# Here, tell people how to use TODO:

# PROBLEMS: 
1. can't make bc of this: 
```mpicc *.c -D MIGSELECTRELEASE -o MigSelect -lm
autoc.c: In function 'set_autoc_vals':
autoc.c:242:11: error: implicit declaration of function 'whichiscoldchain' [-Wimplicit-function-declaration]
  242 |   int z = whichiscoldchain();
      |           ^~~~~~~~~~~~~~~~
ima_main.c: In function 'scan_commandline':
ima_main.c:567:24: warning: assignment discards 'const' qualifier from pointer target type [-Wdiscarded-qualifiers]
  567 |               for (opt =  (pstr, sep); opt; opt = strtok (NULL, sep))
      |                        ^
ima_main.c: In function 'qupdate':
ima_main.c:1525:7: error: implicit declaration of function 'whichiscoldchain' [-Wimplicit-function-declaration]
 1525 |   z = whichiscoldchain();
      |       ^~~~~~~~~~~~~~~~
In file included from imamp.h:11,
                 from update_t_RY.c:8:
update_t_RY.c: In function 'changet_RY1':
update_t_RY.c:479:21: warning: statement will never be executed [-Wswitch-unreachable]
  479 |       assert (pdgnew[ui] == 0);
      |               ~~~~~~^~~~
make: *** [Makefile:7: all] Error 1
```
to fix, maybe add whichiscoldchain() to a header file and import it in `autoc.c`

1. how to set up env

`conda env create -f environment.yml`

`conda activate migselect_env`
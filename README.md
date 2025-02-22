
# MDC: Classical Molecular Dynamics Simulation Package in C

##  Please refer to the Sec. 7.1 of
      http://web.cc.iwate-u.ac.jp/~nisidate/main.pdf
  for more information.	    

###                   kazume.nishidate@gmail.com


  * The 'examples' directory contains following systems.
  
**    NaCl.hm:  NaCl with Huggins-Mayer potentail (default system)
**    NaCl.sx1: NaCl with sx1 potentail proposed by K. Kawamura
**    AgI:      Super ionic conductor with the soft-core potential
**    CaF2:     Super ionic conductor with the soft-core potential
**    SW.Si:    Silicon with Stillinger Weber potential
**    AT.C60:   C60 with Abell-Tersoff potential
**    AT.Ge:    Germanium with Abell-Tersoff potential
**    AT.C60:   Si with Abell-Tersoff potential    

  * Make the code using the following command
      make
    and run the code with
      ./md
    command.
    
  * To make the X-Window version, use the following commands.
      make clean   
      make -f Makefile.X



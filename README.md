

# MDC: Classical Molecular Dynamics Simulation Package in C

                 kazume.nishidate@gmail.com

### Manual (in Japanese)
     Please refer to the Sec. 7.1 of [the manual page (Japanese)](http://web.cc.iwate-u.ac.jp/~nisidate/main.pdf]or more information.	    
### The 'examples' directory contains following systems.
  
- NaCl.hm:  NaCl with Huggins-Mayer potentail (default system)
- NaCl.sx1: NaCl with sx1 potentail proposed by K. Kawamura
- AgI:      Super ionic conductor with the soft-core potential
- CaF2:     Super ionic conductor with the soft-core potential
- SW.Si:    Silicon with Stillinger Weber potential
- AT.C60:   C60 with Abell-Tersoff potential
- AT.Ge:    Germanium with Abell-Tersoff potential
- AT.C60:   Si with Abell-Tersoff potential    

### How to compile.
    Make the code using the following command
      make
    and run the code with
      ./md
    
    
### X window.
    To make the X-Window version, use the following commands.
      make clean   
      make -f Makefile.X



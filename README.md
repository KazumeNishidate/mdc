

# MDC: Classical Molecular Dynamics Simulation Package in C
[Kazume_NISHIDATE](mailto:kazume.nishidate@gmail.com)
			  
-------------------------------------------------------------------------------

### Manual (in Japanese)
Please take a look at the Sec. 7.1 of the book
 [C言語によるコンピュータシミュレーション](http://web.cc.iwate-u.ac.jp/~nisidate/main.pdf), for more information.

	 
### The 'examples' directory contains the following systems.
  
NaCl.hm:  NaCl with Huggins-Mayer potential (default system)

- NaCl.sx1: NaCl with sx1 potential proposed by K. Kawamura
- AgI:      Superionic conductor with a soft-core potential
- CaF2:     Superionic conductor with the soft-core potential
- SW.Si:    Silicon with Stillinger Weber potential
- AT.C60:   C60 with Abell-Tersoff potential
- AT.cnt55: (5,5) Carbon nanotube with Abell-Tersoff potential
- AT.Ge:    Germanium with Abell-Tersoff potential
- AT.Si:   Si with Abell-Tersoff potential

### How to compile.

   Make the code using the following command
  
  	make

   and run the code with

  	./md

If you encounter a 'multiple definition error' during the GCC compilation, you can avoid it by adding the option '-z muldefs' in the Makefile.

### X window.

To make the X-Window version, use the following commands.

      make clean   
      make -f Makefile.X


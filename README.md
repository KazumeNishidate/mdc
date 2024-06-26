
  MDC: Classical Molecular Dynamics Simulation Package in C

  kazume.nishidate@gmail.com
-------------------------------------------------------------

  * By default, it calculates a microcanonical ensemble.

  * The example directory contains other model systems.

  * To make the X-Window version, make the code using the following
    command.
    
      make clean   
      make -f Makefile.X

-------------------------------------------------------------
  * References 

  [1] Y. Hiwatari, Solid State Physics (KOTAIBUTSURI: Japanese),
       <Vol. 17, No. 3, 141 (1982)>, <Vol. 17, No. 4, 197 (1982)>,
       <Vol. 17, No. 6, 317 (1982)>, <Vol. 17, No. 8 453 (1982)>,
       <Vol. 24, No. 277, pp108(242)-118(252) (1989)>.

  [2] D. W. Heermann, "Computer Simulation Methods in Theoretical
       Physics" (2nd edn.), Springer-Verlag, (1990).

  [3] K. Kawamura, "PASOKON BUNSHI SIMULATION", (Japanese), KAIBUNDO,
       (1990).

  [4] K. Kawamura, Solid State Physics (KOTAIBUTSURI: Japanese),
       Vol. 29, No. 12, 909 (1994).

  [5] A. Ueda, "COMPUTER SIMULATION", (Japanese), ASAKUEA SYOTEN,
      (1990).

  [6] M. P. Tosi and F. G. Fumi, J. Phys. Chem. Solids 23, 359 (1962).

  [7] M. Mori, "A Method for Evaluation of the Error Function of Real
       and Complex Valiable with High Relative Accuracy", RIMS, Kyoto
       University Japan, 19, 1081 (1983).

  [9] S. Okazaki, "COMPUTER SIMULATION NO KISO" (Japanese),
       Kagaku-dojin, (2000).
-------------------------------------------------------------

<< MD calculational unit >>---------------------------------------------

   [mass]        =>  1.6605402 x 10^{-27}  [Kg]:  atomic mass unit [m_u] 
   [length]      =>  1 x 10^{-10} [m] = 1 [A]                                
   [second']     =>  4.07497263794 x 10^{-15} [s]                          
   [energy]      =>  1 x 10^{-18} [J]                                        
   [temperature] =>  1 [K]                                                 
   [pressure]    =>  1000 [Giga Pascal] = 1 [Tera Pascal]

<< time unit conversion >>---------------------------------------------

  [J] = [N][m] = [kg][m^2]/[s^2]

  [s^2] = [kg][m^2]/[J]

        = {1/{1.6605402 * 10^{-27}}} * 10^{20} / 10^{18} 

                                         [ [mass][length^2]/[energy] ]

        = {1/{1.6605402 * 10^{-27}}} * 10^{20} / 10^{18} 

                                                       [ [second']^2 ]

  [second'] = 

     In[1]:= Sqrt[ 1/{ {1/{1.6605402 * 10^{-27}}} * 10^{20} / 10^{18} }]

                           -15
     Out[1]= {{{{4.07497 10   }}}}
     
     In[2]:= N[%,20]
     
                                     -15
     Out[2]= {{{{4.074972637944947 10   }}}}  [s]

           = 4.074972637944947 [fs]

<< Pressure unit conversion >>------------------------------------------

   [standard atmosphere] = 760 [mmHg] = 101325 [Pa]                      
                         = 1.01325 x 10^{-4} [GPa]                        
                                                                         
   1 [Pa] = [N]/[m^2],  [N] = [kg][m]/[s^2]                              
   1 [Pa] = ( 4.074972637944947^{2}/1.6605402 )                                
                x 10^{-13}   [mass]/([second'^{2}][length])              
                                                                         
   where, 1 [mass]/([second'^{2}][length]) = 1 [Pa'],                     
   and [Giga] = 10^{9}, [Tera] = 10^{12}.                                

                    1.6605402                                               
   1 [Pa'] = ------------------------ x 10^{4}  [GPa]                          
               4.074972637944947^{2}                                           
                    1.6605402                                               
           = ------------------------ x 10000  [GPa]
               4.074972637944947^{2}                                           

           = (1.6605402 * 10000)/(4.074972637944947^{2})  [GPa]
           = (1.6605402 * 10000)/(4.074972637944947^{2})  [GPa]
           = (1.6605402 * 10000)/(4.074972637944947^{2})  [GPa] 
           = {1000.}  [GPa] 
           = N[%,20]  [GPa] 
           = {1000.}  [GPa] 
           =  1.  [TPa]

  [GPa] -> [Pa'] unit converter;
   (see "get_control_param()" in "program/control.c")

     press_unit_conversion = 1.0/( (1.6605402*10000.0)/
                              ((4.074972637944947)*(4.074972637944947)));
                           = 1/sys.pp2gpa

  [Pa'] -> [GPa] unit converter;
     sys.pp2gpa = (1.6605402*10000.0)/
                     ((4.074972637944947)*(4.074972637944947));

<< Physical constants >>------------------------------------------------

   * Boltzmann constant                                                    
        kB  = 1.380658e-23  [J/K]                                        
            => 1.380658e-5  [energy/temperature]                         
   * Elementally charge                                                    
        e   =  1.60217733e-19  [C]                                       
            => 1.0 (included in sqrt(sys.kk))                            

<< Coulomb force >>----------------------------------------------------

    * Dielectric constant in vacuum                                         
        epsilon0 = 8.854187817e-12 [F/m]                                  
                 = 8.854187817e-21 [C^2/(dyne cm^2)]                     

    * [C^2]/[F] = [J]
                                                                         
           1               e^2  Zi Zj  [C^2]                             
   U = --------------  --------------------------                        
       4 Pi Eps [F/m]    r [A] x 10^{-10}  [m]                           
                                                                         
         e^2        Zi Zj       [C^2]                                    
     = --------- ------------- -------                                   
       4 Pi Eps    r 10^{-10}    [F]                                     
                                                                         
                     Zi Zj                                               
     = 2.307079556 ---------  [energy]                                   
                      r                                                  
                                                                         
  In:=(1.60217733*10^(-19))^2*10^(18)/(4Pi 8.854187817*10^(-12)*10^(-10))  

       7.2479                                                              
  Out= ------                                                              
         Pi                                                                

  In:= N[%,10]                                                             

  Out= 2.307079556                                                         

  sys.kk = 1/(4Pi*Epsiron) [m/F]  <MKS-unit>
         = 1 [cm*erg/(stat.C^2)]  <cgs-unit>

  See also "program/rcprcl.c" and "program/real.c".

<< Energy unit conversion >>--------------------------------------------

   [energy] -> [energy][mol^{-1}] -> [J][mol^{-1}]                       
    N_A [avogadro constant] = 6.0221367 x 10^{23}                          
                            = 6.0221367e5 [in IEMD unit] 
                                                                         
   [energy] = 10^{-18} [J]                                               
   [energy][mol^{-1}] = 10^{-18} x 6.0221367 x 10^{23} [J][mol^{-1}]     
                      = 6.0221367 x 10^{23-18} [J][mol^{-1}]             
                      = 6.0221367 x 10^{5} [J][mol^{-1}]                 
                                                                         
   In NaCl crystal system, there are two atoms per one unit cell.  In    
   this case, the total number of NaCl unit (= avogadro unit = Mol. unit) 
   in the system is [sys.N/2].

  [energy] -> [J][mol^{-1}]  unit conversion 
  sys.perMol = AVOGADRO/((double)(sys.N/ctl.natoms_in_mol_unit)); 
  * AVOGADRO = 6.0221367e5 [in IEMD unit] 

  [energy] -> [temperature] = [K]                         
   sys.e2t = 2.0 / (3.0 * (double)(sys.N) * sys.kB);

--------------------------------------------------------------------------


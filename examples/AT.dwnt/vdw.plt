# set term pdfcairo enhanced size 8in, 6in
# set output "vdw.pdf"

set tics font "Times New Roman,18"  
set xlabel font "Times New Roman,18"
set ylabel font "Times New Roman,18"
set zlabel font "Times New Roman,18"
set key font "Times New Roman,18"   

set lmargin 10
set rmargin 10
set tmargin 5
set bmargin 5

# set  size   0.6, 1.0
set size 1.0,1.0
set  origin 0, 0

set title "vdW potential U and forces, F_{i} and F_{j}"
set title font"Arial,18"
set grid x
set key spacing 1.3

set xlabel "Distance of r_{j} from r_{i} at the origin (angstrom)" offset -1,0

plot [2:][-5:5]"vdw.dat" using 1:3 w l title "U/{/Symbol e}", \
     "vdw.dat" using 1:4 w l title "F_{i} {/Symbol s}/{/Symbol e}", \
     "vdw.dat" using 1:5 w l title "F_{j} {/Symbol s}/{/Symbol e}"

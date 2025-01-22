# set term pdfcairo enhanced size 8in, 6in
# set output "band.pdf"

set tics font "Times New Roman,16"   
set xlabel font "Times New Roman,20" 
set ylabel font "Times New Roman,20" 
set key font "Times New Roman,16"    
set lmargin 5
set rmargin 0
set tmargin 0
set bmargin 5

set  size   0.8, 0.7
set  origin 0.12, 0.1

## *** Plot range ***
set yrange [-770.5:-769.5]

set xzeroaxis
set grid x
set ylabel "Total Energy (kJ/mol)" offset -4,0
set xlabel "Time (psec)" offset 0,-2

unset key
unset grid

plot "out.csv" u ($2/1000):($5/1000) w l




pause -1
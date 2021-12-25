../../bin/corr.py -c 2 -t 1000 correlated_data.txt > tmp.txt
echo plot \"tmp.txt\" u 1:2:3 w e | gnuplot -p

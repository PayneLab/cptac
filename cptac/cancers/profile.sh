python -m cProfile -o test.pstats $1
echo "outputting to profile_output/$2"
gprof2dot -f pstats test.pstats | dot -Tpng -o profile_output/$2
eog profile_output/$2

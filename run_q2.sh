rm -f q2
g++ -O2 q2.cpp -o q2
((time ./q2) 2>&1) | tee out.txt
rm -f q2

g++ -O2 -DLOCAL q2.cpp -o q2
time ( ./q2 2>&1 | tee out.txt )
rm q2
ip1=$1":~/tri/"

echo $ip1
scp ./CMakeLists.txt ${ip1}
scp ./COPYING ${ip1}
scp -r ./src ${ip1}
scp -r ./include ${ip1}
scp -r ./test ${ip1}
scp ./CMakeLists.txt ${ip1}
scp ./COPYING ${ip1}

ssh $1 "cd ~/tri;cmake -DCMAKE_BUILD_TYPE=Release .;make -j4"
ssh $1 -t 'cd ~/run;bash --login'

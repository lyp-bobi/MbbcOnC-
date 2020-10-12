ip1="root@"$1":~/tri/"

echo $ip1
scp ./copy.sh ${ip1}
scp ./CMakeLists.txt ${ip1}
scp ./COPYING ${ip1}
scp -r ./src ${ip1}
scp -r ./include ${ip1}
scp -r ./test ${ip1}
scp ./CMakeLists.txt ${ip1}
scp ./COPYING ${ip1}

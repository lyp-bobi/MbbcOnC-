ip1=$1":~/tri/"

echo $ip1
ssh $1 "killall test-testx*"
rsync -v ./CMakeLists.txt ${ip1}
rsync -v ./COPYING ${ip1}
rsync -v -r ./src ${ip1}
rsync -v -r ./include ${ip1}
rsync -v -r ./test ${ip1}
rsync -v -r ./stxxl ${ip1}
rsync -v -r ./sqlite ${ip1}
rsync -v ./CMakeLists.txt ${ip1}
rsync -v ./COPYING ${ip1}

ssh $1 "mkdir ~/tri;mkdir ~/run;mkdir ~/build"
ssh $1 "rm ~/build/CMakeCache.txt"
ssh $1 "cd ~/build;cmake -DCMAKE_BUILD_TYPE=Release ../tri/;make -j4"
#ssh $1 "cd ~/build;cmake -DCMAKE_BUILD_TYPE=Debug ../tri/;make -j4"
# -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg .
ssh $1 -t 'cd ~/run;bash --login'

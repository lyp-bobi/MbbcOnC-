//
// Created by Chuang on 2020/11/24.
//

#include "testFuncs.h"

int main() {
    auto trajs = loadGLToTrajs("D://simp.csv");
    dumpToFile(trajs,"dumpedtraj.txt",500);
    return 0;
}
//
// Created by Chuang on 2021/1/20.
//

#include "testFuncs.h"

int main() {
    auto trajs = loadDumpedFiledToTrajs(testFileName("tdfilter.txt"));
    dumpToFile(trajs,testFileName("tdfilter.txt")+"ex");
    for(int i=1;i<10;i++) {
        auto repli = trajs;
        affine_transform(repli, xPoint(116.3972282409668, 39.90960456049752, 0),
                         0.2*i, xPoint(0, 0, 604800*i));
        dumpToFile_append(repli, testFileName("tdfilter.txt")+"ex");
    }
    return 0;
}
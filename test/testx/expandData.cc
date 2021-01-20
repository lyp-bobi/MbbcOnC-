//
// Created by Chuang on 2021/1/20.
//

#include "testFuncs.h"

int main() {
    auto trajs = loadDumpedFiledToTrajs("/root/tdfilter.txt");
    dumpToFile(trajs,"/root/tdexpand.txt");
    for(int i=1;i<10;i++) {
        auto repli = trajs;
        affine_transform(repli, xPoint(116.3972282409668, 39.90960456049752, 0),
                         0.2*i, xPoint(5*i, 5*i, 0));
        dumpToFile_append(repli, "/root/tdexpand.txt");
    }
    return 0;
}
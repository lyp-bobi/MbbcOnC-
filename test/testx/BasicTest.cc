//
// Created by Chuang on 2020/10/28.
//
#include "storagemanager/xStore.h"
using namespace xRTreeNsp;
int main(){
    xStore x("test", "D://TRI-framework/dumpedtraj.txt",true);
    auto stat = trajStat::instance();
    stat->bt = 100000;
    CUTFUNC f=xTrajectory::GLL;
    buildMBRRTree(&x,f);
}
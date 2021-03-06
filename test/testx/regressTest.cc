//
// Created by Chuang on 2020/10/28.
//
#include "testFuncs.h"

int main(){
    struct str{
        char a[100];
    };
    cerr<<sizeof(str)<<endl;
    string target = "tdexpand.datas";
    xStore x(target, testFileName(target), true,true);
    double seglen[]={600,900,1200,1800};
    xCylinder query(xPoint(40,116.327,6516),0.0001,6516,9516,2);
    xTrajectory tj,tj2;
    x.loadTraj(tj, xStoreEntry(8,0,1000));
    x.loadTraj(tj2,xStoreEntry(1858,0, 10000));
//    tj.getMinimumDistance(tj2);
//    x.loadTraj(tj, xStoreEntry(4337,0,1000));
//    tj.intersectsxCylinder(query);
//    cerr<<tj.intersectsxCylinder(query);
    tjstat->bt=50;
    xSBB s;
    s.loadFromString("0 1 0 39.984438 116.347402 9194.000000 39.998213 116.365543 10194.000000 0.043406 0.000149");
    query.checkRel(s);
    MyVisitor vis;
//    tj.loadFromString("117.110025,40.151262,1513476.000000 117.110610,40.151180,1513827.000000 117.110370,40.151440,1514427.000000 117.109900,40.151280,1515026.000000 117.110030,40.151520,1515626.000000 117.109830,40.151620,1516226.000000 117.109810,40.151600,1516826.000000 117.109840,40.151600,1517426.000000 117.109834,40.151588,1517776.000000");
//    x.loadTraj(tj2,xStoreEntry(23411,0,100000));
//    xTrajectory tj3;
//    tj2.getPartialxTrajectory(1513476,1517776,tj3);
//    cerr<<tj3<<endl;
//    cerr<<tj.getMinimumDistance(tj3)<<endl;
//    tj2.getPartialxTrajectory(1513476,1513613,tj3);
//    cerr<<tj3.getMinimumDistance(tj)<<endl;
//    tj2.getPartialxTrajectory(1513613,1514813,tj3);
//    cerr<<tj3.getMinimumDistance(tj)<<endl;
//    tj2.getPartialxTrajectory(1514813,1516013,tj3);
//    cerr<<tj3.getMinimumDistance(tj)<<endl;
//    tj2.getPartialxTrajectory(1516013,1517213,tj3);
//    cerr<<tj3.getMinimumDistance(tj)<<endl;
//    tj2.getPartialxTrajectory(1517213,1517776,tj3);
//    cerr<<tj3.getMinimumDistance(tj)<<endl;

    try {
        for(int i=0;i<testtime;i++) {
//            tj=x.randomSubtraj(4300);
            tj.loadFromString("39.999976,116.326565,8016.000000 40.000032,116.326175,8021.000000 40.000038,116.326154,8026.000000 40.000035,116.326081,8031.000000 39.996877,116.326645,8151.000000 39.996825,116.326536,8156.000000");
            double res1 = 0;
            double res2 = 0;
            for (auto querylen:seglen) {
                tjstat->bt = querylen;
                auto r = buildMBCRTreeWP(&x, xTrajectory::ISS, querylen);
//            r->intersectsWithQuery(query, vis);
                r->nearestNeighborQuery(6, tj, vis);

                std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
                if(res1==0&&res2==0) {
                    res1 = vis.m_resultGet;
                    res2 = vis.m_lastResult;
                }else{
                    if(vis.m_resultGet!=res1 || vis.m_lastResult!=res2){
                        cerr<<tj<<endl;
                        return 1;
                    }
                }
                vis.clear();
//            auto r2 = buildMBRRTreeWP(&x, xTrajectory::ISS,querylen);
//            r2->intersectsWithQuery(query, vis);
//            r2->nearestNeighborQuery(6, tj, vis);
//            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
//            vis.clear();
            }
//        {
//            auto r = buildTBTreeWP(&x);
//            r->intersectsWithQuery(query, vis);
//            r->nearestNeighborQuery(6, tj, vis);
//            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
//            vis.clear();
//        }
        {
            auto r = buildSTRTreeWP(&x);
//            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(6, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
//        for (auto querylen:querylens) {
//            tjstat->bt = querylen;
//            auto r = buildMBCRTreeWP(&x, xTrajectory::ISS,querylen);
//            r->intersectsWithQuery(query, vis);
//            r->nearestNeighborQuery(6, tj, vis);
//            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
//            vis.clear();
//            auto r2 = buildMBRRTreeWP(&x, xTrajectory::ISS,querylen);
//            r2->intersectsWithQuery(query, vis);
//            r2->nearestNeighborQuery(6, tj, vis);
//            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
//            vis.clear();
//        }
//        {
//            auto r = buildTBTreeWP(&x);
//            r->intersectsWithQuery(query, vis);
//            r->nearestNeighborQuery(6, tj, vis);
//            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
//            vis.clear();
//        }
//        {
//            auto r = buildSTRTreeWP(&x);
//            r->intersectsWithQuery(query, vis);
//            r->nearestNeighborQuery(6, tj, vis);
//            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
//            vis.clear();
//        }
        }
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
}
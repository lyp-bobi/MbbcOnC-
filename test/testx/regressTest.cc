//
// Created by Chuang on 2020/10/28.
//
#include "testFuncs.h"
int main(){
    xStore x("test", "D://TRI-framework/dumpedtraj.txt",true);

    double querylens[]={50,100,200,400,800,1600,10000};
    xCylinder query(xPoint(40,116.327,6516),0.0001,6516,9516,2);
    xTrajectory tj,tj2;
    x.loadTraj(tj, xStoreEntry(1850,0,1000));
    x.loadTraj(tj2,xStoreEntry(1858,0, 10000));
    tj.getMinimumDistance(tj2);
//    x.loadTraj(tj, xStoreEntry(4337,0,1000));
    tj.intersectsxCylinder(query);
//    cerr<<tj.intersectsxCylinder(query);
    tjstat->bt=50;
    xSBB s;
    s.loadFromString("0 1 0 39.984438 116.347402 9194.000000 39.998213 116.365543 10194.000000 0.043406 0.000149");
    query.checkRel(s);
    MyVisitor vis;
    try {
        bUsingSBBD = true;
        for (auto querylen:querylens) {
            tjstat->bt = querylen;
            auto r = buildMBCRTreeWP(&x, xTrajectory::ISS,querylen);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
            auto r2 = buildMBRRTreeWP(&x, xTrajectory::ISS,querylen);
            r2->intersectsWithQuery(query, vis);
            r2->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        {
            auto r = buildTBTreeWP(&x);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        {
            auto r = buildSTRTreeWP(&x);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        bUsingSBBD = false;
        for (auto querylen:querylens) {
            tjstat->bt = querylen;
            auto r = buildMBCRTreeWP(&x, xTrajectory::ISS,querylen);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
            auto r2 = buildMBRRTreeWP(&x, xTrajectory::ISS,querylen);
            r2->intersectsWithQuery(query, vis);
            r2->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        {
            auto r = buildTBTreeWP(&x);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
        {
            auto r = buildSTRTreeWP(&x);
            r->intersectsWithQuery(query, vis);
            r->nearestNeighborQuery(5, tj, vis);
            std::cerr << vis.m_resultGet << " " << vis.m_lastResult << endl;
            vis.clear();
        }
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
}
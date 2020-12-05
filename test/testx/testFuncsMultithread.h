//
// Created by Chuang on 2020/11/30.
//

#ifndef SPATIALINDEX_TESTFUNCSMULTITHREAD_H
#define SPATIALINDEX_TESTFUNCSMULTITHREAD_H

#include "testFuncs.h"
#include <thread>

#define NUMCORE 4
#define NUMTHREAD (2*NUMCORE)

using namespace std;

struct queryRet{
    double time=0;
    double indexVisit=0;
    double leafVisit=0;
    double leaf1 = 0;
    double leaf2 = 0;
    double indexIO=0;
    double trajIO=0;
    queryRet operator+(const queryRet& r) const{
        queryRet res;
        res.time= time+r.time;
        res.indexVisit= indexVisit+r.indexVisit;
        res.leafVisit= leafVisit+r.leafVisit;
        res.leaf1= leaf1+r.leaf1;
        res.leaf2= leaf2+r.leaf2;
        res.indexIO= indexIO+r.indexIO;
        res.trajIO= trajIO+r.trajIO;
        return res;
    }
    string toString() const{
        stringstream s;
        s<<time<<"\t"<<indexVisit<<"\t"<<leafVisit<<"\t"<<indexIO<<"\t"<<trajIO<<"\n";
        return s.str();
    }
};


queryRet average(vector<queryRet> &sum){
    queryRet res;
    for(auto &s:sum) res=res+s;
    res.time/=sum.size();
    res.indexVisit/=sum.size();
    res.leafVisit/=sum.size();
    res.leaf1/=sum.size();
    res.leaf2/=sum.size();
    res.indexIO/=sum.size();
    res.trajIO/=sum.size();
    return res;
}

struct queryInp{
    xRTree * tree = nullptr;
    vector<xTrajectory> knn_queries;
    int nnk =5;
};


void kNNQueryBatchThread(queryInp inp, queryRet *res) {
    xStore * ts = inp.tree->m_ts;
    ts->cleanStatistic();
    ts->flush();
    drop_cache(3);
    int num = inp.knn_queries.size();
    MyVisitor vis;
    vis.ts = ts;
    auto start = std::chrono::system_clock::now();
    double rad = 0;
    int indio = 0;
    std::vector<int> indios;
    for (int i = 0; i < inp.knn_queries.size(); i++) {
        vis.m_query = (IShape *) &(inp.knn_queries.at(i));
        inp.tree->nearestNeighborQuery(inp.nnk, inp.knn_queries.at(i), vis);
        rad += vis.m_lastDist;
        indios.emplace_back(ts->m_indexIO - indio);
        indio = ts->m_indexIO;
    }
    rad /= inp.knn_queries.size();
    double time;
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    time = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    res->time = time/num;
    res->indexVisit = 1.0*vis.m_indexvisited/num;
    res->leafVisit = 1.0 * vis.m_leafvisited/ num;
    res->leaf1 = 1.0 * ts->m_leaf1/num;
    res->leaf2 = 1.0 * ts->m_leaf2 / num;
    res->indexIO = 1.0 * ts->m_indexIO / num;
    res->trajIO = 1.0 * ts->m_trajIO / num;
    return;
}

#endif //SPATIALINDEX_TESTFUNCSMULTITHREAD_H

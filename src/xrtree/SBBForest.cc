//
// Created by Chuang on 2021/1/30.
//
#include "SBBForest.h"
using namespace xRTreeNsp;
xRTree * SBBForest::chooseTree(double qt) {
    for(auto &s:m_trees){
        if(s.first.first<=qt && s.first.second>=qt){
            return s.second;
        }
    }
    cerr<<"Not supported qt "<<qt<<endl;
    return NULL;
}

void SBBForest::freeall() {
    for(auto &s:m_trees){
        delete s.second;
    }
}

void SBBForest::insert(double l, double h, xRTree *r) {
    m_trees[make_pair(l,h)] = r;
}

void SBBForest::intersectsWithQuery(const xCylinder &query, IVisitor &v) {
    xRTree* r = chooseTree(query.m_endTime - query.m_startTime);
    r->intersectsWithQuery(query,v); 
}

void SBBForest::nearestNeighborQuery(uint32_t k, const xTrajectory &query, IVisitor &v) {
    xRTree* r = chooseTree(query.m_endTime() - query.m_startTime());
    r->nearestNeighborQuery(k,query,v);
}

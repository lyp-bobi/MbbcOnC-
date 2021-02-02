//
// Created by Chuang on 2021/1/30.
//

#ifndef SPATIALINDEX_SBBFOREST_H
#define SPATIALINDEX_SBBFOREST_H
#include "xRTree.h"
namespace SpatialIndex {
    namespace xRTreeNsp {
        class SBBForest {
        public:
            map<pair<double, double>, xRTree *> m_trees;

            xRTree *chooseTree(double qt);
            void insert(double l,double h, xRTree* r);
            void freeall();
            void intersectsWithQuery(const xCylinder& query, IVisitor& v);
            void nearestNeighborQuery(uint32_t k, const xTrajectory& query, IVisitor& v);
        };
    }
}

#endif //SPATIALINDEX_SBBFOREST_H

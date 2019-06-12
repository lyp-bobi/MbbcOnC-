//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "Rtree.h"
#include "storagemanager/TrajStore.h"

class rtreeSegmentStream:public baseSegmentStream{
public:
    rtreeSegmentStream(TrajStore *ts): baseSegmentStream(ts){}
    IData* constructData(id_type id,Region mbr,MBC mbc) override {
        RTree::Data* d=new RTree::Data(0, nullptr, mbr, id);
        return d;
    }
};


ISpatialIndex* SpatialIndex::RTree::createAndBulkLoadNewRTreeWithTrajStore(IStorageManager *tsm,
                                                                           uint32_t indexCapacity, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    TrajStore *ts= dynamic_cast<TrajStore*>(tsm);
    assert(ts!= nullptr);
    auto dataStream=new rtreeSegmentStream(ts);
    ISpatialIndex* tree= createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR,*dataStream,*ts,0.5,indexCapacity,2,dimension,SpatialIndex::RTree::RV_RSTAR,indexIdentifier);
    RTree* r= dynamic_cast<RTree*>(tree);
    r->m_DataType=TrajectoryType;
    r->m_bUsingTrajStore=true;
    r->m_ts=ts;
    return r;
}
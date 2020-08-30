//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "RTree.h"
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
                                                                           uint32_t pageSize, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    TrajStore *ts= static_cast<TrajStore*>(tsm);
    assert(ts!= nullptr);
    auto dataStream=new rtreeSegmentStream(ts);
    int indexCapacity=std::floor((pageSize-56)/56);
    int leafCapacity=std::floor((pageSize-56)/88);
    ISpatialIndex* tree;
    if (ts->m_stream!= nullptr) {
        tree = createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, *ts->m_stream, *ts, 0.9,
                                            indexCapacity, leafCapacity, dimension,
                                            SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
    } else {
        auto dataStream = new rtreeSegmentStream(ts);
        tree = createAndBulkLoadNewRTree(SpatialIndex::RTree::BLM_STR, *dataStream, *ts, 0.9,
                                            indexCapacity, leafCapacity, dimension,
                                            SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
        delete dataStream;
    }

    delete dataStream;
    RTree* r= static_cast<RTree*>(tree);
    r->m_DataType=TrajectoryType;
    r->m_bUsingTrajStore=true;
    r->m_ts=ts;
    return r;
}
//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "MBCRtree.h"
#include "storagemanager/TrajStore.h"

class mbcrtreeSegmentStream:public baseSegmentStream{
public:
    mbcrtreeSegmentStream(TrajStore *ts): baseSegmentStream(ts){}
    IData* constructData(id_type id,Region mbr,MBC mbc) override {
        MBCRTree::Data* d=new MBCRTree::Data(0, nullptr, mbc, id);
        return d;
    }
};


ISpatialIndex* SpatialIndex::MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(IStorageManager *tsm,
                                                                           uint32_t indexCapacity, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    TrajStore *ts= dynamic_cast<TrajStore*>(tsm);
    assert(ts!= nullptr);
    auto dataStream=new mbcrtreeSegmentStream(ts);
    ISpatialIndex* tree= createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR,*dataStream,*ts,0.5,indexCapacity,2,dimension,SpatialIndex::MBCRTree::RV_RSTAR,indexIdentifier);
    MBCRTree* r= dynamic_cast<MBCRTree*>(tree);
    r->m_DataType=TrajectoryType;
    r->m_bUsingTrajStore=true;
    r->m_ts=ts;
    return r;
}
//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "MBCRTree.h"
#include "storagemanager/TrajStore.h"

class mbcrtreeSegmentStream:public baseSegmentStream{
public:
    mbcrtreeSegmentStream(TrajStore *ts): baseSegmentStream(ts){}
    IData* constructData(id_type id,Region mbr,MBC mbc) override {
        MBCRTree::Data* d=new MBCRTree::Data(0, nullptr, mbc,mbr, id);
        return d;
    }
};


ISpatialIndex* SpatialIndex::MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(IStorageManager *tsm,
                                                                           uint32_t indexCapacity, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    TrajStore *ts= static_cast<TrajStore*>(tsm);
    auto dataStream=new mbcrtreeSegmentStream(ts);
    ISpatialIndex* tree= createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR,*dataStream,*ts,0.9,indexCapacity,indexCapacity,dimension,SpatialIndex::MBCRTree::RV_RSTAR,indexIdentifier);
    MBCRTree* r= static_cast<MBCRTree*>(tree);
    return r;
}

ISpatialIndex* SpatialIndex::MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(IStorageManager *tsm,
                                                                                 uint32_t indexCapacity, uint32_t dimension,
                                                                                 id_type &indexIdentifier) {
    TrajStore *ts= static_cast<TrajStore*>(tsm);
    auto dataStream=new mbcrtreeSegmentStream(ts);
    ISpatialIndex* tree= createAndBulkLoadNewMBCRTree(SpatialIndex::MBCRTree::BLM_STR,*dataStream,*ts,0.9,indexCapacity,indexCapacity,dimension,SpatialIndex::MBCRTree::RV_RSTAR,indexIdentifier,true);
    MBCRTree* r= static_cast<MBCRTree*>(tree);
    return r;
}
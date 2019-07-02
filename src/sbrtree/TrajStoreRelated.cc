//
// Created by Chuang on 2019/6/12.
//

#include <spatialindex/SpatialIndex.h>
#include "SBRTree.h"
#include "storagemanager/TrajStore.h"

class SBRTreeSegmentStream:public baseSegmentStream{
public:
    SBRTreeSegmentStream(TrajStore *ts): baseSegmentStream(ts){}
    IData* constructData(id_type id,Region mbr,MBC mbc) override {
        SBR sbr;
        sbr.getFromMBC(mbc,0,50);
//        std::cerr<<sbr<<"\n";
        SBRTree::Data* d=new SBRTree::Data(0, nullptr, mbc,sbr, id);
        return d;
    }
};


ISpatialIndex* SpatialIndex::SBRTree::createAndBulkLoadNewSBRTreeWithTrajStore(IStorageManager *tsm,
                                                                           uint32_t indexCapacity, uint32_t dimension,
                                                                           id_type &indexIdentifier) {
    TrajStore *ts= static_cast<TrajStore*>(tsm);
    auto dataStream=new SBRTreeSegmentStream(ts);
    ISpatialIndex* tree= createAndBulkLoadNewSBRTree(SpatialIndex::SBRTree::BLM_STR,*dataStream,*ts,0.9,indexCapacity,indexCapacity,dimension,SpatialIndex::SBRTree::RV_RSTAR,indexIdentifier);
    SBRTree* r= static_cast<SBRTree*>(tree);
    r->m_DataType=TrajectoryType;
    r->m_bUsingTrajStore=true;
    r->m_ts=ts;
    return r;
}
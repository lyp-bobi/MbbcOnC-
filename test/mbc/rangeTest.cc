//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main(){
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs();
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        auto stat=trajStat::instance();
        int maxseg = 0;
        double avgSegLen=100;
//        double segLenParas[]={10,20,50,80,100,200,300,400,500,750,1000,1500,2000,2500,3000};
        double segLenParas[]={60,70,80,90,100,140,150,170,180,200};
        double queryLenParas[]={0,0,3600};
        std::cerr<<"Starting range test\n"<<"Segmentation lengths are:";
        for(auto p:segLenParas) std::cerr<<p<<"\t";
        std::cerr<<"\nQuery lengths are:";
        for(auto p:queryLenParas) std::cerr<<p<<"\t";
        std::cerr<<"\n";
        vector<vector<IShape *>> querySet;
        for(auto queryLen:queryLenParas) {
            vector<IShape *> queries;
            for (int i = 0; i < 1000; i++) {
                double t = int(random(stat->mint, stat->maxt - queryLen));
                double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
//                double pHigh[3] = {pLow[0] + random(0.05,0.1),
//                                   pLow[1] + random(0.05,0.1), t + queryLen};
//                Region *rg = new Region(pLow, pHigh, 3);
                Cylinder *rg = new Cylinder(pLow,random(0.1,0.5),t,t+queryLen,2);
                queries.emplace_back(rg);
            }
            querySet.emplace_back(queries);
            queries.clear();
        }
        for (double segLen:segLenParas) {
            maxseg=300;
            segs.clear();
            int totallen = 0, totalseg = 0;
            string name0 ="name0", name1 ="name1", name2 = "name2";
            id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                    *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                    *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
            // Create a new storage manager with the provided base name and a 4K page size.
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

            for (const auto &traj:trajs) {
                totallen += traj.second.m_points.size();
//                if(traj.second.m_endTime()-traj.second.m_startTime()<10){
//                    std::cerr<<traj.second<<endl;
//                }
                auto seg = traj.second.getSegments(segLen);
                totalseg += seg.size();
                maxseg = std::max(int(seg.size()), maxseg);
                segs.emplace_back(make_pair(traj.first, seg));
            }
            std::cerr<<"seg is"<<totalseg<<"\n";
            avgSegLen=double(totallen)/totalseg;
            std::cerr<<"segments' average length is "<<totallen*1.0/totalseg<<"\n";

            TrajStore *ts1 = new TrajStore(diskfile1, 4096, maxseg+1);
            TrajStore *ts2 = new TrajStore(diskfile2, 4096, maxseg+1);

            ts1->loadSegments(segs,true);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            ts2->loadSegments(segs,true);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            segs.clear();
            segs.swap(emptyseg);
            std::cerr<<"Seg len:"<<segLen<<"\n";
            for(const auto &qs:querySet) {
                rangeQueryBatch(rc, qs, ts2);
                rangeQueryBatch(r, qs, ts1);
            }
//            std:cerr<<*r<<*rc<<endl;

            delete r;
            delete ts1;
            delete rc;
            delete ts2;
            delete file0;
            delete file1;
            delete file2;
            delete diskfile0;
            delete diskfile1;
            delete diskfile2;
        }
        for(auto &qs:querySet){
            for(auto &shape:qs){
                delete(shape);
            }
        }
    }
    catch (Tools::Exception& e)
    {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    catch(...){

    }
    return 0;
}
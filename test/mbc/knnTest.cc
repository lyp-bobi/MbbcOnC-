//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main(){
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        auto stat=trajStat::instance();
        int maxseg = 0;
        double avgSegLen=100;
        double segLenParas[]={300,500,800,1000,1500,2000,2500,3000};
        double queryLenParas[]={900,3600,18000};
        std::cerr<<"Starting knn test\n"<<"Segmentation lengths are:";
        for(auto p:segLenParas) std::cerr<<p<<"\t";
        std::cerr<<"\nQuery lengths are:";
        for(auto p:queryLenParas) std::cerr<<p<<"\t";
        std::cerr<<"\n";
        vector<vector<IShape *>> querySet;
        for (auto queryLen:queryLenParas) {
            vector<IShape *> queries;
            for (int i = 0; i < 100; i++) {
                auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
                Trajectory *concate = new Trajectory();
                double ts = std::max(ori->m_startTime(), random(ori->m_startTime(), ori->m_endTime() - queryLen));
                ori->getPartialTrajectory(ts, ts + queryLen, *concate);
                if (!concate->m_points.empty())
                    queries.emplace_back(concate);
            }
            querySet.emplace_back(queries);
            queries.clear();
        }
        for (double segLen:segLenParas) {
                maxseg=300;
                segs.clear();
                int totallen = 0, totalseg = 0;
            for (const auto &traj:trajs) {
                totallen += traj.second.m_points.size();
                auto seg = traj.second.getSegments(segLen);
                totalseg += seg.size();
                maxseg = std::max(int(seg.size()), maxseg);
                segs.emplace_back(make_pair(traj.first, seg));
            }
            avgSegLen=double(totallen)/totalseg;
            std::cerr<<"segments' average length is "<<totallen*1.0/totalseg<<"\n";

            string name0 ="name0", name1 ="name1", name2 = "name2";
            id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                    *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                    *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
            // Create a new storage manager with the provided base name and a 4K page size.
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 100, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

            TrajStore *ts1 = new TrajStore(file1, 4096, maxseg+1);
            ts1->loadSegments(segs, true);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(file2, 4096, maxseg+1);
            ts2->loadSegments(segs,true);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            //kNN
            segs.clear();
            segs.swap(emptyseg);
            disttype=0;
            std::cerr<<"Seg len:"<<segLen<<"\n";
            for(const auto &qs:querySet) {
                kNNQueryBatch(r, qs, ts1);
                kNNQueryBatch(rc, qs, ts2);
            }
////                    simpli= false;
//                    disttype=1;
//                    kNNQueryBatch(r, queries, ts1);
//                    kNNQueryBatch(rc, queries, ts2);


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
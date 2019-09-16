//
// Created by Chuang on 2019/9/3.
//

//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main(){
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTToTrajs("D://00.txt");
        vector<pair<id_type, vector<Trajectory>>> segs1,segs2,segs3,segs4;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        for (double segpara = 100; segpara<=2000 ; segpara+=300) {
            maxseg=100;
            segs1.clear();segs2.clear();segs3.clear();segs4.clear();
            segs1.swap(emptyseg);segs2.swap(emptyseg);segs3.swap(emptyseg);segs4.swap(emptyseg);
            std::cerr<<segpara<<"\n";
            for (const auto &traj:trajs) {
                auto seg = traj.second.getSegments(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            for (const auto &traj:trajs) {
                auto seg = traj.second.getStaticSegments(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs2.emplace_back(make_pair(traj.first, seg));
            }
            for (const auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments(segpara/5);
                maxseg = std::max(int(seg.size()), maxseg);
                segs3.emplace_back(make_pair(traj.first, seg));
            }
            for (const auto &traj:trajs) {
                auto seg = traj.second.getRDPSegments(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs4.emplace_back(make_pair(traj.first, seg));
            }
            string name0 ="name0", name1 ="name1", name2 = "name2",name3="name3",name4="name4";
            id_type indexIdentifier0, indexIdentifier1, indexIdentifier2,indexIdentifier3,indexIdentifier4;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                    *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                    *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096),
                    *diskfile3 = StorageManager::createNewDiskStorageManager(name3, 4096),
                    *diskfile4 = StorageManager::createNewDiskStorageManager(name4, 4096);
            // Create a new storage manager with the provided base name and a 4K page size.
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 100, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false),
                    *file3 = StorageManager::createNewRandomEvictionsBuffer(*diskfile3, 10, false),
                    *file4 = StorageManager::createNewRandomEvictionsBuffer(*diskfile4, 10, false);

            TrajStore *ts1 = new TrajStore(file1, 4096, maxseg+1);
            ts1->loadSegments(segs1,true);
            ISpatialIndex *r1 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(file2, 4096, maxseg+1);
            ts2->loadSegments(segs2,true);
            ISpatialIndex *r2 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            TrajStore *ts3 = new TrajStore(file3, 4096, maxseg+1);
            ts3->loadSegments(segs3,true);
            ISpatialIndex *r3 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts3, 4096, 3, indexIdentifier3);

            TrajStore *ts4 = new TrajStore(file4, 4096, maxseg+1);
            ts4->loadSegments(segs4,true);
            ISpatialIndex *r4 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts4, 4096, 3, indexIdentifier4);

            std::cerr<<"traj Store Loaded\n";
            //kNN
            segs1.clear();segs2.clear();segs3.clear();segs4.clear();
            segs1.swap(emptyseg);segs2.swap(emptyseg);segs3.swap(emptyseg);segs4.swap(emptyseg);

            vector<IShape *> queries;
//            double segattri[]={900,3600,1000000};
            double segattri[]={3600};
            auto stat=trajStat::instance();
            for (auto queryLen:segattri) {
//                for(int thek=1;thek<=21;thek+=5){
                for (int thek = 5; thek == 5; thek++) {
                    for (int i = 0; i < 200; i++) {
                        auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
                        Trajectory *concate = new Trajectory();
                        double ts = std::max(ori->m_startTime(),random(ori->m_startTime(), ori->m_endTime() - queryLen));
                        ori->getPartialTrajectory(ts, ts + queryLen, *concate);
                        if (!concate->m_points.empty())
                            queries.emplace_back(concate);
                    }
//                    cerr << "=================\n\n";
//                    std::cerr << "Querying with segmenting len " << seglen <<
//                              ", querying len " << queryLen << ", NN's k" << thek << "\n";
                    cerr<< segpara<<"\t"<<queryLen<<"\n";
//                    simpli=true;
                    disttype=0;
                    kNNQueryBatch(r1, queries, ts1);
                    kNNQueryBatch(r2, queries, ts2);
                    kNNQueryBatch(r3, queries, ts3);
                    kNNQueryBatch(r4, queries, ts4);

                    cerr << "\n";
                    for(auto &tras:queries){
                        delete tras;
                    }
                    queries.clear();
                }
            }
            delete r1;
            delete r2;
            delete r3;
            delete r4;
            delete ts1;
            delete ts2;
            delete ts3;
            delete ts4;
            delete file0;
            delete file1;
            delete file2;
            delete file3;
            delete file4;
            delete diskfile0;
            delete diskfile1;
            delete diskfile2;
            delete diskfile3;
            delete diskfile4;
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
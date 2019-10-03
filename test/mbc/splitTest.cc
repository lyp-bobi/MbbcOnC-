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
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTToTrajs("D://00.txt");
        vector<pair<id_type, vector<Trajectory>>> segs0,segs1,segs2,segs3,segs4;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        double para[]={300,600,900,1200,1500,2000,2500,3000};
        auto stat=trajStat::instance();
        double queryLen=3600;
        vector<IShape *> queries;
        for (int i = 0; i < 200; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = std::max(ori->m_startTime(),random(ori->m_startTime(), ori->m_endTime() - queryLen));
            ori->getPartialTrajectory(ts, ts + queryLen, *concate);
            if (!concate->m_points.empty())
                queries.emplace_back(concate);
        }
//        for (int i = 0; i < 1000; i++) {
//            double t = int(random(stat->mint, stat->maxt - queryLen));
//            double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
//            double pHigh[3] = {pLow[0] + random(stat->Dx / 40, stat->Dx * 3 / 40),
//                               pLow[1] + random(stat->Dy / 40, stat->Dy * 3 / 40), t + queryLen};
//            Region *rg = new Region(pLow, pHigh, 3);
//            queries.emplace_back(rg);
//        }
        for (double segpara:para) {
            maxseg=500;
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

            segs0.clear();segs1.clear();segs2.clear();segs3.clear();segs4.clear();
            segs0.swap(emptyseg);segs1.swap(emptyseg);segs2.swap(emptyseg);segs3.swap(emptyseg);segs4.swap(emptyseg);
            std::cerr<<segpara<<"\n";
            std::cerr<<"start splitting\n";

            for (const auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments(int(std::ceil(segpara/stat->tl)+2));
                maxseg = std::max(int(seg.size()), maxseg);
                segs3.emplace_back(make_pair(traj.first, seg));
            }
            TrajStore *ts3 = new TrajStore(file3, 4096, maxseg+1);
            ts3->loadSegments(segs3,true);
            segs3.clear();segs3.swap(emptyseg);
            ISpatialIndex *r3 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts3, 4096, 3, indexIdentifier3);
            kNNQueryBatch(r3, queries, ts3);
            delete r3;delete ts3;delete file3;delete diskfile3;

            for (const auto &traj:trajs) {
                auto seg = traj.second.getStaticSegments(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            TrajStore *ts1 = new TrajStore(file1, 4096, maxseg+1);
            ts1->loadSegments(segs1,true);
            segs1.clear();segs1.swap(emptyseg);
            ISpatialIndex *r1 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);
            kNNQueryBatch(r1, queries, ts1);
            delete r1;delete ts1;delete file1;delete diskfile1;


            for (const auto &traj:trajs) {
                auto seg = traj.second.getStaticSegmentsCut(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs2.emplace_back(make_pair(traj.first, seg));
            }
            TrajStore *ts2 = new TrajStore(file2, 4096, maxseg+1);
            ts2->loadSegments(segs2,true);
            segs2.clear();segs2.swap(emptyseg);
            ISpatialIndex *r2 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            kNNQueryBatch(r2, queries, ts2);
//            rangeQueryBatch(r2, queries, ts2);
            delete r2;delete ts2;delete file2;delete diskfile2;

            for (const auto &traj:trajs) {
                auto seg = traj.second.getSegments(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs0.emplace_back(make_pair(traj.first, seg));
            }
            TrajStore *ts0 = new TrajStore(file0, 4096, maxseg+1);
            ts0->loadSegments(segs0,true);
            segs0.clear();segs0.swap(emptyseg);
            ISpatialIndex *r0 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts0, 4096, 3, indexIdentifier0);

            kNNQueryBatch(r0, queries, ts0);
//            rangeQueryBatch(r0, queries, ts0);
            delete r0;delete ts0;delete file0;delete diskfile0;



            for (const auto &traj:trajs) {
                auto seg = traj.second.getRDPSegments(segpara);
                maxseg = std::max(int(seg.size()), maxseg);
                segs4.emplace_back(make_pair(traj.first, seg));
            }
            TrajStore *ts4 = new TrajStore(file4, 4096, maxseg+1);
            ts4->loadSegments(segs4,true);
            segs4.clear();segs4.swap(emptyseg);
            ISpatialIndex *r4 = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts4, 4096, 3, indexIdentifier4);
            kNNQueryBatch(r4, queries, ts4);
            delete r4;delete ts4;delete file4;delete diskfile4;

            std::cerr<<"split finished\n";
//            kNNQueryBatch(r0, queries, ts0);
//            kNNQueryBatch(r1, queries, ts1);
//            kNNQueryBatch(r2, queries, ts2);
//            kNNQueryBatch(r3, queries, ts3);
//            kNNQueryBatch(r4, queries, ts4);
//            delete r0;delete ts0;delete file0;delete diskfile0;
//            delete r1;delete ts1;delete file1;delete diskfile1;
//            delete r2;delete ts2;delete file2;delete diskfile2;
//            delete r2;delete ts3;delete file3;delete diskfile3;
//            delete r3;delete ts4;delete file4;delete diskfile4;
        }
        for(auto &tras:queries){
            delete tras;
        }
        queries.clear();
    }
    catch (Tools::Exception& e)
    {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}
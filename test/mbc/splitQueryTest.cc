//
// Created by Chuang on 2019/9/19.
//

//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main() {
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        double avgSegLen = 100;
//        for (double segpara = 0.1; avgSegLen>10 ; segpara/=2) {
//        double para[]={100,300,500,800,1000,1500,2000,2500,3000};
//        double para[]={1200};

        vector<IShape *> queries;
//            double segattri[]={900,3600,1000000};
//        double segattri[]={50,1000};
        double queryLen = 3600;
        auto stat = trajStat::instance();
        for (int i = 0; i < 1; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - queryLen/2;
            ori->getPartialTrajectory(ts, ts + queryLen, *concate);
            if (!concate->m_points.empty())
                queries.emplace_back(concate);
            std::cerr << "query is" << *concate << "\n";
        }
//        double para[]={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
        double para[] = {300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000};
        for (double segpara:para) {
            maxseg = 300;
            segs.clear();
            int totallen = 0, totalseg = 0;
            for (const auto &traj:trajs) {
                totallen += traj.second.m_points.size();
                auto seg = traj.second.getSegments(segpara);
                totalseg += seg.size();
                maxseg = std::max(int(seg.size()), maxseg);
                segs.emplace_back(make_pair(traj.first, seg));
            }
            avgSegLen = double(totallen) / totalseg;
            std::cerr << "segments' average length is " << totallen * 1.0 / totalseg << "\n";

            string name0 = "name0", name1 = "name1", name2 = "name2";
            id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                    *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                    *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
            // Create a new storage manager with the provided base name and a 4K page size.
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

            TrajStore *ts1 = new TrajStore(name1, file1, 4096, maxseg + 1);
            ts1->loadSegments(segs, true);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(name2, file2, 4096, maxseg + 1);
//            ts2->loadSegments(segs,true);
//            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            //kNN
            segs.clear();
            segs.swap(emptyseg);
            for (int thek = 5; thek == 5; thek++) {
                disttype = 0;
                kNNQueryBatch(r, queries, ts1, 5, true);
//                kNNQueryBatch(rc, queries, ts2);
                cerr << "\n";
            }
//            for(int i=0;i<1;i++){
//                double t = int(random(stat->mint, stat->maxt));
//                double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
//                double pHigh[3] = {pLow[0] + random(stat->Dx/40,stat->Dx*3/40), pLow[1] + random(stat->Dy/40,stat->Dy*3/40), t};
//                Region *rg = new Region(pLow, pHigh, 3);
//                queries.emplace_back(rg);
//            }
//            rangeQueryBatch(r,queries,ts1);
//            rangeQueryBatch(rc,queries,ts2);
            delete r;
            delete ts1;
//            delete rc;
            delete ts2;
            delete file0;
            delete file1;
            delete file2;
            delete diskfile0;
            delete diskfile1;
            delete diskfile2;
        }
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    catch (...) {

    }
    return 0;
}
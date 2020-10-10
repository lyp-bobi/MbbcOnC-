//
// Created by Chuang on 2019/9/3.
//

//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main() {
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs();
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGTToTrajs("D://00.txt");
        int maxseg = 0;
        double para[] = {1800};
        auto stat = trajStat::instance();
        double queryLen = 3600;
        vector<IShape *> queries;
        for (int i = 0; i < 1000; i++) {
            auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
            Trajectory *concate = new Trajectory();
            double ts = ori->randomPoint().m_time - queryLen/2;
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
            maxseg = 500;
//            string name[6]={"name1","name2","name3","name4","name5","name6"};
//            id_type indexIdentifier[6];
//            IStorageManager* diskfile[6];
//            StorageManager::IBuffer* file[6];
//            for(int i=0;i<6;i++){
//                diskfile[i]= StorageManager::createNewDiskStorageManager(name[i], 4096);
//                file[i] = StorageManager::createNewRandomEvictionsBuffer(*diskfile[i], 10, false),
//            }

            std::cerr << segpara << "\n";
            std::cerr << "start splitting\n";

            for (int i = 0; i < 6; i++) {
                string name = "name";
                id_type id;
                IStorageManager *diskfile = StorageManager::createNewDiskStorageManager(name, 4096);
                StorageManager::IBuffer *file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 10, false);
                vector<pair<id_type, vector<Trajectory>>> segs;
                for (const auto &traj:trajs) {
                    std::vector<Trajectory> seg;
                    switch (i) {
                        case 0:
                            seg = traj.second.getFixedSegments(int(std::ceil(segpara / stat->tl) + 1));
                            break;
                        case 1:
                            seg = traj.second.getStaticSegments(segpara);
                            break;
                        case 2:
                            seg = traj.second.getStaticSegmentsCut(segpara);
                            break;
                        case 3:
                            seg = traj.second.getGlobalSegmentsCut(segpara);
                            break;
                        case 4:
                            seg = traj.second.getRDPSegments(segpara);
                            break;
                        case 5:
                            seg = traj.second.getHybridSegments(segpara);
                            break;
                        default:
                            break;
                    }
                    maxseg = std::max(int(seg.size()), maxseg);
                    segs.emplace_back(make_pair(traj.first, seg));
                }
                TrajStore *ts = new TrajStore(name, file, 4096, maxseg + 1);
                ts->loadSegments(segs, true, false);
//                drop_cache(1);
                segs.clear();
                ISpatialIndex *r = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts, 4096, 3, id);
                kNNQueryBatch(r, queries, ts);
                delete r;
                delete ts;
                delete file;
                delete diskfile;
            }

            std::cerr << "split finished\n";
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
        for (auto &tras:queries) {
            delete tras;
        }
        queries.clear();
    }
    catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}
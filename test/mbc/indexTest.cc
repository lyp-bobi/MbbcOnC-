//
// Created by Chuang on 2020/8/19.
//

#include "testFuncs.h"


int main() {
    auto stat = trajStat::instance();
    stat->set(0, 3.35294e+12, 691571277, 764988, 4848.29, 4.383e+06, 8.75481e-07, -90, 90, -180, 179.783, -6.21356e+10,
              2.22148e+11, 180, 359.783, 2.84283e+11, 2.93544e+06);
    knncost(1, 5, 3600, 40, true, 0.012);
    double a = biSearchMax(5, 3600, 40, false, 0.012, 1000, 1000000);
    std::cerr << a;


    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("/root/TD.csv");
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
        vector<pair<id_type, Trajectory> > trajs = loadDumpedFiledToTrajs("/root/idfb.txt");
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
        auto stat = trajStat::instance();
        int maxseg = 0;
        //double segLenParas[]={100,200,300,400,500,750,1000,1500,2000,2500,3000,3500,4000};
        double segLenParas[] = {1000};

//        double segLenParas[]={60,70,80,90,100,140,150,170,180,200};
        //double queryLenParas[]={100,200,300,400,500,750,1000,1500,2000,2500,3000,3500,4000};
        double queryLenParas[] = {1000};
        std::cerr << "\nQuery lengths are:";
        for (auto p:queryLenParas) std::cerr << p << "\t";
        std::cerr << "\n";


        for (double queryLen:queryLenParas) {
            double segLen = biSearchMax(5, 3600, 40, false, 0.012);
            maxseg = 3000;
            string name0 = "name0", name1 = "name1", name2 = "name2";
            id_type indexIdentifier0, indexIdentifier1, indexIdentifier2;
            IStorageManager *diskfile0 = StorageManager::createNewDiskStorageManager(name0, 4096),
                    *diskfile1 = StorageManager::createNewDiskStorageManager(name1, 4096),
                    *diskfile2 = StorageManager::createNewDiskStorageManager(name2, 4096);
            // Create a new storage manager with the provided base name and a 4K page size.
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

            maxseg = std::max(3000, Trajectory::cutTrajsIntoFile(trajs, segLen));
            std::cerr << "maxseg is " << maxseg << "\n";

            TrajStore *ts1 = new TrajStore(name1, diskfile1, 4096, maxseg + 1);

            ts1->loadSegments(subTrajFile, true);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier2);

            SpatialIndex::MBCRTree::MBCRTree *mbcr = static_cast<SpatialIndex::MBCRTree::MBCRTree *>(r);

            ts1->flush();
            std::cerr << "Query len:" << queryLen << "\n";
            std::cerr << "Seg len:" << segLen << "\n";

            vector<IShape *> queries;
            for (int i = 0; i < testtime; i++) {
                auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
                Trajectory *concate = new Trajectory();
                double ts = ori->randomPoint().m_time - queryLen/2;
                ori->getPartialTrajectory(ts, ts + queryLen, *concate);
                if (!concate->m_points.empty())
                    queries.emplace_back(concate);
            }
            for (int i = 0; i < 5; i++) {
                bUsingSBBD = false;
                kNNQueryBatch(r, queries, ts1);
                bUsingSBBD = true;
                kNNQueryBatch(r, queries, ts1);
            }
            for (auto &shape:queries) {
                delete (shape);
            }
            queries.clear();
//            std:cerr<<*r<<*rc<<endl;

            delete r;
            delete ts1;
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
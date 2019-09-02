//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main(){
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs("D://simp.csv");
        vector<pair<id_type, vector<Trajectory>>> segs;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        double avgSegLen=100;
//        for (double segpara = 0.1; avgSegLen>10 ; segpara/=2) {
        for (double segpara = 50; segpara<=2000 ; segpara+=50) {
            maxseg=0;
            segs.clear();
            int totallen = 0, totalseg = 0;
            for (const auto &traj:trajs) {
                totallen += traj.second.m_points.size();
                auto seg = traj.second.getSegments(segpara);
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
            ts1->loadSegments(segs);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(file2, 4096, maxseg+1);
            ts2->loadSegments(segs);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            //kNN
            segs.clear();
            segs.swap(emptyseg);
            vector<IShape *> queries;
//            double segattri[]={900,3600,1000000};
            double segattri[]={900};
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
                    kNNQueryBatch(r, queries, ts1);
                    kNNQueryBatch(rc, queries, ts2);
////                    simpli= false;
                    disttype=1;
                    kNNQueryBatch(r, queries, ts1);
                    kNNQueryBatch(rc, queries, ts2);

                    cerr << "\n";
                    for(auto &tras:queries){
                        delete tras;
                    }
                    queries.clear();
//                    for(int i=0;i<200;i++){
//                        double t = int(random(stat->mint, stat->maxt));
//                        double pLow[3] = {random(stat->minx, stat->maxx), random(stat->miny, stat->maxy), t};
//                        double pHigh[3] = {pLow[0] + random(stat->Dx/40,stat->Dx*3/40), pLow[1] + random(stat->Dy/40,stat->Dy*3/40), t};
//                        Region *rg = new Region(pLow, pHigh, 3);
//                        queries.emplace_back(rg);
//                    }
//                    rangeQueryBatch(r,queries,ts1);
//                    rangeQueryBatch(rc,queries,ts2);
//                    queries.clear();
                }
            }
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
//
// Created by Chuang on 2019/10/24.
//
//
// Created by Chuang on 2019/8/19.
//

#include "testFuncs.h"

int main(){
    try {
        calcuTime[0] = 0;
        srand((int) time(NULL));
//        vector<pair<id_type, Trajectory> > trajs = loadGTFolder();
        auto stat=trajStat::instance();
        int maxseg = 0;
        double avgSegLen=100;
//        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs();
        double segLenParas[]={20,50,80,100,200,300,500,700,900,1100,1300,1500
                              ,1700,1900,2100,2500,3000};
//        double queryLenParas[]={900,3600};
        vector<pair<id_type, Trajectory> > trajs = loadGLToTrajs();;
        double queryLenParas[]={900,3600};
        std::cerr<<"Starting bUsingSimp test\n"<<"Segmentation lengths are:";
        for(auto p:segLenParas) std::cerr<<p<<"\t";
        std::cerr<<"\nQuery lengths are:";
        for(auto p:queryLenParas) std::cerr<<p<<"\t";
        std::cerr<<"\n";
        vector<vector<IShape *>> querySet;
        for (auto queryLen:queryLenParas) {
            vector<IShape *> queries;
            for (int i = 0; i < 200; i++) {
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
//        for (double segLen=300;segLen<3000;segLen+=200) {
            vector<pair<id_type, vector<Trajectory>>> segs;
            vector<pair<id_type, vector<Trajectory>>> emptyseg;
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
            StorageManager::IBuffer *file0 = StorageManager::createNewRandomEvictionsBuffer(*diskfile0, 10, false),
                    *file1 = StorageManager::createNewRandomEvictionsBuffer(*diskfile1, 10, false),
                    *file2 = StorageManager::createNewRandomEvictionsBuffer(*diskfile2, 10, false);

            TrajStore *ts1 = new TrajStore(name1, file1, 4096, maxseg+1);
            ts1->loadSegments(segs, true);
            ISpatialIndex *r = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            //kNN
            segs.clear();
            segs.swap(emptyseg);
            std::cerr<<"Seg len:"<<segLen<<"\n";
            disttype=0;
            for(const auto &qs:querySet) {
                bUsingSimp=true;
                kNNQueryBatch(r, qs, ts1);
                bUsingSimp=false;
                kNNQueryBatch(r, qs, ts1);
            }


            delete r;
            delete ts1;
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

//
// Created by Chuang on 2019/9/2.
//

//
// Created by Chuang on 2019/9/1.
//

#include "testFuncs.h"


int main(){
    try {
        vector<string> files;
        struct dirent *ptr;
        DIR *dir;
        string PATH = fileFolder;
        dir=opendir(PATH.c_str());
        while((ptr=readdir(dir))!=NULL)
        {
            if(ptr->d_name[0] == '.')
                continue;
            //cout << ptr->d_name << endl;
            files.emplace_back(PATH+ptr->d_name);
        }
//        files.emplace_back("D://00.txt");
//        files.emplace_back("D://01.txt");
        calcuTime[0] = 0;
        srand((int) time(NULL));
        vector<pair<id_type, Trajectory> > trajs;
        for(const auto &f:files){
            auto tmptj=loadGTToTrajs(f);
            trajs.insert(trajs.end(),tmptj.begin(),tmptj.end());
        }
        vector<pair<id_type, vector<Trajectory>>> segs1,segs2;
        vector<pair<id_type, vector<Trajectory>>> emptyseg;
        int maxseg = 0;
        for (double queryLen=10;queryLen<=200;queryLen+=10) {
            maxseg=0;
            segs1.clear();
            segs2.clear();
            emptyseg.clear();
            for (auto &traj:trajs) {
                auto seg = traj.second.getFixedSegments();
                maxseg = std::max(int(seg.size()), maxseg);
                segs1.emplace_back(make_pair(traj.first, seg));
            }
            double segpara1=biSearchMax(5,queryLen,40,false,0.012);
            std::cerr<<"query len:"<<queryLen<<",partial traj len:"<<segpara1<<"\n";
            for (auto &traj:trajs) {
                auto seg = traj.second.getSegments(segpara1);
                maxseg = std::max(int(seg.size()), maxseg);
                segs2.emplace_back(make_pair(traj.first, seg));
            }

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
            ts1->loadSegments(segs1);
            ISpatialIndex *r = RTree::createAndBulkLoadNewRTreeWithTrajStore(ts1, 4096, 3, indexIdentifier1);

            TrajStore *ts2 = new TrajStore(file2, 4096, maxseg+1);
            ts2->loadSegments(segs2);
            ISpatialIndex *rc = MBCRTree::createAndBulkLoadNewMBCRTreeWithTrajStore(ts2, 4096, 3, indexIdentifier2);

            //kNN
            segs1.clear();
            segs1.swap(emptyseg);
            segs2.clear();
            segs2.swap(emptyseg);
            emptyseg.clear();
            vector<IShape *> queries;
            for (int i = 0; i < 200; i++) {
                auto ori = &trajs[(int(random(0, trajs.size()))) % trajs.size()].second;
                Trajectory *concate = new Trajectory();
                double ts = std::max(ori->m_startTime(),random(ori->m_startTime(), ori->m_endTime() - queryLen));
                ori->getPartialTrajectory(ts, ts + queryLen, *concate);
                if (!concate->m_points.empty())
                    queries.emplace_back(concate);
            }
            kNNQueryBatch(r, queries, ts1);
            kNNQueryBatch(rc, queries, ts2);
            cerr << "\n";
            for(auto &tras:queries){
                delete tras;
            }
            queries.clear();
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
    return 0;
}
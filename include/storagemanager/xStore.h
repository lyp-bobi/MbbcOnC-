//
// Created by Chuang on 2019/6/10.
//

#pragma once
#define CUTENTRY pair<pair<int,int>,xSBB>
#define CUTFUNC function<queue<CUTENTRY>(xTrajectory&)>
#define CUTFUNC_PARA function<queue<CUTENTRY>(xTrajectory&, double)>
#include <spatialindex/SpatialIndex.h>
#include <cstring>
#if !WIN32
#define filedirprefix "/root/run/"
#else
#define filedirprefix "D://run/"
#endif
using namespace SpatialIndex;
using namespace SpatialIndex::StorageManager;
using std::vector;

extern bool bSecondaryIndex;
#include "xTrajIdx.h"
namespace SpatialIndex
{
    struct xStoreEntry{
        id_type m_id;
        uint32_t m_s,m_e;
        xStoreEntry(id_type id,uint32_t s,uint32_t e)
        :m_id(id),m_s(s),m_e(e){}
        xStoreEntry():m_id(0),m_s(0),m_e(0){}
    };
    namespace StorageManager {
        class xSBBData:public IData{
        public:
            xStoreEntry m_se;
            xSBB m_b;
            id_type m_sbbid;
            bool m_hasPrev = true, m_hasNext=true;
            xSBBData()=default;
            xSBBData(const id_type sbbid,const xStoreEntry &entry,const xSBB &b, bool pv=true, bool nt=true)
                :m_sbbid(sbbid),m_se(entry),m_b(b), m_hasPrev(pv), m_hasNext(nt){}
            virtual xSBBData* clone(){
                return new xSBBData(m_sbbid,m_se,m_b,m_hasPrev, m_hasNext);
            }
            virtual id_type getIdentifier() const{
                return m_sbbid;
            }
            virtual void getShape(IShape** out) const{
            }
            virtual void getData(uint32_t& len, uint8_t** data) const{
            }
        };

        class xStore;
        class SIDX_DLL xStore:public IStorageManager{
        public:
            xStore(){};
            ~xStore();
            xStore(string myname, string file, bool subtrajs=true, bool forceNew=false);
            xStore(xStore &r);
            void loadFile(string filename);
            void flush();
            void loadByteArray(const id_type page, uint32_t& len, uint8_t** data){
                auto start = std::chrono::system_clock::now();
                m_pStorageManager->loadByteArray(page,len,data);
                auto end = std::chrono::system_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                m_IOtime+=double(duration.count()) * std::chrono::microseconds::period::num
                        / std::chrono::microseconds::period::den;
            }//for inner nodes
            void storeByteArray(id_type& page, const uint32_t len, const uint8_t* const data){
                m_pStorageManager->storeByteArray(page,len,data);
            }//for inner nodes
            void deleteByteArray(const id_type page){
                m_pStorageManager->deleteByteArray(page);
            }//should not be used in bulkload mode i guess
            void cleanStatistic(){
                m_trajIO=0;m_indexIO=0;
                m_leaf1=0;m_leaf2=0;
                m_loadedTraj=0;
                m_boundingVisited=0;
                m_IOtime=0;
            }

            void loadTraj(xTrajectory &out, const xStoreEntry &e);
            xTrajectory randomSubtraj(double len);
            xPoint randomPoint();

            json m_property;
            xTrajIdx *m_trajIdx=nullptr;
            IStorageManager* m_pStorageManager=nullptr;
            std::string m_name;
            bool m_bSubTraj=false;
            uint32_t m_pageSize;
            uint32_t m_maxTrajSegs=100;

            /*statistic*/
            uint32_t m_leaf1=0,m_leaf2=0;
            uint32_t m_trajIO=0,m_indexIO=0;
            uint32_t m_loadedTraj=0;
            uint32_t m_boundingVisited=0;
            double m_IOtime=0;
        };

        class xSBBStream:public IDataStream{
        public:
            CUTFUNC m_cutFunc;
            queue<CUTENTRY> m_buf;
            xStore *m_pstore;
            //tmp recorder for split sbbs
            id_type m_id;
//#ifdef WIN32
            std::map<id_type,xTrajEntry>::iterator m_it;
//#else
//            id_entry_map::iterator m_it;
//#endif
            id_type m_count=0;
            xSBBStream(xStore *p, CUTFUNC f);
            bool hasNext();
            uint32_t size() override;
            xSBBData* getNext() override ;
            void rewind() override;
        };
    }//namespace StorageManager
}
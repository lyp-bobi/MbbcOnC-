//
// Created by Chuang on 2019/6/10.
//
#ifndef STORAGEMANGER_XSTORE_H
#define STORAGEMANGER_XSTORE_H
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
#include "tjsql.h"
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
            virtual xStore* clone();
            virtual void loadFile(string filename);
            virtual void flush();
            virtual void loadByteArray(const id_type page, uint32_t& len, uint8_t** data);
            virtual void storeByteArray(id_type& page, const uint32_t len, const uint8_t* const data);
            virtual void deleteByteArray(const id_type page){
                m_pStorageManager->deleteByteArray(page);
            }//should not be used in bulkload mode i guess
            virtual void cleanStatistic(){
                m_trajIO=0;m_indexIO=0;
                m_leaf1=0;m_leaf2=0;
                m_loadedTraj=0;
                m_boundingVisited=0;
                m_IOtime=0;
            }

            virtual void loadTraj(xTrajectory &out, const xStoreEntry &e);
            virtual xTrajectory randomSubtraj(double len);
            virtual xPoint randomPoint();



            json m_property;
            /* for file implemntation*/
            IStorageManager* m_pStorageManager=nullptr;
            /* for db implementation */

            class PageCahce
            {
            public:
                PageCahce(uint32_t l, const uint8_t* const d) : m_pData(0), m_length(l), m_bDirty(false)
                {
                    m_pData = new uint8_t[m_length];
                    memcpy(m_pData, d, m_length);
                }

                ~PageCahce() { delete[] m_pData; }

                uint8_t* m_pData;
                uint32_t m_length;
                bool m_bDirty;
                bool m_bUsing = false;
            }; // Entry
            void drop_one_cache();

#define PAGECACHE_DEFAULT false
            bool m_bPageCache = PAGECACHE_DEFAULT;
            int m_buffersize = 1024 * 16;
            std::map<id_type, PageCahce*> m_pagecache;

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
        class SIDX_DLL xStoreDB: public xStore{
        public:
            xStoreDB();
            xStoreDB(string myname, string file, bool subtrajs=true, bool forceNew=false);
            virtual xStore* clone();
            virtual void loadFile(string filename);
            virtual void flush();
            virtual void loadByteArray(const id_type page, uint32_t& len, uint8_t** data);
            virtual void storeByteArray(id_type& page, const uint32_t len, const uint8_t* const data);
            virtual void deleteByteArray(const id_type page);

            virtual void loadTraj(xTrajectory &out, const xStoreEntry &e);
            virtual xTrajectory randomSubtraj(double len);
            virtual xPoint randomPoint();
        };
        class xSBBStream:public IDataStream{
        public:
            CUTFUNC m_cutFunc;
            queue<CUTENTRY> m_buf;
            xStore *m_pstore;
            //tmp recorder for split sbbs
            id_type m_id=0;
//#ifdef WIN32
            std::map<id_type,xTrajEntry>::iterator m_it;
//#else
//            id_entry_map::iterator m_it;
//#endif
            bool m_isdb;
            id_type m_numit=0;
            id_type m_size=-1;
            id_type m_count=0;
            xSBBStream(xStore *p, CUTFUNC f);
            bool hasNext();
            uint32_t size() override;
            xSBBData* getNext() override ;
            void rewind() override;
        };
    }//namespace StorageManager
}
#endif
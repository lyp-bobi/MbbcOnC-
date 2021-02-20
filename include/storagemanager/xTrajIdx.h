//
// Created by Chuang on 2021/2/17.
//

#ifndef SPATIALINDEX_XTRAJIDX_H
#define SPATIALINDEX_XTRAJIDX_H

#pragma once
//#include <stxxl.h>
#include <thread>
#include <mutex>
using namespace SpatialIndex;
using namespace SpatialIndex::StorageManager;
using std::vector;

extern std::mutex mtx;

namespace SpatialIndex {
    namespace StorageManager {
        class xTrajEntry{
        public:
            id_type m_page=-1;
            uint32_t m_npoint=-1;
            xTrajEntry(){};
            xTrajEntry(id_type page,uint32_t len);
            std::string toString();
            xTrajEntry(string &s);
        };
        SIDX_DLL ostream & operator<<( ostream & os,const xTrajEntry & c);
        SIDX_DLL istream & operator>>( istream & is,xTrajEntry & c);
//#ifndef WIN32
//        typedef stxxl::map<id_type,xTrajEntry, CmpIdLess,(uint64_t)16*1024,(uint64_t)128*1024> id_entry_map;
//#endif
        class xTrajIdx{
        public:
//#ifdef WIN32
            std::map<id_type,xTrajEntry> m_idx;
//#else
//            id_entry_map m_idx;
//            xTrajIdx():m_idx(4*16*1024,4*128*1024){};
//#endif
            xTrajEntry& operator[](uint64_t idx){
//                mtx.lock();
                auto res = &(m_idx[idx]);
//                mtx.unlock();
                return *res;
            }
            auto begin(){return m_idx.begin();}
            auto end(){return m_idx.end();}
            auto size(){return m_idx.size();}
            void insert(pair<id_type,xTrajEntry>e){
//                mtx.lock();
                m_idx.insert(e);
//                mtx.unlock();
            }
        };
    }
}


#endif //SPATIALINDEX_XTRAJIDX_H

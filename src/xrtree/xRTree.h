/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#pragma once


#include "storagemanager/xStore.h"

#include "Statistics.h"
#include "Node.h"
#include "PointerPoolNode.h"
#include <cmath>

extern bool bulkloadt;

//#define SKIPMAXSPEED

namespace SpatialIndex
{
	namespace xRTreeNsp
	{
        class poppq;
		class xRTree:public xRTreeQueryObject
		{
                  //class NNEntry;

		public:
			xRTree(IStorageManager&, Tools::PropertySet&);
				// String                   Value     Description
				// ----------------------------------------------
				// IndexIndentifier         VT_LONG   If specified an existing index will be openened from the supplied
				//                          storage manager with the given index id. Behaviour is unspecified
				//                          if the index id or the storage manager are incorrect.
				// Dimension                VT_ULONG  Dimensionality of the data that will be inserted.
				// IndexCapacity            VT_ULONG  The index node capacity. Default is 100.
				// LeafCapactiy             VT_ULONG  The leaf node capacity. Default is 100.
				// FillFactor               VT_DOUBLE The fill factor. Default is 70%
				// TreeVariant              VT_LONG   Can be one of Linear, Quadratic or Rstar. Default is Rstar
				// NearMinimumOverlapFactor VT_ULONG  Default is 32.
				// SplitDistributionFactor  VT_DOUBLE Default is 0.4
				// ReinsertFactor           VT_DOUBLE Default is 0.3
				// EnsureTightMBRs          VT_BOOL   Default is true
				// IndexPoolCapacity        VT_LONG   Default is 100
				// LeafPoolCapacity         VT_LONG   Default is 100
				// xMBRPoolCapacity       VT_LONG   Default is 1000
				// xPointPoolCapacity        VT_LONG   Default is 500
                // todo:DataType                 VT_LONG   Can be BoundingBox and xTrajectory


            virtual ~xRTree();

            
			virtual void getIndexProperties(Tools::PropertySet& out) const;
			virtual void getStatistics(IStatistics** out) const;

            virtual void intersectsWithQuery(const xCylinder& query, IVisitor& v);
            virtual void nearestNeighborQuery(uint32_t k, const xTrajectory& query, IVisitor& v);
            virtual void findid(id_type id);
        public:
			void initNew(Tools::PropertySet&);
			void initOld(Tools::PropertySet& ps);
			void storeHeader();
			void loadHeader();
			

			id_type writeNode(Node*);
			NodePtr readNode(id_type page);
			void deleteNode(Node*);

        public:
            bool m_bUsingMBR=false;//for TBTree
            bool m_bUsingMBC=false;//for SBBRTree
            bool m_bUsingMBL=false;//for STR-Tree

            bool m_bUsingSimp=true;
            bool m_bUsingSBBD=true;

            bool m_bStoringLinks = true;

            xStore *m_ts=nullptr;

            string m_name="";

            double m_bt = 100000;

			IStorageManager* m_pStorageManager;

			id_type m_rootID, m_headerID;

			xRTreeVariant m_treeVariant;

        private:

			double m_fillFactor;

			uint32_t m_indexCapacity;

			uint32_t m_leafCapacity;

			uint32_t m_nearMinimumOverlapFactor;
				// The R*-Tree 'p' constant, for calculating nearly minimum overlap cost.
				// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
				// for xPoints and Rectangles', Section 4.1]

			double m_splitDistributionFactor;
				// The R*-Tree 'm' constant, for calculating spliting distributions.
				// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
				// for xPoints and Rectangles', Section 4.2]

			double m_reinsertFactor;
				// The R*-Tree 'p' constant, for removing entries at reinserts.
				// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
				//  for xPoints and Rectangles', Section 4.3]

			uint32_t m_dimension=3;

			Statistics m_stats;

			bool m_bTightMBRs;

            void insertData_impl(xMBR& mbr, id_type id);
            void insertData_impl(xMBR& mbr, id_type id, uint32_t level, uint8_t* overflowTable);
            bool deleteData_impl(const xMBR& mbr, id_type id);

			Tools::PointerPool<xPoint> m_xPointPool;
			Tools::PointerPool<xMBR> m_xMBRPool;
            Tools::PointerPool<xSBB> m_xSBBPool;
			Tools::PointerPool<Node> m_indexPool;
			Tools::PointerPool<Node> m_leafPool;

#ifdef HAVE_PTHREAD_H
			pthread_mutex_t m_lock;
#endif
		public:
		    class simpleData:public IData{
            public:
                simpleData(id_type id,double dist):m_id(id),m_dist(dist){}
                id_type m_id;
                double m_dist;
                virtual simpleData* clone(){throw Tools::NotSupportedException(".");}
                virtual id_type getIdentifier() const{return m_id;}
                virtual void getShape(IShape** out) const{ *out= nullptr;}
                virtual void getData(uint32_t& len, uint8_t** data) const{len=0;*data= nullptr;}
            };



			class ValidateEntry
			{
			public:
				ValidateEntry(xMBR& r, NodePtr& pNode) : m_parentMBR(r), m_pNode(pNode) {}

				xMBR m_parentMBR;
				NodePtr m_pNode;
			}; // ValidateEntry

			friend class Node;
			friend class Leaf;
			friend class Index;
			friend class BulkLoader;

			friend std::ostream& operator<<(std::ostream& os, const xRTree& t);
		}; // xRTree

		std::ostream& operator<<(std::ostream& os, const xRTree& t);


        struct leafInfo{
            xStoreEntry m_se;
            bool m_hasPrev, m_hasNext;
            double m_ts, m_te;
        };

        class NNEntry
        {
        public:
            uint32_t m_type;
            id_type m_id;
            DISTE m_dist;
            leafInfo* m_pEntry = nullptr;

            NNEntry(id_type id, DISTE f, uint32_t type)
                    : m_id(id), m_dist(f), m_type(type) {
            }
            NNEntry(id_type id, leafInfo* e, double f, uint32_t type)
                    : m_id(id), m_dist(f), m_type(type), m_pEntry(e) {
            }
            ~NNEntry() {
                if(m_pEntry!= nullptr) delete m_pEntry;
            }
        }; // NNEntry

        class NNComparator : public INearestNeighborComparator
        {
        public:
            double getMinimumDistance(const IShape& query, const IShape& entry)
            {
//                    query.getMinimumDistance(entry);
                return query.getMinimumDistance(entry);
            }

            double getMinimumDistance(const IShape& query, const IData& data)
            {
                IShape* pS;
                data.getShape(&pS);
//                    query.getMinimumDistance(*pS);
                double ret = query.getMinimumDistance(*pS);
                delete pS;
                return ret;
            }

        }; // NNComparator

        class EntryMPQ: public MutablePriorityQueue<NNEntry*>{
        public:
            EntryMPQ()
                    :MutablePriorityQueue<NNEntry*>(
                    [](NNEntry* &a, NNEntry* &b) {
                        return std::tie(a->m_dist, b->m_type, a->m_id) < std::tie(b->m_dist, a->m_type, b->m_id);
                    }
            ){}
            void update(const MutablePriorityQueue<NNEntry*>::handle_type &handle,id_type id, DISTE minDist,int type)
            {
                m_vElements[handle]->m_id=id;
                m_vElements[handle]->m_dist=minDist;
                m_vElements[handle]->m_type=type;
//                    std::cerr<<"update"<<minDist<<" "<<type<<"\n";
                size_t index=m_vHandleIndex[handle];
                if (!lower(index))
                {
                    raise(index);
                }
            }
            void updateValue(const MutablePriorityQueue<NNEntry*>::handle_type &handle,id_type id, DISTE minDist,int type)
            {
                m_vElements[handle]->m_id=id;
                m_vElements[handle]->m_dist=minDist;
                m_vElements[handle]->m_type=type;
            }
            void updateOrder(const MutablePriorityQueue<NNEntry*>::handle_type &handle)
            {
                size_t index=m_vHandleIndex[handle];
                if (!lower(index))
                {
                    raise(index);
                }
            }
            void clean(){
                for(long i=1;i<m_vHandleHeap.size();i++){
                    id_type idd=m_vHandleHeap[i];
                    delete m_vElements[idd];
                }
                m_vHandleHeap.clear();
            }
        };


        class poppq{
        public:
            std::list<std::pair<double,id_type>> m_pq;
            int m_len=5;
            void setLen(int k){
                m_len=k;
                while(m_pq.size()>m_len) m_pq.pop_back();
            }
            double threshold(){
                if(m_pq.size()==m_len)
                    return m_pq.back().first;
                return 1e300;
            }
            bool insert(id_type id, double priority){
                bool doi= false;
                if(m_pq.size()<m_len){
                    doi = true;
                }
                else if(priority < m_pq.back().first){
                    doi = true;
                }
                for(auto it = m_pq.begin();it != m_pq.end();it++){
                    if(it->second == id){
                        m_pq.erase(it);
                        doi = true;
                        break;
                    }
                }
                if(doi){
                    for(auto it = m_pq.begin();it != m_pq.end();it++){
                        if (priority<it->first){
                            m_pq.insert(it,std::make_pair(priority,id));
                            if(m_pq.size()>m_len){
                                m_pq.pop_back();
                            }
                            return true;
                        }
                    }
                    if(m_pq.size()<m_len){
                        m_pq.emplace_back(priority,id);
                        return true;
                    }
                }
                return false;
            }

        };


        class PartsStore{ /* for SBB-Driven*/
        protected:
            class Parts{
            public:
                class slab{
                public:
                    double ts,te;
                    DISTE d;
                    xMBR ms,me;
                    slab(double s, double e):
                            ts(s),te(e),d(0,0,true){}
                    slab(double s, double e,DISTE dd):
                            ts(s),te(e),d(dd){}
                };
                PartsStore* m_ps;
                std::set<id_type> m_missingLeaf, m_loadedLeaf;
                std::vector<xSBB> m_UCsbbs;
                std::list<slab> m_line;
                std::map<std::pair<double,double>,DISTE> m_computedDist;
                std::map<double,xStoreEntry> m_ses;
                DISTE m_calcMin;
                double m_mintime=1e300,m_maxtime=-1e300;
                bool m_hasPrev=true,m_hasNext=true;
                double m_computedTime=0,m_loadedTime=0;
                bool is_modified = true;
                Parts(PartsStore* ps= nullptr){
                    m_ps = ps;
                    m_line.emplace_back(slab(ps->m_query.m_startTime(),ps->m_query.m_endTime()));
                }
                void insert(xSBB &r,id_type prev,id_type next,xStoreEntry &entry);
                void putSBB(xSBB &b);
            };

            std::map<id_type ,MutablePriorityQueue<NNEntry>::handle_type > m_handlers;
            EntryMPQ m_mpq;
            EntryMPQ m_nodespq;
            xTrajectory m_query;
            double m_error;
            xStore* m_ts;
            poppq m_pes;
            std::set<id_type> m_except;
            xRTree * m_pTree;
            trajStat* stat = trajStat::instance();
            std::map<id_type ,Parts> m_parts;
            int m_dimension=2;
            std::set<id_type > loadedLeaf;
            bool insert(id_type id, xSBB &b, id_type leafid,id_type prev,id_type next,xStoreEntry &entry);

            DISTE updateValue(id_type id,bool Inqueue=true) ;

        public:
            string explain(id_type id){
                stringstream  ss;
                auto s =m_parts[id];
                ss<<"id is "<<id<<endl;
//                for(auto &b:m_parts[id].m_UCsbbs){
//                    ss<<b.second.toString()<<endl;
//                }
                for(auto &b:m_parts[id].m_computedDist){
                    if(b.second.infer== false){
                        ss<<b.first.first<<"\t"<<b.first.second<<"\t"<<b.second.opt<<endl;
                    }
                }
                return ss.str();
            }
            bool isLoaded(id_type id){ return loadedLeaf.count(id)>0;}
            void loadLeaf(const Node &n, double dist = 0);

            NNEntry* top();

            NNEntry* pop(int type);
            void clean(){
                m_mpq.clean();
                m_nodespq.clean();
            }

            void push(NNEntry* e);

            auto empty(){return m_mpq.empty()&&m_nodespq.empty();}
            id_type getOneMissingPart(id_type id) {
                if(m_parts[id].m_missingLeaf.empty()){
                    throw Tools::IllegalStateException("wrong NNEntry Type");
                }
                return (*m_parts[id].m_missingLeaf.begin());
            }

            auto nodetop(){return m_nodespq.top();}

             PartsStore(xTrajectory &traj,double error, xRTree *r, int nnk)
                    :m_query(traj),m_error(error), m_pTree(r),m_ts(r->m_ts){
                m_pes.setLen(nnk);
            }
            ~PartsStore(){}
            xTrajectory getTraj(id_type id){
                vector<STPoint> buff;
                m_ts->m_loadedTraj+=1;
                auto ses = &m_parts[id].m_ses;
                xTrajectory res;
                m_ts->loadTraj(res,xStoreEntry(ses->begin()->second.m_id,
                                            ses->begin()->second.m_s,
                                            ses->rbegin()->second.m_e));
                return res;
            }
        };//PartStore


        class PartsStoreBFMST{
        protected:
            class Parts{
            public:
                PartsStoreBFMST* m_ps;
                std::set<id_type> m_missingLeaf, m_loadedLeaf;
                std::list<xTrajectory> m_pTrajs;
                std::map<std::pair<double,double>,DISTE> m_computedDist;
                std::map<double,xStoreEntry> m_entries;
                DISTE m_calcMin;
                double m_mintime=1e300,m_maxtime=-1e300;
                bool m_hasPrev=true,m_hasNext=true;
                double m_computedTime=0,m_loadedTime=0;
                Parts(PartsStoreBFMST* ps= nullptr):m_ps(ps){}
                void insert(xTrajectory &r,id_type prev,id_type next,xStoreEntry &entry){
                    if(m_pTrajs.empty()) m_pTrajs.emplace_back(r);
                    else{
                        auto j=m_pTrajs.begin();
                        for(;j!=m_pTrajs.end()&&(*j).m_startTime()<r.m_startTime();j++);
                        m_pTrajs.insert(j,r);
                    }
                    m_loadedTime+=r.m_endTime()-r.m_startTime();
                    if(r.m_endTime()>m_maxtime){
                        m_maxtime=r.m_endTime();
                        if(next==-1) {
                            m_hasNext=false;
                            m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
                        }
                    }
                    if(r.m_startTime()<m_mintime){
                        m_mintime=r.m_startTime();
                        if(prev==-1) {
                            m_hasPrev=false;
                            m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
                        }
                    }
                }
            };

            std::map<id_type ,MutablePriorityQueue<NNEntry>::handle_type > m_handlers;
            EntryMPQ m_mpq;
            EntryMPQ m_nodespq;
            std::set<id_type> m_except;
            poppq m_pes;
            xTrajectory m_query;
            double m_error;
            xRTree * m_pTree;
            xStore* m_ts;
            double m_lastNodeDist=0; // the top might decrease for sbbs, so we use this
            trajStat* stat = trajStat::instance();
            std::map<id_type ,Parts> m_parts;
            int m_dimension=2;
            std::set<id_type > loadedLeaf;
            void insert(id_type id, xTrajectory &r,id_type prev,id_type next,xStoreEntry &entry);
            DISTE update(id_type id) ;
        public:
            void loadPartTraj(id_type id, leafInfo * e, double dist);

            NNEntry* top();
            NNEntry* pop();
            string explain(id_type id);


            void push(NNEntry* e);
            void clean(){
                m_mpq.clean();
                m_nodespq.clean();
            }

            auto empty(){return m_mpq.empty()&&m_nodespq.empty();}
            PartsStoreBFMST(xTrajectory &traj,double error,xRTree* r, int nnk)
                    :m_query(traj),m_error(error), m_pTree(r),m_ts(r->m_ts){
                m_pes.setLen(nnk);
            }
            ~PartsStoreBFMST(){}
        };//PartStoreBFMST
	}
}


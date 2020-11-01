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

class poppq;
class PartsStore;
class PartsStoreBFMST;


namespace SpatialIndex
{
	namespace xRTreeNsp
	{
		class xRTree
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

            virtual void intersectsWithQuery(const IShape& query, IVisitor& v);
            virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v);

		private:
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

            bool m_bStoringLinks = true;

            xStore *m_ts=nullptr;

            string m_name="";

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

            NNEntry(id_type id, double f, uint32_t type)
                    : m_id(id), m_dist(f), m_type(type) {
            }
            NNEntry(id_type id, leafInfo* e, double f, uint32_t type)
                    : m_id(id), m_dist(f), m_type(type), m_pEntry(e) {
            }
            ~NNEntry() {}
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
                    if(m_vElements[idd]->m_pEntry!= nullptr){
                        delete(m_vElements[idd]->m_pEntry);
                    }
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
                        m_pq.push_back(std::make_pair(priority,id));
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
                PartsStore* m_ps;
                std::set<id_type> m_missingLeaf, m_loadedLeaf;
                std::list<xSBB> m_sbbs;
                std::map<std::pair<double,double>,DISTE> m_computedDist;
                std::map<double,xStoreEntry> m_ses;
                DISTE m_calcMin;
                double m_mintime=1e300,m_maxtime=-1e300;
                bool m_hasPrev=true,m_hasNext=true;
                double m_computedTime=0,m_loadedTime=0;
                Parts(PartsStore* ps= nullptr):m_ps(ps){}
                void insert(xSBB &r,id_type prev,id_type next,xStoreEntry &entry){
                    if(m_sbbs.empty()) m_sbbs.emplace_back(r);
                    else{
                        auto j=m_sbbs.begin();
                        for(;j!=m_sbbs.end()&&j->startTime()<r.startTime();j++);
                        m_sbbs.insert(j,r);
                    }
                    m_loadedTime+=r.endTime()-r.startTime();
                    if(r.endTime()>m_maxtime){
                        m_maxtime=r.endTime();
                        if(next==-1) {
                            m_hasNext=false;
                            m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
                        }
                    }
                    if(r.startTime()<m_mintime){
                        m_mintime=r.startTime();
                        if(prev==-1) {
                            m_hasPrev=false;
                            m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
                        }
                    }
                    if(prev>=0&&m_loadedLeaf.count(prev)==0) m_missingLeaf.insert(prev);
                    if(next>=0&&m_loadedLeaf.count(next)==0) m_missingLeaf.insert(next);
                    m_ses[r.startTime()]=entry;
                }
            };

            std::map<id_type ,MutablePriorityQueue<NNEntry>::handle_type > m_handlers;
            EntryMPQ m_mpq;
            EntryMPQ m_nodespq;
//    bool m_useMBR=false;
            xTrajectory m_query;
            double m_error;
            xStore* m_ts;
            poppq m_pes;
            trajStat* stat = trajStat::instance();
            std::map<id_type ,Parts> m_parts;
            int m_dimension=2;
            std::set<id_type > loadedLeaf;
            void insert(id_type id, xSBB &b,id_type prev,id_type next,xStoreEntry &entry){
                if(m_parts.count(id)==0){
                    m_parts[id]=Parts(this);
                }
                m_parts[id].insert(b,prev,next,entry);
            }

            DISTE updateValue(id_type id) {
                Parts *parts = &m_parts[id];
                double computedTime = 0;
                DISTE pd;
                DISTE res;
                std::pair<double, double> timeInterval;
                //inferred distance(front dist, back dist and mid dist) should be stored as negative values
                //front dist
                if (parts->m_mintime > m_query.m_startTime()) {
                    timeInterval.first=m_query.m_startTime();
                    timeInterval.second=parts->m_mintime;
                    if (parts->m_computedDist.count(timeInterval) > 0) {
                        pd= parts->m_computedDist[timeInterval];
                    } else {
                        if (parts->m_hasPrev) {
                            pd = m_query.frontDist(parts->m_sbbs.front(),stat->vmax);
                            parts->m_computedDist[timeInterval] = pd;
                        } else {
                            pd = m_query.frontDistStatic(parts->m_sbbs.front());
                            parts->m_computedDist[timeInterval] = pd;
                            computedTime += timeInterval.second - timeInterval.first;
                        }
                    }
                    if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty())
                        pd.opt=std::max(pd.opt, m_nodespq.top()->m_dist.opt *
                                                (timeInterval.second - timeInterval.first) / (m_query.m_endTime()-m_query.m_startTime()));
                    res = res + pd;
                }
                //mid dist
                const xSBB *prev= nullptr;
                for (const auto &box:parts->m_sbbs) {
                    //this box
                    timeInterval.first=box.startTime();
                    timeInterval.second=box.endTime();
                    if (parts->m_computedDist.count(timeInterval) > 0) {
                        pd = parts->m_computedDist[timeInterval];
                    } else {
                        pd = DISTE(m_query.sbbDist(box));
                        parts->m_computedDist[timeInterval] = pd;
                    }
                    res = res + pd;
                    computedTime += timeInterval.second - timeInterval.first;
                    //the gap
                    if (box.startTime() != parts->m_sbbs.front().startTime()) {//not first
                        if (prev->endTime() < box.startTime()) {
                            timeInterval.first=prev->endTime();
                            timeInterval.second=box.startTime();
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                pd= parts->m_computedDist[timeInterval];
                            } else {
                                pd = m_query.gapDist(*prev, box, stat->vmax);
                                parts->m_computedDist[timeInterval] = pd;
                            }
                            if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty())
                                pd.opt=std::max(pd.opt, m_nodespq.top()->m_dist.opt *
                                                        (timeInterval.second - timeInterval.first) / (m_query.m_endTime()-m_query.m_startTime()));
                            res = res + pd;
                        }
                    }
                    prev = &box;
                }
                //backdist
                if (parts->m_maxtime < m_query.m_endTime()) {
                    timeInterval.first=parts->m_maxtime;
                    timeInterval.second=m_query.m_endTime();
                    if (parts->m_computedDist.count(timeInterval) > 0) {
                        pd= parts->m_computedDist[timeInterval];
                    } else {
                        if (parts->m_hasNext) {
                            pd = m_query.backDist(parts->m_sbbs.back(),stat->vmax);
                            parts->m_computedDist[timeInterval] = pd;
                        } else {
                            pd = m_query.backDistStatic(parts->m_sbbs.back());
                            parts->m_computedDist[timeInterval] = pd;
                            computedTime += timeInterval.second - timeInterval.first;
                        }
                    }
                    if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty())
                        pd.opt=std::max(pd.opt, m_nodespq.top()->m_dist.opt *
                                                (timeInterval.second - timeInterval.first) / (m_query.m_endTime()-m_query.m_startTime()));
                    res = res + pd;
                }
                parts->m_calcMin = res;
                parts->m_computedTime = computedTime;
                int type = 2;
                if (parts->m_missingLeaf.empty()) type = 3;
                res.opt -= m_error;
                res.pes += m_error;
                m_mpq.updateValue(m_handlers[id], id, res, type);
                return res;
            }

        public:
            bool isLoaded(id_type id){ return loadedLeaf.count(id)>0;}
            void loadLeaf(const Node &n, double dist = 0){
//                    std::cerr<<"load leaf"<<n.m_nodeMBR<<"\n";
//                    std::cerr<<"leaf dist"<<m_query.getNodeMinimumDistance(n.m_nodeMBR,100)/(m_query.m_endTime()-m_query.m_startTime())<<"\n";
//                    std::cerr<<"load leaf"<<n.m_identifier<<"\n";
                loadedLeaf.insert(n.m_identifier);
                std::set<id_type > relatedIds;
                for(int i=0;i<n.m_children;i++){
                    id_type trajid=n.m_pIdentifier[i];
                    xStoreEntry entry= n.m_se[i];
                    double bts=n.m_ptrxSBB[i]->startTime(),bte=n.m_ptrxSBB[i]->endTime();
                    if(bts>=m_query.m_endTime()||
                       bte<=m_query.m_startTime()){}
                    else {
                        insert(trajid, *n.m_ptrxSBB[i],
                               (m_query.m_startTime()<bts)?n.m_prevNode[i]:-1
                                , (m_query.m_endTime()>bte)?n.m_nextNode[i]:-1,
                               entry);
                        relatedIds.insert(trajid);
                    }
                }
                for(const auto &rid:relatedIds){
                    m_parts[rid].m_missingLeaf.erase(n.m_identifier);
                    m_parts[rid].m_loadedLeaf.insert(n.m_identifier);
                    if(m_handlers.count(rid)==0){
                        auto handle = m_mpq.push(new NNEntry(rid, dist, 2));
                        m_handlers[rid] = handle;
                    }
                }
            }

            auto top(){
                if(!m_mpq.empty()&&m_mpq.top()->m_type==2){
                    id_type lastid=m_mpq.top()->m_id;
                    updateValue(lastid);
                    m_mpq.updateOrder(m_handlers[lastid]);
                    while(m_mpq.top()->m_type==2&&m_mpq.top()->m_id!=lastid){
                        lastid=m_mpq.top()->m_id;
                        updateValue(lastid);
                        m_mpq.updateOrder(m_handlers[lastid]);
                    }
                }
                if(!m_mpq.empty()&&(m_nodespq.empty()|| m_mpq.top()->m_dist < m_nodespq.top()->m_dist))
                    return m_mpq.top();
                else
                    return m_nodespq.top();
            }

            auto pop(){
                if(!m_mpq.empty()&&(m_nodespq.empty()|| m_mpq.top()->m_dist < m_nodespq.top()->m_dist))
                    return m_mpq.pop();
                else
                    return m_nodespq.pop();
            }
            void clean(){
                m_mpq.clean();
                m_nodespq.clean();
            }

            auto push(NNEntry* e){
                if(e->m_type==0||e->m_type==1){
                    return m_nodespq.push(e);
                }else{
                    return m_mpq.push(e);
                }
            }

            auto empty(){return m_mpq.empty()&&m_nodespq.empty();}
            id_type getOneMissingPart(id_type id) {
                if(m_parts[id].m_missingLeaf.empty()){
                    throw Tools::IllegalStateException("wrong NNEntry Type");
                }
                return (*m_parts[id].m_missingLeaf.begin());
            }
            auto getMissingPart(id_type id){return m_parts[id].m_missingLeaf;}

             PartsStore(xTrajectory &traj,double error,xStore* ts)
                    :m_query(traj),m_error(error),m_ts(ts){}
            ~PartsStore(){}
            xTrajectory getTraj(id_type id){
                vector<STPoint> buff;
                m_ts->m_loadedTraj+=1;
                auto ses = &m_parts[id].m_ses;
                xTrajectory res;
                m_ts->loadTraj(res,xStoreEntry(ses->begin()->second.m_id,
                                            ses->begin()->second.m_s,
                                            ses->end()->second.m_e));
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
            bool m_useMBR=false;
            xTrajectory m_query;
            double m_error;
            xStore* m_ts;
            trajStat* stat = trajStat::instance();
            std::map<id_type ,Parts> m_parts;
            int m_dimension=2;
            std::set<id_type > loadedLeaf;
            void insert(id_type id, xTrajectory &r,id_type prev,id_type next,xStoreEntry &entry){
                if(m_parts.count(id)==0){
                    m_parts[id]=Parts(this);
                }
//                    m_parts[id].insert(r,prev,next,entry);
                xTrajectory tmp;
                tmp.m_points.emplace_back(r.m_points[0]);
                tmp.m_points.emplace_back(r.m_points[r.m_points.size()-1]);
                m_parts[id].m_computedDist[std::make_pair(r.m_startTime(),r.m_endTime())]
                        = DISTE( r.getMinimumDistance(m_query));
                m_parts[id].insert(tmp,prev,next,entry);
//                    std::cerr<<"part"<<id<<"\t"<<m_parts[id].m_loadedTime<<"\t"<<m_query.m_endTime()-m_query.m_startTime()<<"\n";
                update(id);
            }

            DISTE update(id_type id) {
                if(m_except.count(id)>0){
                    return DISTE(1e300);
                }
                Parts *parts = &m_parts[id];
                double computedTime = 0;
                DISTE pd(0);
                DISTE res;
                std::pair<double, double> timeInterval;
                if (disttype == 0) {
                    //inferred distance(front dist, back dist and mid dist) should be stored as negative values
                    //front dist
                    if (parts->m_mintime > m_query.m_startTime()) {
                        timeInterval.first=m_query.m_startTime();
                        timeInterval.second=parts->m_mintime;
                        if (parts->m_computedDist.count(timeInterval) > 0) {
                            pd= parts->m_computedDist[timeInterval];
                        } else {
                            if (parts->m_hasPrev) {
                                pd = m_query.frontDist(parts->m_pTrajs.front().m_points.front(),
                                                         stat->vmax);
                                parts->m_computedDist[timeInterval] = pd;
                            } else {
                                pd = m_query.frontDistStatic(parts->m_pTrajs.front().m_points.front());
                                parts->m_computedDist[timeInterval] = pd;
                                computedTime += timeInterval.second - timeInterval.first;
                            }
                        }
                        if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty())
                            pd.opt=std::max(pd.opt, m_nodespq.top()->m_dist.opt *
                                                    (timeInterval.second - timeInterval.first) / (m_query.m_endTime()-m_query.m_startTime()));
                        res = res+pd;
                    }
                    //mid dist
                    const xTrajectory *prev= nullptr;
                    for (const auto &traj:parts->m_pTrajs) {
                        //this box
                        timeInterval.first=traj.m_startTime();
                        timeInterval.second=traj.m_endTime();
                        if (parts->m_computedDist.count(timeInterval) > 0) {
                            if (!parts->m_computedDist[timeInterval].infer) {
                                pd = parts->m_computedDist[timeInterval];
                            } else {
                                pd = DISTE(traj.getMinimumDistance(m_query));
                                parts->m_computedDist[timeInterval] = pd;
                            }
                        } else {
                            pd = DISTE(traj.getMinimumDistance(m_query));
                            parts->m_computedDist[timeInterval] = pd;
                        }
                        res = res+pd;
                        computedTime += timeInterval.second - timeInterval.first;
                        //the gap
                        if (traj.m_startTime() != parts->m_pTrajs.front().m_startTime()) {//not first
                            if (prev->m_endTime() < traj.m_startTime()) {
                                timeInterval.first=prev->m_endTime();
                                timeInterval.second=traj.m_startTime();
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd= parts->m_computedDist[timeInterval];
                                } else {
                                    pd = m_query.gapDist(prev->m_points.back(), traj.m_points.front(), stat->vmax);
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty())
                                    pd.opt=std::max(pd.opt, m_nodespq.top()->m_dist.opt *
                                                            (timeInterval.second - timeInterval.first) / (m_query.m_endTime()-m_query.m_startTime()));
                                res = res+pd;
                            }
                        }
                        prev = &traj;
                    }
                    //backdist
                    if (parts->m_maxtime < m_query.m_endTime()) {
                        timeInterval.first=parts->m_maxtime;
                        timeInterval.second=m_query.m_endTime();
                        if (parts->m_computedDist.count(timeInterval) > 0) {
                            pd= parts->m_computedDist[timeInterval];
                        } else {
                            if (parts->m_hasNext) {
                                pd = m_query.backDist(parts->m_pTrajs.back().m_points.back(),stat->vmax);
                                parts->m_computedDist[timeInterval] = pd;
                            } else {
                                pd = m_query.backDistStatic(parts->m_pTrajs.back().m_points.back());
                                parts->m_computedDist[timeInterval] = pd;
                                computedTime += timeInterval.second - timeInterval.first;
                            }
                        }
                        if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                            pd.opt=std::max(pd.opt, m_nodespq.top()->m_dist.opt *
                                                    (timeInterval.second - timeInterval.first) / (m_query.m_endTime()-m_query.m_startTime()));
                        res = res+pd;
                    }

                    parts->m_calcMin = res;
                    parts->m_computedTime = computedTime;
                    int type = 2;
                    if (parts->m_loadedTime + 1e-7 >= (m_query.m_endTime() - m_query.m_startTime())) type = 3;
                    res.opt -= m_error;
                    res.pes += m_error;
                    if (m_handlers.count(id) == 0) {
                        auto handle = m_mpq.push(new NNEntry(id,nullptr, res, type));
                        m_handlers[id] = handle;
                    } else {
                        m_mpq.update(m_handlers[id], id, res, type);
                    }
                    m_pes.insert(id, res.pes);
                    if(res.opt>m_pes.threshold()){
                        m_except.insert(id);
                    }
                    return res;
                }
                else if(disttype==1){
                    throw Tools::IllegalStateException("maxsed nolonger supported.");
                }
                else throw Tools::IllegalStateException("");
            }
        public:
            void loadPartTraj(id_type id, leafInfo * e, double dist){
                xTrajectory tmpTraj;
                id_type trajid = e->m_se.m_id;
                m_ts->loadTraj(tmpTraj,e->m_se);
                if(m_except.count(trajid)>0) return;
                xTrajectory inter;
                tmpTraj.getPartialxTrajectory(max(m_query.m_startTime(),e->m_ts),
                                              min(m_query.m_endTime(),e->m_te),inter);
                insert(trajid, inter, e->m_hasPrev ? 1 : -1, e->m_hasNext ? 1 : -1, e->m_se);
            }

            auto top(){
                if(!m_mpq.empty()){
                    id_type prevtop=-1;
                    while(m_mpq.top()->m_id!=prevtop){
                        prevtop=m_mpq.top()->m_id;
                        update(prevtop);
                    }
                }
                if(!m_mpq.empty()&&m_mpq.top()->m_type==3&&(m_nodespq.empty()|| m_mpq.top()->m_dist < m_nodespq.top()->m_dist)){
                    return m_mpq.top();
                }
                if(!m_nodespq.empty())
                    return m_nodespq.top();
                return m_mpq.top();
            }

            auto pop(){
                if(!m_mpq.empty()&&m_mpq.top()->m_type==3&&(m_nodespq.empty()|| m_mpq.top()->m_dist < m_nodespq.top()->m_dist))
                    return m_mpq.pop();
                if(!m_nodespq.empty())
                    return m_nodespq.pop();
                return m_mpq.pop();
            }

            auto push(NNEntry* e){
                if(e->m_type==0||e->m_type==1){
                    return m_nodespq.push(e);
                }else{
                    return m_mpq.push(e);
                }
            }
            void clean(){
                m_mpq.clean();
                m_nodespq.clean();
            }

            auto empty(){return m_mpq.empty()&&m_nodespq.empty();}
            PartsStoreBFMST(xTrajectory &traj,double error,xStore* ts,bool useMBR)
                    :m_query(traj),m_error(error),m_useMBR(useMBR),m_ts(ts){}
            ~PartsStoreBFMST(){}
        };//PartStoreBFMST
	}
}


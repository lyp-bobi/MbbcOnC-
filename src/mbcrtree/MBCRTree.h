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

#include "Statistics.h"
#include "Node.h"
#include "PointerPoolNode.h"
#include "storagemanager/TrajStore.h"
#include <cmath>
#include <thread>
#include <mutex>

extern bool bUsingSimp;
extern bool bUsingSBBD;

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


namespace SpatialIndex
{
	namespace MBCRTree
	{
		class MBCRTree : public ISpatialIndex
		{
                  //class NNEntry;

		public:
			MBCRTree(IStorageManager&, Tools::PropertySet&);
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
				// RegionPoolCapacity       VT_LONG   Default is 1000
				// PointPoolCapacity        VT_LONG   Default is 500
                // todo:DataType                 VT_LONG   Can be BoundingBox and Trajectory


            virtual ~MBCRTree();



			//
			// ISpatialIndex interface
			//
			virtual void insertData(uint32_t len, const uint8_t* pData, const IShape& shape, id_type shapeIdentifier);
			virtual bool deleteData(const IShape& shape, id_type id);
			virtual void containsWhatQuery(const IShape& query, IVisitor& v);
			virtual void intersectsWithQuery(const IShape& query, IVisitor& v);
			virtual void pointLocationQuery(const Point& query, IVisitor& v);
			virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v, INearestNeighborComparator&);
			virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v);
			virtual void selfJoinQuery(const IShape& s, IVisitor& v);
			virtual void queryStrategy(IQueryStrategy& qs);
			virtual void getIndexProperties(Tools::PropertySet& out) const;
			virtual void addCommand(ICommand* pCommand, CommandType ct);
			virtual bool isIndexValid();
			virtual void getStatistics(IStatistics** out) const;

		private:
			void initNew(Tools::PropertySet&);
			void initOld(Tools::PropertySet& ps);
			void storeHeader();
			void loadHeader();

			void insertData_impl(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id);
			void insertData_impl(uint32_t dataLength, uint8_t* pData, Region& mbr, id_type id, uint32_t level, uint8_t* overflowTable);
			bool deleteData_impl(const Region& mbr, id_type id);

			id_type writeNode(Node*);
			NodePtr readNode(id_type page);
			void deleteNode(Node*);

			void rangeQuery(RangeQueryType type, const IShape& query, IVisitor& v);
			void selfJoinQuery(id_type id1, id_type id2, const Region& r, IVisitor& vis);
            void visitSubTree(NodePtr subTree, IVisitor& v);
        public:
            bool m_bUsingMBR=false;
            bool m_bStoringLinks = true;

            TrajStore *m_ts=nullptr;

        private:
			IStorageManager* m_pStorageManager;

			id_type m_rootID, m_headerID;

			MBCRTreeVariant m_treeVariant;

        private:

			double m_fillFactor;

			uint32_t m_indexCapacity;

			uint32_t m_leafCapacity;

			uint32_t m_nearMinimumOverlapFactor;
				// The R*-Tree 'p' constant, for calculating nearly minimum overlap cost.
				// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
				// for Points and Rectangles', Section 4.1]

			double m_splitDistributionFactor;
				// The R*-Tree 'm' constant, for calculating spliting distributions.
				// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
				// for Points and Rectangles', Section 4.2]

			double m_reinsertFactor;
				// The R*-Tree 'p' constant, for removing entries at reinserts.
				// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
				//  for Points and Rectangles', Section 4.3]

			uint32_t m_dimension;

			Region m_infiniteRegion;

			Statistics m_stats;

			bool m_bTightMBRs;

			Tools::PointerPool<Point> m_pointPool;
			Tools::PointerPool<Region> m_regionPool;
            Tools::PointerPool<MBC> m_mbcPool;
			Tools::PointerPool<Node> m_indexPool;
			Tools::PointerPool<Node> m_leafPool;

			std::vector<Tools::SmartPointer<ICommand> > m_writeNodeCommands;
			std::vector<Tools::SmartPointer<ICommand> > m_readNodeCommands;
			std::vector<Tools::SmartPointer<ICommand> > m_deleteNodeCommands;

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
            struct leafInfo{
                id_type m_page;
                uint32_t m_off;
                uint32_t m_len;
                bool m_hasPrev, m_hasNext;
                double m_ts, m_te;
            };
			class NNEntry
			{
			public:
                uint32_t m_type;
				id_type m_id;
				double m_minDist;
                leafInfo* m_pEntry = nullptr;

				NNEntry(id_type id, double f, uint32_t type)
				    : m_id(id), m_minDist(f),m_type(type) {
				}
                NNEntry(id_type id, leafInfo* e, double f, uint32_t type)
                        : m_id(id), m_minDist(f),m_type(type), m_pEntry(e) {
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

			class ValidateEntry
			{
			public:
				ValidateEntry(Region& r, NodePtr& pNode) : m_parentMBR(r), m_pNode(pNode) {}

				Region m_parentMBR;
				NodePtr m_pNode;
			}; // ValidateEntry


			class EntryMPQ: public MutablePriorityQueue<NNEntry*>{
            public:
			    EntryMPQ()
			        :MutablePriorityQueue<NNEntry*>(
			                [](NNEntry* &a, NNEntry* &b) {
			                    return std::tie(a->m_minDist,b->m_type, a->m_id) < std::tie(b->m_minDist, a->m_type, b->m_id);
			                }
			                ){}
                void update(const MutablePriorityQueue<NNEntry*>::handle_type &handle,id_type id, double minDist,int type)
                {
                    m_vElements[handle]->m_id=id;
                    m_vElements[handle]->m_minDist=minDist;
                    m_vElements[handle]->m_type=type;
//                    std::cerr<<"update"<<minDist<<" "<<type<<"\n";
                    size_t index=m_vHandleIndex[handle];
                    if (!lower(index))
                    {
                        raise(index);
                    }
                }
                void updateValue(const MutablePriorityQueue<NNEntry*>::handle_type &handle,id_type id, double minDist,int type)
                {
                    m_vElements[handle]->m_id=id;
                    m_vElements[handle]->m_minDist=minDist;
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
            struct storeEntry{
                id_type m_page;
                uint32_t m_off;
                uint32_t m_len;
                bool operator <( const storeEntry& y) const{
                    return std::tie(m_page, m_off, m_len) < std::tie(y.m_page, y.m_off, y.m_len);
                }
                bool operator==(const storeEntry& y) const{
                    return std::tie(m_page, m_off, m_len) == std::tie(y.m_page, y.m_off, y.m_len);
                }
            };
			class PartsStore{ /* for SBB-Driven*/
            protected:
                class Parts{
                public:
                    PartsStore* m_ps;
                    std::set<id_type> m_missingLeaf, m_loadedLeaf;
                    std::list<Region> m_mbrs;
                    std::list<MBC> m_mbcs;
                    std::map<std::pair<double,double>,DISTE> m_computedDist;
                    std::map<double,storeEntry> m_entries;
                    double m_calcMin=0;
                    double m_mintime=1e300,m_maxtime=-1e300;
                    bool m_hasPrev=true,m_hasNext=true;
                    double m_computedTime=0,m_loadedTime=0;
                    Parts(PartsStore* ps= nullptr):m_ps(ps){}
                    void insert(Region &r,id_type prev,id_type next,storeEntry &entry){
                        if(m_mbrs.empty()) m_mbrs.emplace_back(r);
                        else{
                            auto j=m_mbrs.begin();
                            for(;j!=m_mbrs.end()&&(*j).m_pLow[m_ps->m_dimension]<r.m_pLow[m_ps->m_dimension];j++);
                            m_mbrs.insert(j,r);
                        }
                        m_loadedTime+=r.m_pHigh[m_ps->m_dimension]-r.m_pLow[m_ps->m_dimension];
                        if(r.m_pHigh[m_ps->m_dimension]>m_maxtime){
                            m_maxtime=r.m_pHigh[m_ps->m_dimension];
                            if(next==-1) {
                                m_hasNext=false;
                                m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
                            }
                        }
                        if(r.m_pLow[m_ps->m_dimension]<m_mintime){
                            m_mintime=r.m_pLow[m_ps->m_dimension];
                            if(prev==-1) {
                                m_hasPrev=false;
                                m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
                            }
                        }
                        if(prev>=0&&m_loadedLeaf.count(prev)==0) m_missingLeaf.insert(prev);
                        if(next>=0&&m_loadedLeaf.count(next)==0) m_missingLeaf.insert(next);
                        m_entries[r.m_pLow[m_ps->m_dimension]]=entry;
                    }
                    void insert(MBC &r,id_type prev,id_type next,storeEntry &entry){
                        if(m_mbcs.empty()) m_mbcs.emplace_back(r);
                        else{
                            auto j=m_mbcs.begin();
                            for(;j!=m_mbcs.end()&&(*j).m_startTime<r.m_startTime;j++);
                            m_mbcs.insert(j,r);
                        }
                        m_loadedTime+=r.m_endTime-r.m_startTime;
                        if(r.m_endTime>m_maxtime){
                            m_maxtime=r.m_endTime;
                            if(next==-1) {
                                m_hasNext=false;
                                m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
                            }
                        }
                        if(r.m_startTime<m_mintime){
                            m_mintime=r.m_startTime;
                            if(prev==-1) {
                                m_hasPrev=false;
                                m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
                            }
                        }
                        if(prev>=0&&m_loadedLeaf.count(prev)==0) m_missingLeaf.insert(prev);
                        if(next>=0&&m_loadedLeaf.count(next)==0) m_missingLeaf.insert(next);
                        m_entries[r.m_startTime]=entry;
                    }
                    void removeCalculated(double ts,double te){}
                };

                std::map<id_type ,MutablePriorityQueue<NNEntry>::handle_type > m_handlers;
                EntryMPQ m_mpq;
                EntryMPQ m_nodespq;
                bool m_useMBR=false;
                Trajectory m_query;
                double m_error;
                TrajStore* m_ts;
                poppq m_pes;
                trajStat* stat = trajStat::instance();
                std::map<id_type ,Parts> m_parts;
                int m_dimension=2;
                std::set<id_type > loadedLeaf;
                void insert(id_type id, Region &br,id_type prev,id_type next,storeEntry &entry){
                    if(m_parts.count(id)==0){
                        m_parts[id]=Parts(this);
                    }
                    m_parts[id].insert(br,prev,next,entry);
                }
                void insert(id_type id, MBC &bc,id_type prev,id_type next,storeEntry &entry){
                    if(m_parts.count(id)==0){
                        m_parts[id]=Parts(this);
                    }
                    m_parts[id].insert(bc,prev,next,entry);
                }

                double updateValue(id_type id) {
                    Parts *parts = &m_parts[id];
                    double computedTime = 0;
                    DISTE pd;
                    double sum = 0,pessi = 0;
                    std::pair<double, double> timeInterval;
                    if (disttype == 0) {
                        if (m_useMBR) {
                            //front dist
                            if (parts->m_mintime > m_query.m_startTime()) {
                                timeInterval.first=m_query.m_startTime();
                                timeInterval.second=parts->m_mintime;
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd= parts->m_computedDist[timeInterval];
                                } else {
                                    if (parts->m_hasPrev) {
                                        pd = m_query.getFrontIED(parts->m_mbrs.front(), stat->vmax);
                                        parts->m_computedDist[timeInterval] = pd;
                                    } else {
                                        pd = DISTE(m_query.getStaticIED(parts->m_mbrs.front(), m_query.m_startTime(),
                                                                  parts->m_mintime));
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                }
                                if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                    pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                   (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                sum += pd.opt;
                                pessi +=pd.pes;
                            }
                            //mid dist
                            const Region *prev= nullptr;
                            for (const auto &box:parts->m_mbrs) {
                                //this box
                                timeInterval.first=box.m_pLow[m_dimension];
                                timeInterval.second=box.m_pHigh[m_dimension];
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd = parts->m_computedDist[timeInterval];
                                } else {
                                    pd = DISTE(m_query.getMinimumDistance(box));
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                sum += pd.opt;
                                pessi +=pd.pes;
                                computedTime += timeInterval.second - timeInterval.first;
                                //the gap
                                if (box.m_pLow[m_dimension] != parts->m_mbrs.front().m_pLow[m_dimension]) {
                                    if (prev->m_pHigh[m_dimension] < box.m_pLow[m_dimension]) {
                                        timeInterval.first=prev->m_pHigh[m_dimension];
                                        timeInterval.second=box.m_pLow[m_dimension];
                                        if (parts->m_computedDist.count(timeInterval) > 0) {
                                            pd= parts->m_computedDist[timeInterval];
                                        } else {
                                            pd = m_query.getMidIED(*prev, box, stat->vmax);
                                            parts->m_computedDist[timeInterval] = pd;
                                        }
                                        if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                            pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                           (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                        sum+=pd.opt;
                                        pessi+=pd.pes;
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
                                        pd = m_query.getBackIED(parts->m_mbrs.back(), stat->vmax);
                                        parts->m_computedDist[timeInterval] = pd;
                                    } else {
                                        pd = DISTE(m_query.getStaticIED(parts->m_mbrs.back(), parts->m_maxtime,
                                                                  m_query.m_endTime()));
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                }
                                if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                    pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                   (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                sum+=pd.opt;
                                pessi+=pd.pes;
                            }
                        }
                        else {
                            //inferred distance(front dist, back dist and mid dist) should be stored as negative values
                            //front dist
                            if (parts->m_mintime > m_query.m_startTime()) {
                                timeInterval.first=m_query.m_startTime();
                                timeInterval.second=parts->m_mintime;
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd= parts->m_computedDist[timeInterval];
                                } else {
                                    if (parts->m_hasPrev) {
                                        pd = m_query.getFrontIED(parts->m_mbcs.front().m_pLow[0],
                                                                 parts->m_mbcs.front().m_pLow[1],
                                                                 parts->m_mbcs.front().m_startTime,
                                                                 stat->vmax);
                                        parts->m_computedDist[timeInterval] = pd;
                                    } else {
                                        pd = DISTE(m_query.getStaticIED(parts->m_mbcs.front().m_pLow[0],
                                                                  parts->m_mbcs.front().m_pLow[1],
                                                                  m_query.m_startTime(), parts->m_mintime));
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                }
                                if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                    pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                   (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                sum+=pd.opt;
                                pessi += pd.pes;
                            }
                            //mid dist
                            const MBC *prev= nullptr;
                            for (const auto &box:parts->m_mbcs) {
                                //this box
                                timeInterval.first=box.m_startTime;
                                timeInterval.second=box.m_endTime;
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd = parts->m_computedDist[timeInterval];
                                } else {
                                    pd = DISTE(m_query.getMinimumDistance(box));
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                sum += pd.opt;
                                pessi += pd.pes;
                                computedTime += timeInterval.second - timeInterval.first;
                                //the gap
                                if (box.m_startTime != parts->m_mbcs.front().m_startTime) {//not first
                                    if (prev->m_endTime < box.m_startTime) {
                                        timeInterval.first=prev->m_endTime;
                                        timeInterval.second=box.m_startTime;
                                        if (parts->m_computedDist.count(timeInterval) > 0) {
                                            pd= parts->m_computedDist[timeInterval];
                                        } else {
                                            pd = m_query.getMidIED(*prev, box, stat->vmax);
                                            parts->m_computedDist[timeInterval] = pd;
                                        }
                                        if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                            pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                           (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                        sum+=pd.opt;
                                        pessi+=pd.pes;
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
                                        pd = m_query.getBackIED(parts->m_mbcs.back().m_pHigh[0],
                                                                 parts->m_mbcs.back().m_pHigh[1],
                                                                 parts->m_mbcs.back().m_endTime, stat->vmax);
                                        parts->m_computedDist[timeInterval] = pd;
                                    } else {
                                        pd = DISTE(m_query.getStaticIED(parts->m_mbcs.back().m_pHigh[0],
                                                                  parts->m_mbcs.back().m_pHigh[1], parts->m_maxtime,
                                                                  m_query.m_endTime()));
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                }
                                if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                    pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                   (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                sum+=pd.opt;
                                pessi+=pd.pes;
                            }
                        }

                        parts->m_calcMin = sum;
                        parts->m_computedTime = computedTime;
                        int type = 2;
                        if (parts->m_missingLeaf.empty()) type = 3;
                        sum = std::max(0.0, sum - m_error);
                        m_mpq.updateValue(m_handlers[id], id, sum, type);
                        return sum;
                    }
                    else if(disttype==1){
                        Tools::IllegalStateException("MaxSED is no longer supported now");
                    }
                    else throw Tools::IllegalStateException("");
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
                        id_type trajid=m_ts->getTrajId(n.m_pIdentifier[i]);
                        storeEntry entry= {
                                .m_page=n.m_pageNum[i],
                                .m_off=n.m_pageOff[i],
                                .m_len=n.m_dataLen[i]
                        };
                        if(m_useMBR) {
                            double bts=n.m_ptrMBR[i]->m_pLow[m_dimension],bte=n.m_ptrMBR[i]->m_pHigh[m_dimension];
                            if(bts>=m_query.m_endTime()||
                                bte<=m_query.m_startTime()){}
                            else{
                                insert(trajid, *n.m_ptrMBR[i],
                                       (m_query.m_startTime()<bts)?n.m_prevNode[i]:-1
                                        , (m_query.m_endTime()>bte)?n.m_nextNode[i]:-1,
                                        entry);
                                relatedIds.insert(trajid);
                            }

                        }else{
                            double bts=n.m_ptrMBC[i]->m_startTime,bte=n.m_ptrMBC[i]->m_endTime;
                            if(bts>=m_query.m_endTime()||
                               bte<=m_query.m_startTime()){}
                            else {
                                insert(trajid, *n.m_ptrMBC[i],
                                       (m_query.m_startTime()<bts)?n.m_prevNode[i]:-1
                                        , (m_query.m_endTime()>bte)?n.m_nextNode[i]:-1,
                                        entry);
                                relatedIds.insert(trajid);
                            }
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
//                    std::cerr<<"rid is: ";
//                    for(const auto &rid:relatedIds){
//                        std::cerr<<rid<<"\t";
//                    }
//                    std::cerr<<"\n";

//                    auto it=relatedIds.begin();
//                        for(int i=0;i<relatedIds.size();i++){
//                            updateValue(*it);
//                            it++;
//                        }
//                    for(const auto &rid:relatedIds){
////                        updateValue(rid);
//                        m_mpq.updateOrder(m_handlers[rid]);
//                    }
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
                    if(!m_mpq.empty()&&(m_nodespq.empty()||m_mpq.top()->m_minDist<m_nodespq.top()->m_minDist))
                        return m_mpq.top();
                    else
                        return m_nodespq.top();
                }

                auto pop(){
                    if(!m_mpq.empty()&&(m_nodespq.empty()||m_mpq.top()->m_minDist<m_nodespq.top()->m_minDist))
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

                PartsStore(Trajectory &traj,double error,TrajStore* ts,bool useMBR)
                :m_query(traj),m_error(error),m_useMBR(useMBR),m_ts(ts){}
                ~PartsStore(){}
                Trajectory getTraj(id_type id){
                    vector<STPoint> buff;
                    m_ts->m_loadedTraj+=1;
//                    std::cerr<<"start getting traj\n";
//                    for(auto &bc:m_parts[id].m_mbcs){
//                        std::cerr<<*bc<<"\n";
//                    }
                    std::vector<storeEntry> toload;
                    storeEntry last = m_parts[id].m_entries.begin()->second;
                    for(const auto &pair:m_parts[id].m_entries){
                        auto e=pair.second;
                        if(e.m_page == last.m_page && e.m_off == last.m_off){
                            //pass
                        }else{//todo: replace these magic values
                            if((e.m_page - last.m_page)* 4096 + (e.m_off - last.m_off) ==
                               last.m_len){
                                /*merge two entries*/
                                last.m_len += e.m_len;
                            }
                            else{
                                toload.emplace_back(last);
                                last = e;
                            }
                        }
                    }
                    toload.emplace_back(last);

                    bool fhead, fback;
                    for(const auto &e:toload){
                        m_ts->m_trajIO += std::ceil(e.m_len / 4096.0);
                        uint32_t len = e.m_off + e.m_len;
                        uint8_t *load;
                        m_ts->loadByteArray(e.m_page, len, &load);
                        uint8_t *ptr = load + e.m_off;
                        memcpy(&fhead, ptr, sizeof(bool));
                        ptr += sizeof(bool);
                        memcpy(&fback, ptr, sizeof(bool));
                        ptr += sizeof(bool);
                        if(!buff.empty()&& fhead){
                            buff.pop_back();
                        }
                        unsigned long size;
                        memcpy(&size, ptr, sizeof(unsigned long));
                        ptr += sizeof(unsigned long);
                        STPoint  p;
                        for(int i=0;i<size;i++){
                            p.makeDimension(m_dimension);
                            memcpy(&p.m_time, ptr, sizeof(double));
                            ptr += sizeof(double);
                            memcpy(p.m_pCoords, ptr, m_dimension * sizeof(double));
                            if(i!=size-1){
                                ptr += m_dimension * sizeof(double);
                            }
                            buff.emplace_back(p);
                        }
                        delete[] load;
                    }
                    return Trajectory(buff);
                }
			};//PartStore


            class PartsStoreBFMST{
            protected:
                class Parts{
                public:
                    PartsStoreBFMST* m_ps;
                    std::set<id_type> m_missingLeaf, m_loadedLeaf;
                    std::list<Trajectory> m_pTrajs;
                    std::map<std::pair<double,double>,DISTE> m_computedDist;
                    std::map<double,storeEntry> m_entries;
                    double m_calcMin=0;
                    double m_mintime=1e300,m_maxtime=-1e300;
                    bool m_hasPrev=true,m_hasNext=true;
                    double m_computedTime=0,m_loadedTime=0;
                    Parts(PartsStoreBFMST* ps= nullptr):m_ps(ps){}
                    void insert(Trajectory &r,id_type prev,id_type next,storeEntry &entry){
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
                Trajectory m_query;
                double m_error;
                TrajStore* m_ts;
                trajStat* stat = trajStat::instance();
                std::map<id_type ,Parts> m_parts;
                int m_dimension=2;
                std::set<id_type > loadedLeaf;
                void insert(id_type id, Trajectory &r,id_type prev,id_type next,storeEntry &entry){
                    if(m_parts.count(id)==0){
                        m_parts[id]=Parts(this);
                    }
//                    m_parts[id].insert(r,prev,next,entry);
                    Trajectory tmp;
                    tmp.m_points.emplace_back(r.m_points[0]);
                    tmp.m_points.emplace_back(r.m_points[r.m_points.size()-1]);
                    m_parts[id].m_computedDist[std::make_pair(r.m_startTime(),r.m_endTime())]
                        = DISTE( r.getMinimumDistance(m_query));
                    m_parts[id].insert(tmp,prev,next,entry);
//                    std::cerr<<"part"<<id<<"\t"<<m_parts[id].m_loadedTime<<"\t"<<m_query.m_endTime()-m_query.m_startTime()<<"\n";
                    update(id);
                }

                double update(id_type id) {
                    if(m_except.count(id)>0){
                        return 1e300;
                    }
                    Parts *parts = &m_parts[id];
                    double computedTime = 0;
                    DISTE pd(0);
                    double sum = 0;
                    double pessi = 0;
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
                                    pd = m_query.getFrontIED(parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                             parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                             parts->m_pTrajs.front().m_startTime(),
                                                             stat->vmax);
                                    parts->m_computedDist[timeInterval] = pd;
                                } else {
                                    pd = DISTE(m_query.getStaticIED(parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                              parts->m_pTrajs.front().m_points.front().m_pCoords[1],
                                                              m_query.m_startTime(), parts->m_mintime));
                                    parts->m_computedDist[timeInterval] = pd;
                                    computedTime += timeInterval.second - timeInterval.first;
                                }
                            }
                            if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                               (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                            sum+=pd.opt;
                            pessi+=pd.pes;
                        }
                        //mid dist
                        const Trajectory *prev= nullptr;
                        for (const auto &traj:parts->m_pTrajs) {
                            //this box
                            timeInterval.first=traj.m_startTime();
                            timeInterval.second=traj.m_endTime();
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                if (parts->m_computedDist[timeInterval].infer==false) {
                                    pd = parts->m_computedDist[timeInterval];
                                } else {
                                    pd = DISTE(traj.getMinimumDistance(m_query));
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                            } else {
                                pd = DISTE(traj.getMinimumDistance(m_query));
                                parts->m_computedDist[timeInterval] = pd;
                            }
                            sum += pd.opt;
                            pessi +=pd.pes;
                            computedTime += timeInterval.second - timeInterval.first;
                            //the gap
                            if (traj.m_startTime() != parts->m_pTrajs.front().m_startTime()) {//not first
                                if (prev->m_endTime() < traj.m_startTime()) {
                                    timeInterval.first=prev->m_endTime();
                                    timeInterval.second=traj.m_startTime();
                                    if (parts->m_computedDist.count(timeInterval) > 0) {
                                        pd= parts->m_computedDist[timeInterval];
                                    } else {
                                        pd = m_query.getMidIED(prev->m_points.back(), traj.m_points.front(), stat->vmax);
                                        parts->m_computedDist[timeInterval] = pd;
                                    }
                                    if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                        pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                                       (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                    sum+=pd.opt;
                                    pessi += pd.pes;
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
                                    pd = m_query.getBackIED(parts->m_pTrajs.back().m_points.back().m_pCoords[0],
                                                             parts->m_pTrajs.back().m_points.back().m_pCoords[1],
                                                             parts->m_pTrajs.back().m_endTime(), stat->vmax);
                                    parts->m_computedDist[timeInterval] = pd;
                                } else {
                                    pd = DISTE(m_query.getStaticIED(parts->m_pTrajs.back().m_points.back().m_pCoords[0],
                                                              parts->m_pTrajs.back().m_points.back().m_pCoords[1], parts->m_maxtime,
                                                              m_query.m_endTime()));
                                    parts->m_computedDist[timeInterval] = pd;
                                    computedTime += timeInterval.second - timeInterval.first;
                                }
                            }
                            if(parts->m_computedDist[timeInterval].infer==true&&!m_nodespq.empty())
                                pd.opt=std::max(pd.opt,m_nodespq.top()->m_minDist*
                                               (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                            sum+=pd.opt;
                            pessi+=pd.pes;
                        }

                        parts->m_calcMin = sum;
                        parts->m_computedTime = computedTime;
                        int type = 2;
                        if (parts->m_loadedTime + 1e-7 >= (m_query.m_endTime() - m_query.m_startTime())) type = 3;
                        sum = std::max(0.0, sum - m_error);
                        if (m_handlers.count(id) == 0) {
                            auto handle = m_mpq.push(new NNEntry(id,nullptr, sum, type));
                            m_handlers[id] = handle;
                        } else {
                            m_mpq.update(m_handlers[id], id, sum, type);
                        }
                        m_pes.insert(id, pessi);
                        if(sum>m_pes.threshold()){
                            m_except.insert(id);
                        }
                        return sum;
                    }
                    else if(disttype==1){
                        throw Tools::IllegalStateException("maxsed nolonger supported.");
                    }
                    else throw Tools::IllegalStateException("");
                }
            public:
                void loadPartTraj(id_type id, leafInfo * e, double dist){
                    Trajectory tmpTraj;
                    uint32_t len = e->m_off + e->m_len;
                    m_ts->m_trajIO += std::ceil(len / 4096.0);
                    uint8_t *load;
                    m_ts->loadByteArray(e->m_page, len, &load);
                    uint8_t *data = load + e->m_off;
                    tmpTraj.loadFromByteArray(data);
                    delete[] load;
                    id_type trajid = m_ts->getTrajId(id);
                    if(m_except.count(trajid)>0) return;
                    storeEntry ee;
                    Trajectory inter;
                    tmpTraj.getPartialTrajectory(max(m_query.m_startTime(),e->m_ts),
                                                 min(m_query.m_endTime(),e->m_te),inter);
                    if(inter.m_points.size()==0){
                        auto ee=m_mpq.top();
                        string q = m_query.toString();
                        auto s = m_parts[ee->m_id];
                        std::cerr<<ee->m_id;
                    }
                    insert(trajid, inter, e->m_hasPrev ? 1 : -1, e->m_hasNext ? 1 : -1, ee);
                }

                auto top(){
                    if(!m_mpq.empty()){
                        id_type prevtop=-1;
                        while(m_mpq.top()->m_id!=prevtop){
                            prevtop=m_mpq.top()->m_id;
                            update(prevtop);
                        }
                    }
                    if(!m_mpq.empty()&&m_mpq.top()->m_type==3&&(m_nodespq.empty()||m_mpq.top()->m_minDist<m_nodespq.top()->m_minDist)){
                        return m_mpq.top();
                    }
                    if(!m_nodespq.empty())
                        return m_nodespq.top();
                    return m_mpq.top();
                }

                auto pop(){
                    if(!m_mpq.empty()&&m_mpq.top()->m_type==3&&(m_nodespq.empty()||m_mpq.top()->m_minDist<m_nodespq.top()->m_minDist))
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
                PartsStoreBFMST(Trajectory &traj,double error,TrajStore* ts,bool useMBR)
                        :m_query(traj),m_error(error),m_useMBR(useMBR),m_ts(ts){}
                ~PartsStoreBFMST(){}
            };//PartStoreBFMST

			friend class Node;
			friend class Leaf;
			friend class Index;
			friend class BulkLoader;

			friend std::ostream& operator<<(std::ostream& os, const MBCRTree& t);
		}; // MBCRTree

		std::ostream& operator<<(std::ostream& os, const MBCRTree& t);
	}
}

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
extern bool rsimpli;
namespace SpatialIndex
{
	namespace RTree
	{
		class RTree : public ISpatialIndex
		{
                  //class NNEntry;

		public:
			RTree(IStorageManager&, Tools::PropertySet&);
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


            virtual ~RTree();



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
            bool m_bUsingTrajStore=false;

            TrajStore *m_ts=nullptr;

        private:
			IStorageManager* m_pStorageManager;

			id_type m_rootID, m_headerID;

			RTreeVariant m_treeVariant;


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
			Tools::PointerPool<Node> m_indexPool;
			Tools::PointerPool<Node> m_leafPool;

			std::vector<Tools::SmartPointer<ICommand> > m_writeNodeCommands;
			std::vector<Tools::SmartPointer<ICommand> > m_readNodeCommands;
			std::vector<Tools::SmartPointer<ICommand> > m_deleteNodeCommands;

#ifdef HAVE_PTHREAD_H
			pthread_mutex_t m_lock;
#endif
        public: class simpleData:public IData{
            public:
                simpleData(id_type id,double dist):m_id(id),m_dist(dist){}
                id_type m_id;
                double m_dist;
                virtual Data* clone(){throw Tools::NotSupportedException(".");}
                virtual id_type getIdentifier() const{return m_id;}
                virtual void getShape(IShape** out) const{ *out= nullptr;}
                virtual void getData(uint32_t& len, uint8_t** data) const{len=0;*data= nullptr;}
            };
			class NNEntry
			{
			public:
			    uint32_t m_type=0;
				id_type m_id;
				IEntry* m_pEntry;
				double m_minDist;

                NNEntry(id_type id, IEntry* e, double f, uint32_t type=0)
                        : m_id(id), m_pEntry(e), m_minDist(f),m_type(type) {}
                        
				~NNEntry() {}

				struct ascending : public std::binary_function<NNEntry*, NNEntry*, bool>
				{
					bool operator()(const NNEntry* __x, const NNEntry* __y) const {
                        if(std::isinf(__x->m_minDist)||std::isnan(__x->m_minDist)) return true;
                        if(std::isinf(__y->m_minDist)||std::isnan(__y->m_minDist)) return true;
                        return __x->m_minDist > __y->m_minDist;
                    }
				};
			}; // NNEntry

			class NNComparator : public INearestNeighborComparator
			{
			public:
				double getMinimumDistance(const IShape& query, const IShape& entry)
				{
//				    std::cerr<<"pushed entry with dist"<<query.getMinimumDistance(entry)<<std::endl;
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
            class PartsStore{
            protected:
                class Parts{
                public:
                    PartsStore* m_ps;
                    std::set<id_type> m_missingLeaf, m_loadedLeaf;
                    std::list<Trajectory> m_pTrajs;
                    std::map<std::pair<double,double>,double> m_computedDist;
                    std::map<double,storeEntry> m_entries;
                    double m_calcMin=0;
                    double m_mintime=1e300,m_maxtime=-1e300;
                    bool m_hasPrev=true,m_hasNext=true;
                    double m_computedTime=0,m_loadedTime=0;
                    Parts(PartsStore* ps= nullptr):m_ps(ps){}
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
                bool m_useMBR=false;
                Trajectory m_query;
                double m_error;
                TrajStore* m_ts;
                std::map<id_type ,Parts> m_parts;
                int m_dimension=2;
                std::set<id_type > loadedLeaf;
                void insert(id_type id, Trajectory &r,id_type prev,id_type next,storeEntry &entry){
                    if(m_parts.count(id)==0){
                        m_parts[id]=Parts(this);
                    }
                    m_parts[id].insert(r,prev,next,entry);
                }

                double update(id_type id) {
                    Parts *parts = &m_parts[id];
                    double computedTime = 0;
                    double pd, sum = 0;
                    std::pair<double, double> timeInterval;
                    if (disttype == 0) {
                        //inferred distance(front dist, back dist and mid dist) should be stored as negative values
                        //front dist
                        if (parts->m_mintime > m_query.m_startTime()) {
                            timeInterval.first=m_query.m_startTime();
                            timeInterval.second=parts->m_mintime;
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                if(parts->m_computedDist[timeInterval]>=0)
                                    pd= parts->m_computedDist[timeInterval];
                                else{
                                    pd = -parts->m_computedDist[timeInterval]-1;
                                }
                            } else {
                                if (parts->m_hasPrev) {
                                    pd = m_query.getFrontIED(parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                             parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                             parts->m_pTrajs.front().m_startTime(),
                                                             m_ts->m_maxVelocity);
                                    parts->m_computedDist[timeInterval] = -pd-1;
                                } else {
                                    pd = m_query.getStaticIED(parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                              parts->m_pTrajs.front().m_points.front().m_pCoords[1],
                                                              m_query.m_startTime(), parts->m_mintime);
                                    parts->m_computedDist[timeInterval] = pd;
                                    computedTime += timeInterval.second - timeInterval.first;
                                }
                            }
                            if(parts->m_computedDist[timeInterval]<0&&!m_nodespq.empty())
                                pd=std::max(pd,m_nodespq.top()->m_minDist*
                                               (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                            sum+=pd;
                        }
                        //mid dist
                        const Trajectory *prev= nullptr;
                        for (const auto &traj:parts->m_pTrajs) {
                            //this box
                            timeInterval.first=traj.m_startTime();
                            timeInterval.second=traj.m_endTime();
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                if (parts->m_computedDist[timeInterval] > 0) {
                                    pd = parts->m_computedDist[timeInterval];
                                } else {
                                    pd = traj.getMinimumDistance(m_query);
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                            } else {
                                pd = traj.getMinimumDistance(m_query);
                                parts->m_computedDist[timeInterval] = pd;
                            }
                            sum += pd;
                            computedTime += timeInterval.second - timeInterval.first;
                            //the gap
                            if (traj.m_startTime() != parts->m_pTrajs.front().m_startTime()) {//not first
                                if (prev->m_endTime() < traj.m_startTime()) {
                                    timeInterval.first=prev->m_endTime();
                                    timeInterval.second=traj.m_startTime();
                                    if (parts->m_computedDist.count(timeInterval) > 0) {
                                        if(parts->m_computedDist[timeInterval]>=0)
                                            pd= parts->m_computedDist[timeInterval];
                                        else{
                                            pd = -parts->m_computedDist[timeInterval]-1;
                                        }
                                    } else {
                                        pd = m_query.getMidIED(prev->m_points.back(), traj.m_points.front(), m_ts->m_maxVelocity);
                                        parts->m_computedDist[timeInterval] = -pd-1;
                                    }
                                    if(parts->m_computedDist[timeInterval]<0&&!m_nodespq.empty())
                                        pd=std::max(pd,m_nodespq.top()->m_minDist*
                                                       (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                                    sum+=pd;
                                }
                            }
                            prev = &traj;
                        }
                        //backdist
                        if (parts->m_maxtime < m_query.m_endTime()) {
                            timeInterval.first=parts->m_maxtime;
                            timeInterval.second=m_query.m_endTime();
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                if(parts->m_computedDist[timeInterval]>=0)
                                    pd= parts->m_computedDist[timeInterval];
                                else{
                                    pd = -parts->m_computedDist[timeInterval]-1;
                                }
                            } else {
                                if (parts->m_hasNext) {
                                    pd = m_query.getFrontIED(parts->m_pTrajs.back().m_points.back().m_pCoords[0],
                                                             parts->m_pTrajs.back().m_points.back().m_pCoords[1],
                                                             parts->m_pTrajs.back().m_endTime(), m_ts->m_maxVelocity);
                                    parts->m_computedDist[timeInterval] = -pd-1;
                                } else {
                                    pd = m_query.getStaticIED(parts->m_pTrajs.back().m_points.back().m_pCoords[0],
                                                              parts->m_pTrajs.back().m_points.back().m_pCoords[1], parts->m_maxtime,
                                                              m_query.m_endTime());
                                    parts->m_computedDist[timeInterval] = pd;
                                    computedTime += timeInterval.second - timeInterval.first;
                                }
                            }
                            if(parts->m_computedDist[timeInterval]<0&&!m_nodespq.empty())
                                pd=std::max(pd,m_nodespq.top()->m_minDist*
                                               (timeInterval.second - timeInterval.first)/(m_query.m_endTime()-m_query.m_startTime()));
                            sum+=pd;
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
                        return sum;
                    }
                    else if(disttype==1){
                        double max=0;
                        for(const auto &traj:parts->m_pTrajs){
                            timeInterval = std::make_pair(traj.m_startTime(), traj.m_endTime());
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                pd = parts->m_computedDist[timeInterval];
                            } else {
                                pd = m_query.getMinimumDistance(traj);
                                parts->m_computedDist[timeInterval] = pd;
                            }
                            max=std::max(max,pd);
                        }
                        if(!parts->m_hasPrev){
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                pd = parts->m_computedDist[timeInterval];
                                max=std::max(max,pd);
                            } else {
                                pd = m_query.getStaticMaxSED(parts->m_pTrajs.front().m_points.front().m_pCoords[0],
                                                             parts->m_pTrajs.front().m_points.front().m_pCoords[1],
                                                             m_query.m_startTime(), parts->m_mintime);
                                parts->m_computedDist[timeInterval] = pd;
                                computedTime += timeInterval.second - timeInterval.first;
                                max=std::max(max,pd);
                            }
                        }
                        if(!parts->m_hasNext){
                            if (parts->m_computedDist.count(timeInterval) > 0) {
                                pd = parts->m_computedDist[timeInterval];
                                max=std::max(max,pd);
                            } else {
                                pd = m_query.getStaticMaxSED(parts->m_pTrajs.back().m_points.back().m_pCoords[0],
                                                             parts->m_pTrajs.back().m_points.back().m_pCoords[1], parts->m_maxtime,
                                                             m_query.m_endTime());
                                parts->m_computedDist[timeInterval] = pd;
                                computedTime += timeInterval.second - timeInterval.first;
                                max=std::max(max,pd);
                            }
                        }

                        parts->m_calcMin = max;
                        parts->m_computedTime = computedTime;
                        int type = 2;
                        if (parts->m_loadedTime + 1e-7 >= (m_query.m_endTime() - m_query.m_startTime())) type = 3;
                        max = std::max(0.0, max - m_error);
                        if (m_handlers.count(id) == 0) {
                            auto handle = m_mpq.push(new NNEntry(id,nullptr, max, type));
                            m_handlers[id] = handle;
                        } else {
                            m_mpq.update(m_handlers[id], id, max, type);
                        }
                        return max;
                    }
                    else throw Tools::IllegalStateException("");
                }
            public:
                void loadPartTraj(id_type id, double sTime, double eTime){
                    sTime=std::max(sTime,m_query.m_startTime());
                    eTime=std::min(eTime,m_query.m_endTime());
                    bool hasPrev,hasNext;
                    if(sTime<eTime) {
                        Trajectory ps = m_ts->getTrajByTime(id, sTime, eTime, &hasPrev, &hasNext);
                        storeEntry e;
                        id_type trajid = m_ts->getTrajId(id);
                        insert(trajid, ps, hasPrev ? 1 : -1, hasNext ? 1 : -1, e);
                        double res = update(trajid);
                    }
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
                PartsStore(Trajectory &traj,double error,TrajStore* ts,bool useMBR)
                        :m_query(traj),m_error(error),m_useMBR(useMBR),m_ts(ts){}
                ~PartsStore(){}
            };//PartStore


			friend class Node;
			friend class Leaf;
			friend class Index;
			friend class BulkLoader;

			friend std::ostream& operator<<(std::ostream& os, const RTree& t);
		}; // RTree

		std::ostream& operator<<(std::ostream& os, const RTree& t);
	}
}

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
            bool m_bUsingTrajStore=false;
            bool m_bUsingMBR=true;

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
                uint32_t m_type;
				id_type m_id;
				double m_minDist;

				NNEntry(id_type id, double f, uint32_t type)
				    : m_id(id), m_minDist(f),m_type(type) {
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
			                    if(a->m_minDist < b->m_minDist) return true;
			                    if(a->m_minDist==b->m_minDist) return a->m_type>b->m_type;
			                    return false;
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
			};
            struct storeEntry{
                id_type m_page;
                uint32_t m_off;
                uint32_t m_len;
                bool operator <( const storeEntry& y) const{
                    return std::tie(m_page, m_off, m_len) < std::tie(y.m_page, y.m_off, y.m_len);
                }
            };
			class PartsStore{
            protected:
                class Parts{
                public:
                    PartsStore* m_ps;
                    std::set<id_type> m_missingLeaf, m_loadedLeaf;
                    std::list<RegionPtr> m_mbrs;
                    std::list<MBCPtr> m_mbcs;
                    std::map<std::pair<double,double>,double> m_computedDist;
                    std::map<double,storeEntry> m_entries;
                    double m_calcMin=0;
                    double m_mintime=1e300,m_maxtime=-1e300;
                    bool m_hasPrev=true,m_hasNext=true;
                    double m_computedTime=0,m_loadedTime=0;
                    Parts(PartsStore* ps= nullptr):m_ps(ps){}
                    void insert(RegionPtr &r,id_type prev,id_type next,storeEntry &entry){
                        if(m_mbrs.empty()) m_mbrs.emplace_back(r);
                        else{
                            auto j=m_mbrs.begin();
                            for(;j!=m_mbrs.end()&&(*j)->m_pLow[m_ps->m_dimension]<r->m_pLow[m_ps->m_dimension];j++);
                            m_mbrs.insert(j,r);
                        }
                        m_loadedTime+=r->m_pHigh[m_ps->m_dimension]-r->m_pLow[m_ps->m_dimension];
                        if(r->m_pHigh[m_ps->m_dimension]>m_maxtime){
                            m_maxtime=r->m_pHigh[m_ps->m_dimension];
                            if(next==-1) {
                                m_hasNext=false;
                                m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
                            }
                        }
                        if(r->m_pLow[m_ps->m_dimension]<m_mintime){
                            m_mintime=r->m_pLow[m_ps->m_dimension];
                            if(prev==-1) {
                                m_hasPrev=false;
                                m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
                            }
                        }
                        if(prev>=0&&m_loadedLeaf.count(prev)==0) m_missingLeaf.insert(prev);
                        if(next>=0&&m_loadedLeaf.count(next)==0) m_missingLeaf.insert(next);
                        m_entries[r->m_pLow[m_ps->m_dimension]]=entry;
                    }
                    void insert(MBCPtr &r,id_type prev,id_type next,storeEntry &entry){
                        if(m_mbcs.empty()) m_mbcs.emplace_back(r);
                        else{
                            auto j=m_mbcs.begin();
                            for(;j!=m_mbcs.end()&&(*j)->m_startTime<r->m_startTime;j++);
                            m_mbcs.insert(j,r);
                        }
                        m_loadedTime+=r->m_endTime-r->m_startTime;
                        if(r->m_endTime>m_maxtime){
                            m_maxtime=r->m_endTime;
                            if(next==-1) {
                                m_hasNext=false;
                                m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
                            }
                        }
                        if(r->m_startTime<m_mintime){
                            m_mintime=r->m_startTime;
                            if(prev==-1) {
                                m_hasPrev=false;
                                m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
                            }
                        }
                        if(prev>=0&&m_loadedLeaf.count(prev)==0) m_missingLeaf.insert(prev);
                        if(next>=0&&m_loadedLeaf.count(next)==0) m_missingLeaf.insert(next);
                        m_entries[r->m_startTime]=entry;
                    }
                    void removeCalculated(double ts,double te){}
                };

                std::map<id_type ,MutablePriorityQueue<NNEntry>::handle_type > m_handlers;
                EntryMPQ m_mpq;
                bool m_useMBR=false;
                Trajectory m_query;
                double m_error;
                TrajStore* m_ts;
                std::map<id_type ,Parts> m_parts;
                int m_dimension=2;
                std::set<id_type > loadedLeaf;
                void insert(id_type id, RegionPtr &br,id_type prev,id_type next,storeEntry &entry){
                    if(m_parts.count(id)==0){
                        m_parts[id]=Parts(this);
                    }
                    m_parts[id].insert(br,prev,next,entry);
                }
                void insert(id_type id, MBCPtr &bc,id_type prev,id_type next,storeEntry &entry){
                    if(m_parts.count(id)==0){
                        m_parts[id]=Parts(this);
                    }
                    m_parts[id].insert(bc,prev,next,entry);
                }

                double update(id_type id) {
                    Parts *parts = &m_parts[id];
                    double computedTime = 0;
                    double pd, sum = 0;
                    std::pair<double, double> timeInterval;
                    if (disttype == 0) {
                        if (m_useMBR) {
                            //inferred distance(front dist, back dist and mid dist) should be stored as negative values
                            //front dist
                            if (parts->m_mintime > m_query.m_startTime()) {
                                timeInterval = std::make_pair(m_query.m_startTime(), parts->m_mintime);
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    sum += std::fabs(parts->m_computedDist[timeInterval]);
                                } else {
                                    if (parts->m_hasPrev) {
                                        pd = m_query.getFrontIED(*parts->m_mbrs.front(), m_ts->m_maxVelocity);
                                        parts->m_computedDist[timeInterval] = -pd;
                                    } else {
                                        pd = m_query.getStaticIED(*parts->m_mbrs.front(), m_query.m_startTime(),
                                                                  parts->m_mintime);
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                    sum += pd;
                                }
                            }
                            //mid dist
                            Region *prev;
                            for (const auto &box:parts->m_mbrs) {
                                //this box
                                timeInterval = std::make_pair(box->m_pLow[m_dimension], box->m_pHigh[m_dimension]);
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    if (parts->m_computedDist[timeInterval] > 0) {
                                        pd = parts->m_computedDist[timeInterval];
                                    } else {
                                        pd = m_query.getMinimumDistance(*box);
                                        parts->m_computedDist[timeInterval] = pd;
                                    }
                                } else {
                                    pd = m_query.getMinimumDistance(*box);
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                sum += pd;
                                computedTime += timeInterval.second - timeInterval.first;
                                //the gap
                                if (box->m_pLow[m_dimension] != parts->m_mbrs.front()->m_pLow[m_dimension]) {
                                    if (prev->m_pHigh[m_dimension] < box->m_pLow[m_dimension]) {
                                        timeInterval = std::make_pair(prev->m_pHigh[m_dimension],
                                                                      box->m_pLow[m_dimension]);
                                        if (parts->m_computedDist.count(timeInterval) > 0) {
                                            sum += std::fabs(parts->m_computedDist[timeInterval]);
                                        } else {
                                            pd = m_query.getMidIED(*prev, *box, m_ts->m_maxVelocity);
                                            parts->m_computedDist[timeInterval] = -pd;
                                            sum += pd;
                                        }
                                    }
                                }
                                prev = box.get();
                            }
                            //backdist
                            if (parts->m_maxtime < m_query.m_endTime()) {
                                timeInterval = std::make_pair(parts->m_maxtime, m_query.m_endTime());
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    sum += std::fabs(parts->m_computedDist[timeInterval]);
                                } else {
                                    if (parts->m_hasNext) {
                                        pd = m_query.getFrontIED(*parts->m_mbrs.back(), m_ts->m_maxVelocity);
                                        parts->m_computedDist[timeInterval] = -pd;
                                    } else {
                                        pd = m_query.getStaticIED(*parts->m_mbrs.back(), parts->m_maxtime,
                                                                  m_query.m_endTime());
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                    sum += pd;
                                }
                            }
                        } else {
                            //inferred distance(front dist, back dist and mid dist) should be stored as negative values
                            //front dist
                            if (parts->m_mintime > m_query.m_startTime()) {
                                timeInterval = std::make_pair(m_query.m_startTime(), parts->m_mintime);
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    sum += std::fabs(parts->m_computedDist[timeInterval]);
                                } else {
                                    if (parts->m_hasPrev) {
                                        pd = m_query.getFrontIED(parts->m_mbcs.front()->m_pLow[0],
                                                                 parts->m_mbcs.front()->m_pLow[1],
                                                                 parts->m_mbcs.front()->m_startTime,
                                                                 m_ts->m_maxVelocity);
                                        parts->m_computedDist[timeInterval] = -pd;
                                    } else {
                                        pd = m_query.getStaticIED(parts->m_mbcs.front()->m_pLow[0],
                                                                  parts->m_mbcs.front()->m_pLow[1],
                                                                  m_query.m_startTime(), parts->m_mintime);
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                    sum += pd;
                                }
                            }
                            //mid dist
                            MBC *prev;
                            for (const auto &box:parts->m_mbcs) {
                                //this box
                                timeInterval = std::make_pair(box->m_startTime, box->m_endTime);
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    if (parts->m_computedDist[timeInterval] > 0) {
                                        pd = parts->m_computedDist[timeInterval];
                                    } else {
                                        pd = m_query.getMinimumDistance(*box);
                                        parts->m_computedDist[timeInterval] = pd;
                                    }
                                } else {
                                    pd = m_query.getMinimumDistance(*box);
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                sum += pd;
                                computedTime += timeInterval.second - timeInterval.first;
                                //the gap
                                if (box->m_startTime != parts->m_mbcs.front()->m_startTime) {
                                    if (prev->m_pHigh[m_dimension] < box->m_pLow[m_dimension]) {
                                        timeInterval = std::make_pair(prev->m_pHigh[m_dimension],
                                                                      box->m_pLow[m_dimension]);
                                        if (parts->m_computedDist.count(timeInterval) > 0) {
                                            sum += std::fabs(parts->m_computedDist[timeInterval]);
                                        } else {
                                            pd = m_query.getMidIED(*prev, *box, m_ts->m_maxVelocity);
                                            parts->m_computedDist[timeInterval] = -pd;
                                            sum += pd;
                                        }
                                    }
                                }
                                prev = box.get();
                            }
                            //backdist
                            if (parts->m_maxtime < m_query.m_endTime()) {
                                timeInterval = std::make_pair(parts->m_maxtime, m_query.m_endTime());
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    sum += std::fabs(parts->m_computedDist[timeInterval]);
                                } else {
                                    if (parts->m_hasNext) {
                                        pd = m_query.getFrontIED(parts->m_mbcs.back()->m_pHigh[0],
                                                                 parts->m_mbcs.back()->m_pHigh[1],
                                                                 parts->m_mbcs.back()->m_endTime, m_ts->m_maxVelocity);
                                        parts->m_computedDist[timeInterval] = -pd;
                                    } else {
                                        pd = m_query.getStaticIED(parts->m_mbcs.back()->m_pHigh[0],
                                                                  parts->m_mbcs.back()->m_pHigh[1], parts->m_maxtime,
                                                                  m_query.m_endTime());
                                        parts->m_computedDist[timeInterval] = pd;
                                        computedTime += timeInterval.second - timeInterval.first;
                                    }
                                    sum += pd;
                                }
                            }
                        }
                        parts->m_calcMin = sum;
                        parts->m_computedTime = computedTime;
                        int type = 2;
                        if (parts->m_missingLeaf.empty()) type = 3;
                        sum = std::max(0.0, sum - m_error);
                        if (m_handlers.count(id) == 0) {
                            auto handle = m_mpq.push(new NNEntry(id, sum, type));
                            m_handlers[id] = handle;
                        } else {
                            m_mpq.update(m_handlers[id], id, sum, type);
                        }
                        return sum;
                    }
                    else if(disttype==1){
                        double max=0;
                        if(m_useMBR){
                            for(const auto &box:parts->m_mbrs){
                                timeInterval = std::make_pair(box->m_pLow[m_dimension], box->m_pHigh[m_dimension]);
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd = parts->m_computedDist[timeInterval];
                                } else {
                                    pd = m_query.getMinimumDistance(*box);
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                max=std::max(max,pd);
                            }
                        }else{
                            for(const auto &box:parts->m_mbcs){
                                timeInterval = std::make_pair(box->m_startTime, box->m_endTime);
                                if (parts->m_computedDist.count(timeInterval) > 0) {
                                    pd = parts->m_computedDist[timeInterval];
                                } else {
                                    pd = m_query.getMinimumDistance(*box);
                                    parts->m_computedDist[timeInterval] = pd;
                                }
                                max=std::max(max,pd);
                            }
                        }
                        parts->m_calcMin = max;
                        parts->m_computedTime = computedTime;
                        int type = 2;
                        if (parts->m_missingLeaf.empty()) type = 3;
                        max = std::max(0.0, max - m_error);
                        if (m_handlers.count(id) == 0) {
                            auto handle = m_mpq.push(new NNEntry(id, max, type));
                            m_handlers[id] = handle;
                        } else {
                            m_mpq.update(m_handlers[id], id, max, type);
                        }
                        return max;
                    }
                    else throw Tools::IllegalStateException("");
                }

            public:
			    bool isLoaded(id_type id){ return loadedLeaf.count(id)>0;}
			    void loadLeaf(Node &n){
//                    std::cerr<<"load leaf"<<n.m_nodeMBR<<"\n";
			        loadedLeaf.insert(n.m_identifier);
                    std::vector<id_type > relatedIds;
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
                                insert(trajid, n.m_ptrMBR[i],
                                       (m_query.m_startTime()<bts)?n.m_prevNode[i]:-1
                                        , (m_query.m_endTime()>bte)?n.m_nextNode[i]:-1,
                                        entry);
                                relatedIds.emplace_back(trajid);
                            }

                        }else{
                            double bts=n.m_ptrMBC[i]->m_startTime,bte=n.m_ptrMBC[i]->m_endTime;
                            if(bts>=m_query.m_endTime()||
                               bte<=m_query.m_startTime()){}
                            else {
                                insert(trajid, n.m_ptrMBC[i],
                                       (m_query.m_startTime()<bts)?n.m_prevNode[i]:-1
                                        , (m_query.m_endTime()>bte)?n.m_nextNode[i]:-1,
                                        entry);
                                relatedIds.emplace_back(trajid);
                            }
                        }
                    }
                    for(const auto &rid:relatedIds){
                        m_parts[rid].m_missingLeaf.erase(n.m_identifier);
                        m_parts[rid].m_loadedLeaf.insert(n.m_identifier);
                        update(rid);
                    }
                }
                auto top(){return m_mpq.top();}
                auto pop(){return m_mpq.pop();}
                auto push(NNEntry* e){return m_mpq.push(e);}
                auto empty(){return m_mpq.empty();}
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
                    Trajectory traj;
                    Trajectory tmpTraj;
                    for(const auto &pair:m_parts[id].m_entries){
                        auto e=pair.second;
                        m_ts->m_trajIO+=std::ceil(e.m_len/4096.0);
                        uint32_t len=e.m_off+e.m_len;
                        uint8_t *load = new uint8_t[len];
                        m_ts->loadByteArray(e.m_page,len,&load);
                        uint8_t *data = load+e.m_off;
                        tmpTraj.loadFromByteArray(data);
                        delete[](load);
                        if(traj.m_points.empty()){
                            traj=tmpTraj;
                        }else{
                            traj.linkTrajectory(tmpTraj);
                        }
                    }
                    return traj;
                }
			};


			friend class Node;
			friend class Leaf;
			friend class Index;
			friend class BulkLoader;

			friend std::ostream& operator<<(std::ostream& os, const MBCRTree& t);
		}; // MBCRTree

		std::ostream& operator<<(std::ostream& os, const MBCRTree& t);
	}
}

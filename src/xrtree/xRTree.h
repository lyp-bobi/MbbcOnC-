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
                // todo:DataType                 VT_LONG   Can be BoundingBox and Trajectory


            virtual ~xRTree();

            
			virtual void getIndexProperties(Tools::PropertySet& out) const;
			virtual void getStatistics(IStatistics** out) const;

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

        private:
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
				ValidateEntry(xMBR& r, NodePtr& pNode) : m_parentMBR(r), m_pNode(pNode) {}

				xMBR m_parentMBR;
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


			friend class Node;
			friend class Leaf;
			friend class Index;
			friend class BulkLoader;

			friend std::ostream& operator<<(std::ostream& os, const xRTree& t);
		}; // xRTree




		std::ostream& operator<<(std::ostream& os, const xRTree& t);
	}
}

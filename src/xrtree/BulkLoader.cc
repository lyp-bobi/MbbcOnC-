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

#include <cstring>
#include <stdio.h>
#include <cmath>

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <spatialindex/SpatialIndex.h>

#include "xRTree.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"


using namespace SpatialIndex;
using namespace SpatialIndex::xRTree;

//
// ExternalSorter::Record
//
ExternalSorter::Record::Record()
:{}
ExternalSorter::Record::Record(const xSBBData &shape, uint32_t s, uint32_t level)
: m_s(s),m_level(level),m_b(shape)
{}

ExternalSorter::Record::~Record()
{}

bool ExternalSorter::Record::operator<(const Record& r) const
{
    if (m_s != r.m_s)
        throw Tools::IllegalStateException("ExternalSorter::Record::operator<: Incompatible sorting dimensions.");

    if (m_b.m_b.br->m_pLow(2-m_s)+m_b.m_b.br->m_pHigh(2-m_s) < r.m_b.m_b.br->m_pLow(2-m_s)+r.m_b.m_b.br->m_pHigh(2-m_s))
        return true;
    else
        return false;
}

void ExternalSorter::Record::storeToFile(Tools::TemporaryFile& f)
{
	f.write(m_s);
	f.write(m_level);
    f.write(static_cast<uint64_t>(m_b.m_sbbid));
	f.write(static_cast<uint64_t>(m_b.m_se.m_id));
    f.write(m_b.m_se.m_s);
    f.write(m_b.m_se.m_e);
    f.write(m_b.m_b.toString());
}

void ExternalSorter::Record::loadFromFile(Tools::TemporaryFile& f)
{
	m_s = f.readUInt32();
	m_level=f.readUInt32();
    m_b.m_sbbid =static_cast<id_type>(f.readUInt64());
	m_b.m_se.m_id =static_cast<id_type>(f.readUInt64());
	m_b.m_se.m_s = f.readUInt32();
    m_b.m_se.m_e = f.readUInt32();
	m_b.m_b.loadFromString(f.readString());
}

//
// ExternalSorter
//
ExternalSorter::ExternalSorter(uint32_t u32PageSize, uint32_t u32BufferPages)
: m_bInsertionPhase(true), m_u32PageSize(u32PageSize),
  m_u32BufferPages(u32BufferPages), m_u64TotalEntries(0), m_stI(0)
{
}

ExternalSorter::~ExternalSorter()
{
	for (m_stI = 0; m_stI < m_buffer.size(); ++m_stI) delete m_buffer[m_stI];
}

void ExternalSorter::insert(Record* r)
{
	if (m_bInsertionPhase == false)
		throw Tools::IllegalStateException("ExternalSorter::insert: Input has already been sorted.");

	m_buffer.emplace_back(r);
	++m_u64TotalEntries;

	// this will create the initial, sorted buckets before the
	// external merge sort.
	if (m_buffer.size() >= m_u32PageSize * m_u32BufferPages)
	{
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
		Tools::TemporaryFile* tf = new Tools::TemporaryFile();
		for (size_t j = 0; j < m_buffer.size(); ++j)
		{
			m_buffer[j]->storeToFile(*tf);
			delete m_buffer[j];
		}
		m_buffer.clear();
		tf->rewindForReading();
		m_runs.emplace_back(Tools::SmartPointer<Tools::TemporaryFile>(tf));
	}
}

void ExternalSorter::sort()
{
	if (m_bInsertionPhase == false)
		throw Tools::IllegalStateException("ExternalSorter::sort: Input has already been sorted.");

	if (m_runs.empty())
	{
		// The data fits in main memory. No need to store to disk.
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
		m_bInsertionPhase = false;
		return;
	}

	if (m_buffer.size() > 0)
	{
		// Whatever remained in the buffer (if not filled) needs to be stored
		// as the final bucket.
		std::sort(m_buffer.begin(), m_buffer.end(), Record::SortAscending());
		Tools::TemporaryFile* tf = new Tools::TemporaryFile();
		for (size_t j = 0; j < m_buffer.size(); ++j)
		{
			m_buffer[j]->storeToFile(*tf);
			delete m_buffer[j];
		}
		m_buffer.clear();
		tf->rewindForReading();
		m_runs.emplace_back(Tools::SmartPointer<Tools::TemporaryFile>(tf));
	}

	if (m_runs.size() == 1)
	{
		m_sortedFile = m_runs.front();
	}
	else
	{
		Record* r = 0;

		while (m_runs.size() > 1)
		{
			Tools::SmartPointer<Tools::TemporaryFile> tf(new Tools::TemporaryFile());
			std::vector<Tools::SmartPointer<Tools::TemporaryFile> > buckets;
			std::vector<std::queue<Record*> > buffers;
			std::priority_queue<PQEntry, std::vector<PQEntry>, PQEntry::SortAscending> pq;

			// initialize buffers and priority queue.
			std::list<Tools::SmartPointer<Tools::TemporaryFile> >::iterator it = m_runs.begin();
			for (uint32_t i = 0; i < (std::min)(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages); ++i)
			{
				buckets.push_back(*it);
				buffers.emplace_back(std::queue<Record*>());

				r = new Record();
				r->loadFromFile(**it);
					// a run cannot be empty initially, so this should never fail.
				pq.push(PQEntry(r, i));

				for (uint32_t j = 0; j < m_u32PageSize - 1; ++j)
				{
					// fill the buffer with the rest of the page of records.
					try
					{
						r = new Record();
						r->loadFromFile(**it);
						buffers.back().push(r);
					}
					catch (Tools::EndOfStreamException)
					{
						delete r;
						break;
					}
				}
				++it;
			}

			// exhaust buckets, buffers, and priority queue.
			while (! pq.empty())
			{
				PQEntry e = pq.top(); pq.pop();
				e.m_r->storeToFile(*tf);
				delete e.m_r;

				if (! buckets[e.m_u32Index]->eof() && buffers[e.m_u32Index].empty())
				{
					for (uint32_t j = 0; j < m_u32PageSize; ++j)
					{
						try
						{
							r = new Record();
							r->loadFromFile(*buckets[e.m_u32Index]);
							buffers[e.m_u32Index].push(r);
						}
						catch (Tools::EndOfStreamException)
						{
							delete r;
							break;
						}
					}
				}

				if (! buffers[e.m_u32Index].empty())
				{
					e.m_r = buffers[e.m_u32Index].front();
					buffers[e.m_u32Index].pop();
					pq.push(e);
				}
			}

			tf->rewindForReading();

			// check if another pass is needed.
			uint32_t u32Count = std::min(static_cast<uint32_t>(m_runs.size()), m_u32BufferPages);
			for (uint32_t i = 0; i < u32Count; ++i)
			{
				m_runs.pop_front();
			}

			if (m_runs.size() == 0)
			{
				m_sortedFile = tf;
				break;
			}
			else
			{
				m_runs.emplace_back(tf);
			}
		}
	}

	m_bInsertionPhase = false;
}

ExternalSorter::Record* ExternalSorter::getNextRecord()
{
	if (m_bInsertionPhase == true)
		throw Tools::IllegalStateException("ExternalSorter::getNextRecord: Input has not been sorted yet.");

	Record* ret;

	if (m_sortedFile.get() == 0)
	{
		if (m_stI < m_buffer.size())
		{
			ret = m_buffer[m_stI];
			m_buffer[m_stI] = 0;
			++m_stI;
		}
		else
			throw Tools::EndOfStreamException("");
	}
	else
	{
		ret = new Record();
		ret->loadFromFile(*m_sortedFile);
	}

	return ret;
}

inline uint64_t ExternalSorter::getTotalEntries() const
{
	return m_u64TotalEntries;
}

//
// BulkLoader
//
void BulkLoader::bulkLoadUsingSTR(
	SpatialIndex::xRTree::xRTree* pTree,
	IDataStream& stream,
	uint32_t bindex,
	uint32_t bleaf,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	if (! stream.hasNext())
		throw Tools::IllegalArgumentException(
			"xRTree::BulkLoader::bulkLoadUsingSTR: Empty data stream given."
		);

	NodePtr n = pTree->readNode(pTree->m_rootID);
	pTree->deleteNode(n.get());

	#ifndef NDEBUG
	std::cerr << "xRTree::BulkLoader: Sorting data." << std::endl;
	#endif

	Tools::SmartPointer<ExternalSorter> es = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));

	calcuTime[0]=calcuTime[1]=calcuTime[2]=0;
    int scount=0;
	while (stream.hasNext())
	{
	    scount++;
	    auto d=stream.getNext();
        xSBBData* d2 = dynamic_cast<xSBBData*>(d);
        if (d2 == nullptr) {
            throw Tools::IllegalArgumentException(
                    "bulkLoadUsingSTR: xRTree bulk load expects SpatialIndex::xRTree::Data entries."
            );
        }
        es->insert(new ExternalSorter::Record(*d2,0,0));
        delete d2;
	}
	std::cerr<<"start sorting the"<<es->getTotalEntries()<<"entries\n";
	es->sort();
    std::cerr<<"sorted.\n";
	pTree->m_stats.m_u64Data = es->getTotalEntries();

	// create index levels.
	uint32_t level = 0;

	while (true)
	{
		#ifndef NDEBUG
		std::cerr << "xRTree::BulkLoader: Building level " << level << std::endl;
		#endif

		pTree->m_stats.m_nodesInLevel.emplace_back(0);

		Tools::SmartPointer<ExternalSorter> es2 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
		createLevel(pTree, es, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
		es = es2;

		if (es->getTotalEntries() == 1) break;
		es->sort();
	}
#ifdef Tristat
	std::cerr<<"Lx is "<<calcuTime[0]/calcuTime[1]<<"\n";
    std::cerr<<"Lt is "<<calcuTime[2]/calcuTime[1]<<"\n";
#endif
	pTree->m_stats.m_u32TreeHeight = level;
	pTree->storeHeader();
}

void BulkLoader::createLevel(
	SpatialIndex::xRTree::xRTree* pTree,
	Tools::SmartPointer<ExternalSorter> es,
	uint32_t dimension,
	uint32_t bleaf,
	uint32_t bindex,
	uint32_t level,
	Tools::SmartPointer<ExternalSorter> es2,
	uint32_t pageSize,
	uint32_t numberOfPages
) {
	uint64_t b = (level == 0) ? bleaf : bindex;
    uint64_t P = static_cast<uint64_t>(std::ceil(static_cast<double>(es->getTotalEntries()) / static_cast<double>(b)));
    int remainDim;
    remainDim = pTree->m_dimension - dimension;
    uint64_t S;
    double ltc;
    auto stat=trajStat::instance();
    if(dimension==0){
        double nt = pow(P * sq(stat->vv())*sq(stat->P) / sq((M_PI * sq(stat->Sr))),1.0/3);
        std::cerr<<"nt is "<<nt<<"\t rest is "<<sqrt(P/nt)<<endl;
        S = max(1,int(floor(nt)));
        if(S>P) S=1;
    }
    else {
        S = static_cast<uint64_t>(ceil(pow(static_cast<double>(P), 1.0 / remainDim)));
    }
//    std::cerr<<"crtlvl at "<<level<<"\t"<<dimension<<"\t"<<P<<"\t"<<S<<"\n";
	if (S == 1 || remainDim==1 || S * b == es->getTotalEntries())
	{
		std::vector<ExternalSorter::Record*> node;
		ExternalSorter::Record* r;

		while (true)
		{
			try { r = es->getNextRecord(); } catch (Tools::EndOfStreamException) { break; }
			node.emplace_back(r);

			if (node.size() == b)
			{
				Node* n = createNode(pTree, node, level);
				node.clear();
				pTree->writeNode(n);
				es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, 0, 0,level));
				pTree->m_rootID = n->m_identifier;
                if(level==0){
                    //state the storage place of bounding boxes
                    for(int i=0;i<n->m_children;i++){
                        pTree->m_ts->m_part2node[n->m_pIdentifier[i]]=n->m_identifier;
                    }
                }
                else if(level==1){
                    //linking the bounding boxes
                    for(int i=0;i<n->m_children;i++){
                        NodePtr child=pTree->readNode(n->m_pIdentifier[i]);
                        for(int j=0;j<child->m_children;j++){
                            id_type id=child->m_pIdentifier[j];
                            xStore::Entry *entry;
                            if(pTree->m_ts->m_stream!= nullptr){
                                auto it=pTree->m_ts->
                                        m_dentries.search(std::to_string(id));
                                DiskMultiMap::MultiMapTuple tuple = *it;
                                entry=new xStore::Entry(tuple.value);
                            }else{
                                entry=pTree->m_ts->m_entries[id];
                            }
                            id_type pvId=entry->m_pvId;
                            if(pvId>=0){
                                auto store=pTree->m_ts->m_part2node[pvId];
                                child->m_prevNode[j]=store;
//                                std::cerr<<"linked"<<id<<"to"<<pvId<<"\n";
                            }else{
                                child->m_prevNode[j]=-1;
                            }
                            id_type ntId=entry->m_ntId;
                            if(ntId>=0){
                                auto store=pTree->m_ts->m_part2node[ntId];
                                child->m_nextNode[j]=store;
//                                std::cerr<<"linked"<<id<<"to"<<ntId<<"\n";
                            }else{
                                child->m_nextNode[j]=-1;
                            }
                            child->m_pageNum[j]=entry->m_page;
                            child->m_pageOff[j]=entry->m_start;
                            child->m_dataLen[j]=entry->m_len;
//                            std::cerr<<store->m_page<<" "<<store->m_start<<" "<<store->m_len<<"\n";
                        }
                        pTree->writeNode(child.get());
                    }
                }
					// special case when the root has exactly bindex entries.
				delete n;
			}
		}

		if (! node.empty())
		{
			Node* n = createNode(pTree, node, level);
			pTree->writeNode(n);
			es2->insert(new ExternalSorter::Record(n->m_nodeMBR, n->m_identifier, 0, 0, 0,level));
			pTree->m_rootID = n->m_identifier;
            if(level==0){
                //state the storage place of bounding boxes
                for(int i=0;i<n->m_children;i++){
                    pTree->m_ts->m_part2node[n->m_pIdentifier[i]]=n->m_identifier;
                }
            }
            else if(level==1){
                //linking the bounding boxes
                for(int i=0;i<n->m_children;i++){
                    NodePtr child=pTree->readNode(n->m_pIdentifier[i]);
                    for(int j=0;j<child->m_children;j++){
                        id_type id=child->m_pIdentifier[j];
                        xStore::Entry *entry;
                        if(pTree->m_ts->m_stream!= nullptr){
                            auto it=pTree->m_ts->
                                    m_dentries.search(std::to_string(id));
                            DiskMultiMap::MultiMapTuple tuple = *it;
                            entry=new xStore::Entry(tuple.value);
                        }else{
                            entry=pTree->m_ts->m_entries[id];
                        }
                        id_type pvId=entry->m_pvId;
                        if(pvId>=0){
                            auto store=pTree->m_ts->m_part2node[pvId];
                            child->m_prevNode[j]=store;
//                                std::cerr<<"linked"<<id<<"to"<<pvId<<"\n";
                        }else{
                            child->m_prevNode[j]=-1;
                        }
                        id_type ntId=entry->m_ntId;
                        if(ntId>=0){
                            auto store=pTree->m_ts->m_part2node[ntId];
                            child->m_nextNode[j]=store;
//                                std::cerr<<"linked"<<id<<"to"<<ntId<<"\n";
                        }else{
                            child->m_nextNode[j]=-1;
                        }
                        child->m_pageNum[j]=entry->m_page;
                        child->m_pageOff[j]=entry->m_start;
                        child->m_dataLen[j]=entry->m_len;
//                            std::cerr<<store->m_page<<" "<<store->m_start<<" "<<store->m_len<<"\n";
                    }
                    pTree->writeNode(child.get());
                }
            }
			delete n;
		}
	}
	else
	{
	    auto stat=trajStat::instance();
	    double curt=stat->mint+ltc;
		bool bMore = true;
//        int count1=0;
        ExternalSorter::Record* pR;
        try { pR = es->getNextRecord(); }
        catch (Tools::EndOfStreamException) {
            bMore = false;
        }
		while (bMore)
		{
		    double lastTime;
			Tools::SmartPointer<ExternalSorter> es3 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
            pR->m_s = dimension + 1;
            es3->insert(pR);
			for (uint64_t i = 1; i < b * ceil(1.0 * P / S); ++i)
            {
                try { pR = es->getNextRecord(); }
                catch (Tools::EndOfStreamException) {
                    bMore = false;
                    break;
                }
                pR->m_s = dimension + 1;
                lastTime=(pR->m_r.m_pLow[2]+pR->m_r.m_pHigh[2])/2;
                es3->insert(pR);
            }
            try { pR = es->getNextRecord(); }
            catch (Tools::EndOfStreamException) {
                bMore = false;
            }
//            while(bMore&&(pR->m_r.m_pLow[2]+pR->m_r.m_pHigh[2])/2==lastTime){
//                pR->m_s = dimension + 1;
//                es3->insert(pR);
//                try { pR = es->getNextRecord(); }
//                catch (Tools::EndOfStreamException) {
//                    bMore = false;
//                    break;
//                }
//            }

            if(es3->getTotalEntries()>0) {
//                if(level==0&&dimension==0) {
//                    std::cerr << "entries is"<<es3->getTotalEntries()<<"\n";
//                }
                es3->sort();
                createLevel(pTree, es3, dimension + 1, bleaf, bindex, level, es2, pageSize, numberOfPages);
            }

		}
	}
}

Node* BulkLoader::createNode(SpatialIndex::xRTree::xRTree* pTree, std::vector<ExternalSorter::Record*>& e, uint32_t level)
{
	Node* n;
	if (level == 0) n = new Leaf(pTree, -1);
	else n = new Index(pTree, -1, level);

	for (size_t cChild = 0; cChild < e.size(); ++cChild)
	{
        if (level == 0) n->insertEntry(e[cChild]->m_len, e[cChild]->m_pData,e[cChild]->m_r, e[cChild]->m_xMBC, e[cChild]->m_id);
        else n->insertEntry(e[cChild]->m_len, e[cChild]->m_pData, e[cChild]->m_r, e[cChild]->m_id);
	    e[cChild]->m_pData = 0;
		delete e[cChild];
	}
#ifdef Tristat
    if(level==0){
        calcuTime[0]+=(n->m_nodeMBR.m_pHigh[0]-n->m_nodeMBR.m_pLow[0]);
        calcuTime[2]+=(n->m_nodeMBR.m_pHigh[2]-n->m_nodeMBR.m_pLow[2]);
        calcuTime[1]+=1;
    }
#endif
	return n;
}
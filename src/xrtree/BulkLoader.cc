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
#include "Node.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"

bool bulkloadt= true;
using namespace SpatialIndex;
using namespace SpatialIndex::xRTreeNsp;


//
// ExternalSorter::Record
//
ExternalSorter::Record::Record(){}
ExternalSorter::Record::Record(const xSBBData &shape, uint32_t s, uint32_t level)
: m_s(s), m_level(level), m_sbbd(shape)
{
    if(!shape.m_b.hasbr) throw Tools::IllegalStateException("SBBs should have MBR");
    m_r = shape.m_b.br;
    m_id = shape.m_sbbid;
}
ExternalSorter::Record::Record(const id_type &id, const xMBR &r, uint32_t s, uint32_t level)
        : m_s(s),m_level(level),m_id(id), m_r(r)
{}

ExternalSorter::Record::~Record()
{}

bool ExternalSorter::Record::operator<(const Record& r) const
{
    if (m_s != r.m_s)
        throw Tools::IllegalStateException("ExternalSorter::Record::operator<: Incompatible sorting dimensions.");
    Record *ll = const_cast<Record*>(this), *rr = const_cast<Record*>(&r);
    return ll->m_r.m_pLow(2 - m_s) + ll->m_r.m_pHigh(2 - m_s) < rr->m_r.m_pLow(2 - m_s) + rr->m_r.m_pHigh(2 - m_s);
}

void ExternalSorter::Record::storeToFile(Tools::TemporaryFile& f)
{
    f.write(static_cast<uint64_t>(m_id));
	f.write(m_s);
	f.write(m_level);
	f.write(m_r.m_xmin);
    f.write(m_r.m_xmax);
    f.write(m_r.m_ymin);
    f.write(m_r.m_ymax);
    f.write(m_r.m_tmin);
    f.write(m_r.m_tmax);
	if(m_level==0){//leafnode
        f.write(static_cast<uint64_t>(m_sbbd.m_sbbid));
        f.write(static_cast<uint64_t>(m_sbbd.m_se.m_id));
        f.write(static_cast<uint32_t>(m_sbbd.m_hasPrev ? 1 : 0));
        f.write(static_cast<uint32_t>(m_sbbd.m_hasNext ? 1 : 0));
        f.write(m_sbbd.m_se.m_s);
        f.write(m_sbbd.m_se.m_e);
        f.write(m_sbbd.m_b.toString());
	}
}

void ExternalSorter::Record::loadFromFile(Tools::TemporaryFile& f)
{
    m_id =static_cast<id_type>(f.readUInt64());
	m_s = f.readUInt32();
	m_level=f.readUInt32();
    m_r.m_xmin = f.readDouble();
    m_r.m_xmax = f.readDouble();
    m_r.m_ymin = f.readDouble();
    m_r.m_ymax = f.readDouble();
    m_r.m_tmin = f.readDouble();
    m_r.m_tmax = f.readDouble();
	if(m_level==0) {//if leaf node
        m_sbbd.m_sbbid = static_cast<id_type>(f.readUInt64());
        m_sbbd.m_se.m_id = static_cast<id_type>(f.readUInt64());
        m_sbbd.m_hasPrev = (f.readUInt32() == 1);
        m_sbbd.m_hasNext = (f.readUInt32() == 1);
        m_sbbd.m_se.m_s = f.readUInt32();
        m_sbbd.m_se.m_e = f.readUInt32();
        m_sbbd.m_b.loadFromString(f.readString());
    }
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
	SpatialIndex::xRTreeNsp::xRTree* pTree,
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

	std::cerr << "xRTree::BulkLoader: Sorting data." << std::endl;

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
		if(bulkloadt&&level==0)
            createLevelT(pTree, es, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
        else
            createLevel(pTree, es, 0, bleaf, bindex, level++, es2, pageSize, numberOfPages);
		es = es2;

		if (es->getTotalEntries() == 1) break;
		es->sort();
	}
    std::cerr<<"Lx is "<<calcuTime[0]/calcuTime[1]<<"\n";
    std::cerr<<"Lt is "<<calcuTime[2]/calcuTime[1]<<"\n";
//	m_part2node.clear();
//    std::map<id_type,id_type> empty;
//    m_part2node.swap(empty);
	pTree->m_stats.m_u32TreeHeight = level;
	pTree->storeHeader();
}

void BulkLoader::createLevel(
	SpatialIndex::xRTreeNsp::xRTree* pTree,
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

    if(dimension==0){
        double nt = pow(P * sq(tjstat->vv()) * sq(tjstat->P) / sq((M_PI * sq(tjstat->Sr))), 1.0 / 3);
        /*debug*/
        cerr<<"P is "<<P<<"vv is "<< tjstat->vv()<<"Time is"<<tjstat->P <<"Sr is"<<tjstat->Sr<<endl;
        std::cerr<<"nt is "<<nt<<"\t rest is "<<sqrt(P/nt)<<endl;
        S = max(1,int(floor(nt)));
        if(S>P) S=1;
    }
    else {
        S = static_cast<uint64_t>(ceil(pow(static_cast<double>(P), 1.0 / remainDim)));
    }
//    std::cerr<<"crtlvl at "<<level<<"\t"<<dimension<<"\t"<<P<<"\t"<<S<<"\n";
	if (remainDim==1||P<=S ||S * b > es->getTotalEntries())//||P<=S|| P/S * b == es->getTotalEntries()
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
				es2->insert(new ExternalSorter::Record(n->m_identifier,n->m_nodeMBR,0,level));
				pTree->m_rootID = n->m_identifier;
				if(pTree->m_bStoringLinks) {
                    if (level == 0) {
                        //state the storage place of bounding boxes
                        for (int i = 0; i < n->m_children; i++) {
                            m_part2node[n->m_pIdentifier[i]] = n->m_identifier;
                        }
                    } else if (level == 1) {
                        //linking the bounding boxes
                        for (int i = 0; i < n->m_children; i++) {
                            NodePtr child = pTree->readNode(n->m_pIdentifier[i]);
                            for (int j = 0; j < child->m_children; j++) {
                                id_type id = child->m_pIdentifier[j];
                                pair<bool, bool> pvnt = make_pair(child->m_prevNode[j], child->m_nextNode[j]);
                                if (pvnt.first) {
                                    auto store = m_part2node[id - 1];
                                    child->m_prevNode[j] = store;
//                                std::cerr<<"linked"<<id<<"to"<<pvId<<"\n";
                                } else {
                                    child->m_prevNode[j] = -1;
                                }
                                if (pvnt.second) {
                                    auto store = m_part2node[id + 1];
                                    child->m_nextNode[j] = store;
//                                std::cerr<<"linked"<<id<<"to"<<ntId<<"\n";
                                } else {
                                    child->m_nextNode[j] = -1;
                                }
//                            std::cerr<<store->m_page<<" "<<store->m_start<<" "<<store->m_len<<"\n";
                            }
                            pTree->writeNode(child.get());
                        }
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
			es2->insert(new ExternalSorter::Record( n->m_identifier,n->m_nodeMBR, 0,level));
			pTree->m_rootID = n->m_identifier;
			if(pTree->m_bStoringLinks) {
                if (level == 0) {
                    //state the storage place of bounding boxes
                    for (int i = 0; i < n->m_children; i++) {
                        m_part2node[n->m_pIdentifier[i]] = n->m_identifier;
                    }
                } else if (level == 1) {
                    //linking the bounding boxes
                    for (int i = 0; i < n->m_children; i++) {
                        NodePtr child = pTree->readNode(n->m_pIdentifier[i]);
                        for (int j = 0; j < child->m_children; j++) {
                            id_type id = child->m_pIdentifier[j];
                            pair<bool, bool> pvnt = make_pair(child->m_prevNode[j], child->m_nextNode[j]);
                            if (pvnt.first) {
                                auto store = m_part2node[id - 1];
                                child->m_prevNode[j] = store;
//                                std::cerr<<"linked"<<id<<"to"<<pvId<<"\n";
                            } else {
                                child->m_prevNode[j] = -1;
                            }
                            if (pvnt.second) {
                                auto store = m_part2node[id + 1];
                                child->m_nextNode[j] = store;
//                                std::cerr<<"linked"<<id<<"to"<<ntId<<"\n";
                            } else {
                                child->m_nextNode[j] = -1;
                            }
//                            std::cerr<<store->m_page<<" "<<store->m_start<<" "<<store->m_len<<"\n";
                        }
                        pTree->writeNode(child.get());
                    }
                }
            }
			delete n;
		}
	}
	else
	{
		bool bMore = true;
        ExternalSorter::Record* pR;
		while (bMore)
		{
			Tools::SmartPointer<ExternalSorter> es3 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
			for (uint64_t i = 0; i < b * ceil(1.0 * P / S); ++i)
            {
                try { pR = es->getNextRecord(); }
                catch (Tools::EndOfStreamException) {
                    bMore = false;
                    break;
                }
                pR->m_s = dimension + 1;
                es3->insert(pR);
            }

            if(es3->getTotalEntries()>0) {
                es3->sort();
                createLevel(pTree, es3, dimension + 1, bleaf, bindex, level, es2, pageSize, numberOfPages);
            }

		}
	}
}

void BulkLoader::createLevelT(
        SpatialIndex::xRTreeNsp::xRTree* pTree,
        Tools::SmartPointer<ExternalSorter> es,
        uint32_t dimension,
        uint32_t bleaf,
        uint32_t bindex,
        uint32_t level,
        Tools::SmartPointer<ExternalSorter> es2,
        uint32_t pageSize,
        uint32_t numberOfPages
) {
    assert(dimension==0&&level==0);
    double curt, lastt;

    bool bMore = true;
    ExternalSorter::Record* pR;
    int total = es->getTotalEntries();
    try { pR = es->getNextRecord(); }
    catch (Tools::EndOfStreamException) {
        bMore = false;
    }
    curt = (pR->m_sbbd.m_b.br.m_tmax + pR->m_sbbd.m_b.br.m_tmin)/2;
    pR->m_s = dimension + 1;
    lastt = curt + pTree->m_bt;
    while (bMore)
    {
        Tools::SmartPointer<ExternalSorter> es3 = Tools::SmartPointer<ExternalSorter>(new ExternalSorter(pageSize, numberOfPages));
        while(es3->getTotalEntries()< (bleaf*bindex) && bMore) {
            lastt = lastt + pTree->m_bt;
            while (curt <= lastt) {
                es3->insert(pR);
                try { pR = es->getNextRecord(); }
                catch (Tools::EndOfStreamException) {
                    bMore = false;
                    break;
                }
                pR->m_s = dimension + 1;
                if(curt > (pR->m_sbbd.m_b.br.m_tmax + pR->m_sbbd.m_b.br.m_tmin)/2)
                    break;
                curt = (pR->m_sbbd.m_b.br.m_tmax + pR->m_sbbd.m_b.br.m_tmin)/2;
            }

        }

        if(es3->getTotalEntries()>0) {
            total -= es3->getTotalEntries();
//            std::cerr<<"slab"<<lastt<<"\t"<<es3->getTotalEntries()<<"\t"<<total<<endl;
            es3->sort();
            createLevel(pTree, es3, dimension + 1, bleaf, bindex, level, es2, pageSize, numberOfPages);
        }

    }
}

Node* BulkLoader::createNode(SpatialIndex::xRTreeNsp::xRTree* pTree, std::vector<ExternalSorter::Record*>& e, uint32_t level)
{
	Node* n;
	if (level == 0) n = new Leaf(pTree, -1);
	else n = new Index(pTree, -1, level);
	for (size_t cChild = 0; cChild < e.size(); ++cChild)
	{
        if (level == 0) n->insertEntry(e[cChild]->m_sbbd.m_b.br, e[cChild]->m_sbbd.m_sbbid,
                                       &(e[cChild]->m_sbbd));
        else n->insertEntry(e[cChild]->m_r,e[cChild]->m_id);
		delete e[cChild];
	}
    if(level==0){
        calcuTime[0]+=(n->m_nodeMBR.m_xmax-n->m_nodeMBR.m_xmin);
        calcuTime[2]+=(n->m_nodeMBR.m_tmax-n->m_nodeMBR.m_tmin);
        calcuTime[1]+=1;
    }
	return n;
}

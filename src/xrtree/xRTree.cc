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
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>
#include "Node.h"
#include "Leaf.h"
#include "Index.h"
#include "BulkLoader.h"
#include "xRTree.h"


bool bUsingSimp=true;
bool bUsingSBBD=false;
int sbb=0,sb=0;

using namespace SpatialIndex::xRTree;
using namespace SpatialIndex;

SpatialIndex::xRTree::Data::Data(uint32_t len, uint8_t* pData, xMBC& r, xMBR &rg, id_type id)
	: m_id(id), m_xMBC(r),m_mbr(rg), m_pData(0), m_dataLength(len)
{
	if (m_dataLength > 0)
	{
		m_pData = new uint8_t[m_dataLength];
		memcpy(m_pData, pData, m_dataLength);
	}
}

SpatialIndex::xRTree::Data::Data(uint32_t len, uint8_t* pData, xMBC& r, id_type id)
        : m_id(id), m_xMBC(r), m_pData(0), m_dataLength(len)
{
    m_mbr.makeInfinite(2);
    if (m_dataLength > 0)
    {
        m_pData = new uint8_t[m_dataLength];
        memcpy(m_pData, pData, m_dataLength);
    }
}

SpatialIndex::xRTree::Data::~Data()
{
	delete[] m_pData;
}

SpatialIndex::xRTree::Data* SpatialIndex::xRTree::Data::clone()
{
	return new Data(m_dataLength, m_pData, m_xMBC,m_mbr, m_id);
}

id_type SpatialIndex::xRTree::Data::getIdentifier() const
{
	return m_id;
}

void SpatialIndex::xRTree::Data::getShape(IShape** out) const
{
	*out = new xMBC(m_xMBC);
}

void SpatialIndex::xRTree::Data::getData(uint32_t& len, uint8_t** data) const
{
	len = m_dataLength;
	*data = 0;

	if (m_dataLength > 0)
	{
		*data = new uint8_t[m_dataLength];
		memcpy(*data, m_pData, m_dataLength);
	}
}

uint32_t SpatialIndex::xRTree::Data::getByteArraySize() const
{
	return
		sizeof(id_type) +
		sizeof(uint32_t) +
		m_dataLength +
		m_mbr.getByteArraySize()+
		m_xMBC.getByteArraySize();
}

void SpatialIndex::xRTree::Data::loadFromByteArray(const uint8_t* ptr)
{
	memcpy(&m_id, ptr, sizeof(id_type));
	ptr += sizeof(id_type);

	delete[] m_pData;
	m_pData = 0;

	memcpy(&m_dataLength, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	if (m_dataLength > 0)
	{
		m_pData = new uint8_t[m_dataLength];
		memcpy(m_pData, ptr, m_dataLength);
		ptr += m_dataLength;
	}
    m_mbr.loadFromByteArray(ptr);
	ptr+=m_mbr.getByteArraySize();
	m_xMBC.loadFromByteArray(ptr);
}

void SpatialIndex::xRTree::Data::storeToByteArray(uint8_t** data, uint32_t& len)
{
	// it is thread safe this way.
	uint32_t xMBRsize,xMBRsize2;
	uint8_t* xMBRdata = 0,*xMBRdata2=0;
    m_mbr.storeToByteArray(&xMBRdata, xMBRsize);
	m_xMBC.storeToByteArray(&xMBRdata2, xMBRsize2);

	len = sizeof(id_type) + sizeof(uint32_t) + m_dataLength + xMBRsize;

	*data = new uint8_t[len];
	uint8_t* ptr = *data;

	memcpy(ptr, &m_id, sizeof(id_type));
	ptr += sizeof(id_type);
	memcpy(ptr, &m_dataLength, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	if (m_dataLength > 0)
	{
		memcpy(ptr, m_pData, m_dataLength);
		ptr += m_dataLength;
	}

	memcpy(ptr, xMBRdata, xMBRsize);
    ptr += xMBRsize;
    memcpy(ptr, xMBRdata2, xMBRsize2);
	delete[] xMBRdata;
	delete[] xMBRdata2;
	// ptr += xMBRsize;
}

SpatialIndex::ISpatialIndex* SpatialIndex::xRTree::returnxRTree(SpatialIndex::IStorageManager& sm, Tools::PropertySet& ps)
{
	SpatialIndex::ISpatialIndex* si = new SpatialIndex::xRTree::xRTree(sm, ps);
	return si;
}

SpatialIndex::ISpatialIndex* SpatialIndex::xRTree::createNewxRTree(
	SpatialIndex::IStorageManager& sm,
	double fillFactor,
	uint32_t indexCapacity,
	uint32_t leafCapacity,
	uint32_t dimension,
	xRTreeVariant rv,
	id_type& indexIdentifier)
{
	Tools::Variant var;
	Tools::PropertySet ps;

	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = fillFactor;
	ps.setProperty("FillFactor", var);

	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = indexCapacity;
	ps.setProperty("IndexCapacity", var);

	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = leafCapacity;
	ps.setProperty("LeafCapacity", var);

	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = dimension;
	ps.setProperty("Dimension", var);

	var.m_varType = Tools::VT_LONG;
	var.m_val.lVal = rv;
	ps.setProperty("TreeVariant", var);

	ISpatialIndex* ret = returnxRTree(sm, ps);

	var.m_varType = Tools::VT_LONGLONG;
	var = ps.getProperty("IndexIdentifier");
	indexIdentifier = var.m_val.llVal;

	return ret;
}

SpatialIndex::ISpatialIndex* SpatialIndex::xRTree::createAndBulkLoadNewxRTree(
	BulkLoadMethod m,
	IDataStream& stream,
	SpatialIndex::IStorageManager& sm,
	double fillFactor,
	uint32_t indexCapacity,
	uint32_t leafCapacity,
	uint32_t dimension,
	SpatialIndex::xRTree::xRTreeVariant rv,
	id_type& indexIdentifier,
	bool useMBR //false;
	)
{
	SpatialIndex::ISpatialIndex* tree = createNewxRTree(sm, fillFactor, indexCapacity, leafCapacity, dimension, rv, indexIdentifier);
    xRTree* r= static_cast<xRTree*>(tree);
    xStore *ts= dynamic_cast<xStore*>(&sm);
    if(ts!= nullptr){
        r->m_DataType=TrajectoryType;
        r->m_ts=ts;
    }
    r->m_bUsingMBR=useMBR;
	uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity)));
	uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity)));

	SpatialIndex::xRTree::BulkLoader bl;

	switch (m)
	{
	case BLM_STR:
		bl.bulkLoadUsingSTR(static_cast<xRTree*>(tree), stream, bindex, bleaf, 10000, 10000);
		break;
	default:
		throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Unknown bulk load method.");
		break;
	}

	return tree;
}

SpatialIndex::ISpatialIndex* SpatialIndex::xRTree::createAndBulkLoadNewxRTree(
	BulkLoadMethod m,
	IDataStream& stream,
	SpatialIndex::IStorageManager& sm,
	Tools::PropertySet& ps,
	id_type& indexIdentifier)
{
	Tools::Variant var;
	xRTreeVariant rv(RV_LINEAR);
	double fillFactor(0.0);
	uint32_t indexCapacity(0);
	uint32_t leafCapacity(0);
	uint32_t dimension(0);
	uint32_t pageSize(0);
	uint32_t numberOfPages(0);

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			var.m_val.lVal != RV_QUADRATIC &&
			var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property TreeVariant must be Tools::VT_LONG and of xRTreeVariant type");

		rv = static_cast<xRTreeVariant>(var.m_val.lVal);
	}

	// fill factor
	// it cannot be larger than 50%, since linear and quadratic split algorithms
	// require assigning to both nodes the same number of entries.
	var = ps.getProperty("FillFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
	    if (var.m_varType != Tools::VT_DOUBLE)
            throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property FillFactor was not of type Tools::VT_DOUBLE");

        if (var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property FillFactor was less than 0.0");

        if (((rv == RV_LINEAR || rv == RV_QUADRATIC) && var.m_val.dblVal > 0.5))
            throw Tools::IllegalArgumentException( "createAndBulkLoadNewxRTree: Property FillFactor must be in range (0.0, 0.5) for LINEAR or QUADRATIC index types");
        if ( var.m_val.dblVal >= 1.0)
            throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property FillFactor must be in range (0.0, 1.0) for RSTAR index type");
		fillFactor = var.m_val.dblVal;
	}

	// index capacity
	var = ps.getProperty("IndexCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property IndexCapacity must be Tools::VT_ULONG and >= 4");

		indexCapacity = var.m_val.ulVal;
	}

	// leaf capacity
	var = ps.getProperty("LeafCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG )// || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property LeafCapacity must be Tools::VT_ULONG and >= 4");

		leafCapacity = var.m_val.ulVal;
	}

	// dimension
	var = ps.getProperty("Dimension");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property Dimension must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property Dimension must be greater than 1");

		dimension = var.m_val.ulVal;
	}

	// page size
	var = ps.getProperty("ExternalSortBufferPageSize");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property ExternalSortBufferPageSize must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property ExternalSortBufferPageSize must be greater than 1");

		pageSize = var.m_val.ulVal;
	}

	// number of pages
	var = ps.getProperty("ExternalSortBufferTotalPages");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property ExternalSortBufferTotalPages must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Property ExternalSortBufferTotalPages must be greater than 1");

		numberOfPages = var.m_val.ulVal;
	}

	SpatialIndex::ISpatialIndex* tree = createNewxRTree(sm, fillFactor, indexCapacity, leafCapacity, dimension, rv, indexIdentifier);

	uint32_t bindex = static_cast<uint32_t>(std::floor(static_cast<double>(indexCapacity * fillFactor)));
	uint32_t bleaf = static_cast<uint32_t>(std::floor(static_cast<double>(leafCapacity * fillFactor)));

	SpatialIndex::xRTree::BulkLoader bl;

	switch (m)
	{
	case BLM_STR:
		bl.bulkLoadUsingSTR(static_cast<xRTree*>(tree), stream, bindex, bleaf, pageSize, numberOfPages);
		break;
	default:
		throw Tools::IllegalArgumentException("createAndBulkLoadNewxRTree: Unknown bulk load method.");
		break;
	}

	return tree;
}

SpatialIndex::ISpatialIndex* SpatialIndex::xRTree::loadxRTree(IStorageManager& sm, id_type indexIdentifier)
{
	Tools::Variant var;
	Tools::PropertySet ps;

	var.m_varType = Tools::VT_LONGLONG;
	var.m_val.llVal = indexIdentifier;
	ps.setProperty("IndexIdentifier", var);

	return returnxRTree(sm, ps);
}

SpatialIndex::xRTree::xRTree::xRTree(IStorageManager& sm, Tools::PropertySet& ps) :
	m_pStorageManager(&sm),
	m_rootID(StorageManager::NewPage),
	m_headerID(StorageManager::NewPage),
	m_treeVariant(RV_RSTAR),
	m_fillFactor(0.7),
	m_indexCapacity(100),
	m_leafCapacity(100),
	m_nearMinimumOverlapFactor(32),
	m_splitDistributionFactor(0.4),
	m_reinsertFactor(0.3),
	m_dimension(2),
	m_bTightMBRs(true),
	m_xPointPool(500),
	m_xMBRPool(1000),
    m_xMBCPool(1000),
	m_indexPool(100),
	m_leafPool(100)
{
#ifdef HAVE_PTHREAD_H
	pthread_mutex_init(&m_lock, NULL);
#endif

	Tools::Variant var = ps.getProperty("IndexIdentifier");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType == Tools::VT_LONGLONG) m_headerID = var.m_val.llVal;
		else if (var.m_varType == Tools::VT_LONG) m_headerID = var.m_val.lVal;
			// for backward compatibility only.
		else throw Tools::IllegalArgumentException("xRTree: Property IndexIdentifier must be Tools::VT_LONGLONG");

		initOld(ps);
	}
	else
	{
		initNew(ps);
		var.m_varType = Tools::VT_LONGLONG;
		var.m_val.llVal = m_headerID;
		ps.setProperty("IndexIdentifier", var);
	}
}

SpatialIndex::xRTree::xRTree::~xRTree()
{
#ifdef HAVE_PTHREAD_H
	pthread_mutex_destroy(&m_lock);
#endif

	storeHeader();
}

//
// ISpatialIndex interface
//

void SpatialIndex::xRTree::xRTree::insertData(uint32_t len, const uint8_t* pData, const IShape& shape, id_type id)
{
	if (shape.getDimension() != m_dimension) throw Tools::IllegalArgumentException("insertData: Shape has the wrong number of dimensions.");

#ifdef HAVE_PTHREAD_H
	Tools::LockGuard lock(&m_lock);
#endif

	// convert the shape into a xMBR (R-Trees index xMBRs only; i.e., approximations of the shapes).
	xMBRPtr mbr = m_xMBRPool.acquire();
	shape.getMBR(*mbr);

	uint8_t* buffer = 0;

	if (len > 0)
	{
		buffer = new uint8_t[len];
		memcpy(buffer, pData, len);
	}

	insertData_impl(len, buffer, *mbr, id);
		// the buffer is stored in the tree. Do not delete here.
}

bool SpatialIndex::xRTree::xRTree::deleteData(const IShape& shape, id_type id)
{
	if (shape.getDimension() != m_dimension) throw Tools::IllegalArgumentException("deleteData: Shape has the wrong number of dimensions.");

#ifdef HAVE_PTHREAD_H
	Tools::LockGuard lock(&m_lock);
#endif

	xMBRPtr mbr = m_xMBRPool.acquire();
	shape.getMBR(*mbr);
	bool ret = deleteData_impl(*mbr, id);

	return ret;
}


void SpatialIndex::xRTree::xRTree::containsWhatQuery(const IShape& query, IVisitor& v)
{
	throw;
}
void SpatialIndex::xRTree::xRTree::intersectsWithQuery(const IShape& query, IVisitor& v)
{
//	if (query.getDimension() != m_dimension)
//	    throw Tools::IllegalArgumentException("intersectsWithQuery: Shape has the wrong number of dimensions.");
	rangeQuery(IntersectionQuery, query, v);
}

void SpatialIndex::xRTree::xRTree::xPointLocationQuery(const xPoint& query, IVisitor& v)
{
	if (query.m_dimension != m_dimension) throw Tools::IllegalArgumentException("xPointLocationQuery: Shape has the wrong number of dimensions.");
	xMBR r(query, query);
	rangeQuery(IntersectionQuery, r, v);
}



void SpatialIndex::xRTree::xRTree::nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v, INearestNeighborComparator& nnc)
{
//	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("nearestNeighborQuery: Shape has the wrong number of dimensions.");
    const Trajectory *queryTraj;
    if(m_DataType==TrajectoryType) {
        queryTraj = dynamic_cast<const Trajectory *>(&query);
        if (queryTraj == nullptr||queryTraj->m_xPoints.size()<2) {
            std::cerr << "bad query traj\n";
            return;
        }
    }
    Trajectory simpleTraj;
    Trajectory ssTraj;
    double delta=0, ssdelta= 0;
    if(bUsingSBBD==true && bUsingSimp == true) {
//        auto stat=trajStat::instance();
        int segnum = std::ceil((queryTraj->m_endTime() - queryTraj->m_startTime()) / (m_ts->m_timeCount/m_ts->m_segCount));
        segnum=std::max(segnum,10);
        vector<vector<STxPoint>> simpseg;
        try {
            simpseg = Trajectory::simplifyWithRDPN(queryTraj->m_xPoints,
                    std::min(segnum,int(std::sqrt(queryTraj->m_xPoints.size()))));
        }
        catch (...){
            std::cerr<<"RDPN query has some problem with\n"<< queryTraj<<"with"<<std::min(segnum,int(std::sqrt(queryTraj->m_xPoints.size())))<<"\n";
        }
        vector<STxPoint> simpp;
        for (const auto &s:simpseg) {
            simpp.emplace_back(s.front());
        }
        simpp.emplace_back(simpseg.back().back());
        simpleTraj=Trajectory(simpp);
        delta = queryTraj->getMinimumDistance(simpleTraj);
        simpp.clear();
        simpp.emplace_back(queryTraj->m_xPoints[0]);
        simpp.emplace_back(queryTraj->m_xPoints[queryTraj->m_xPoints.size()-1]);
        ssTraj=Trajectory(simpp);
        ssdelta = queryTraj->getMinimumDistance(ssTraj);
    }else{
        simpleTraj=*queryTraj;
    }

#ifdef HAVE_PTHREAD_H
	Tools::LockGuard lock(&m_lock);
#endif

    double knearest = 0.0;
    int iternum = 0;
    /*SBB-Driven*/
    if(bUsingSBBD == true) {
        PartsStore ps(simpleTraj, delta, m_ts, m_bUsingMBR);
        ps.push(new NNEntry(m_rootID, 0, 0));

        uint32_t count = 0;

        std::map<id_type, int> insertedTrajId;
        while (!ps.empty()) {
            iternum++;
            NNEntry *pFirst = ps.top();

            // report all nearest neighbors with equal greatest distances.
            // (neighbors can be more than k, if many happen to have the same greatest distance).
            if (count >= k && pFirst->m_minDist > knearest) {
//            std::cerr<<"find minDist"<<knearest<<"\n";
                break;
            }
            switch (pFirst->m_type) {
                case 0: {//inner node
                    ps.pop();
                    NodePtr n = readNode(pFirst->m_id);
                    v.visitNode(*n);
                    for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                        double pd;
                        if(n->m_level>=2 && bUsingSimp){
                            pd = std::max(0.0, ssTraj.getNodeMinimumDistance(*(n->m_ptrMBR[cChild]),
                                                                                 m_ts->m_maxVelocity) - ssdelta);
                        }else{
                            pd = std::max(0.0, simpleTraj.getNodeMinimumDistance(*(n->m_ptrMBR[cChild]),
                                                                                   m_ts->m_maxVelocity) - delta);
                        }
                        if (pd < 1e300) {
                            if (n->m_level == 1)
                                ps.push(new NNEntry(n->m_pIdentifier[cChild], pd, 1));
                            else
                                ps.push(new NNEntry(n->m_pIdentifier[cChild], pd, 0));
                        }
                    }
                    delete pFirst;
//                n.relinquish();
                    break;
                }
                case 1: {//leaf node
                    ps.pop();
                    if (!ps.isLoaded(pFirst->m_id)) {
                        NodePtr n = readNode(pFirst->m_id);
                        m_ts->m_leaf1 += 1;
                        v.visitNode(*n);
                        ps.loadLeaf(*n);
//                    n.relinquish();
                    }
                    delete pFirst;
                    break;
                }
                case 2: {//incomplete bounding
                    id_type missing = ps.getOneMissingPart(pFirst->m_id);
                    NodePtr n = readNode(missing);
                    v.visitNode(*n);
                    ps.loadLeaf(*n,pFirst->m_minDist);
                    m_ts->m_leaf2 += 1;
//                n.relinquish();
                    break;
                }
                case 3: {//complete bounding
                    ps.pop();
                    Trajectory traj = ps.getTraj(pFirst->m_id);
//                std::cerr<<"trajIO"<<m_ts->m_trajIO<<"\n";
//                std::cerr<<"getTraj"<<traj<<"\n";
//                Trajectory traj = m_ts->getTrajByTime(pFirst->m_id, queryTraj->m_startTime(), queryTraj->m_endTime());
                    ps.push(new NNEntry(pFirst->m_id, queryTraj->getMinimumDistance(traj), 4));
                    delete pFirst;
                    break;
                }
                case 4: {//exact traj
                    ps.pop();
                    ++(m_stats.m_u64QueryResults);
                    ++count;
                    knearest = pFirst->m_minDist;
                    simpleData d(pFirst->m_id, pFirst->m_minDist);
                    v.visitData(d);
                    delete pFirst;
                    break;
                }
                default:
                    throw Tools::IllegalStateException("illegal NNEntry state");
            }
        }
        ps.clean();
    }
    else{/*BFMST*/
        PartsStoreBFMST ps(simpleTraj,0,m_ts,true);
        string str = queryTraj->toString();
        ps.push(new NNEntry(m_rootID, 0, 0));

        uint32_t count = 0;
        double knearest = 0.0;
        int iternum=0;
        std::map<id_type ,int> insertedTrajId;

        while (! ps.empty()) {
            iternum++;
            NNEntry *pFirst = ps.top();
//            std::cerr<<"pfirst\t"<<pFirst->m_minDist<<"\n";
            // report all nearest neighbors with equal greatest distances.
            // (neighbors can be more than k, if many happen to have the same greatest distance).
            if (count >= k && pFirst->m_minDist > knearest) {
//            std::cerr<<"find minDist"<<knearest<<"\n";
                break;
            }
            switch (pFirst->m_type) {
                case 0: {//inner and leaf node
                    ps.pop();
                    NodePtr n = readNode(pFirst->m_id);
                    v.visitNode(*n);
                    for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                        if (n->m_level == 0) {
                            double pd;
                            double ts, te;
                            if (m_bUsingMBR) {
                                pd = std::max(0.0, simpleTraj.getNodeMinimumDistance(*(n->m_ptrMBR[cChild]),0));
                                ts = n->m_ptrMBR[cChild]->m_pLow[2];
                                te = n->m_ptrMBR[cChild]->m_pHigh[2];
                            } else {
                                double d = simpleTraj.getMinimumDistance(*(n->m_ptrxMBC[cChild]))/
                                ((n->m_ptrxMBC[cChild])->m_endTime- (n->m_ptrxMBC[cChild])->m_startTime)
                                *(simpleTraj.m_endTime()-simpleTraj.m_startTime());
                                pd = std::max(0.0, d);
                                ts = n->m_ptrxMBC[cChild]->m_startTime;
                                te = n->m_ptrxMBC[cChild]->m_endTime;
                            }
                            leafInfo *e = new leafInfo();
                            e->m_page = n->m_pageNum[cChild];
                            e->m_off = n->m_pageOff[cChild];
                            e->m_len = n->m_dataLen[cChild];
                            e->m_hasPrev = (n->m_prevNode[cChild] != -1);
                            e->m_hasNext = (n->m_nextNode[cChild] != -1);
                            e->m_ts = ts;
                            e->m_te = te;
                            ps.push(new NNEntry(n->m_pIdentifier[cChild], e, pd, 1));
                        } else {
                            double pd = std::max(0.0, simpleTraj.getNodeMinimumDistance(*(n->m_ptrMBR[cChild]),
                                                                                        m_ts->m_maxVelocity) - delta);
                            ps.push(new NNEntry(n->m_pIdentifier[cChild], nullptr, pd, 0));
                        }
                    }
                    delete pFirst;
                    break;
                }
                case 1: {//leaf node
                    ps.pop();
                    ps.loadPartTraj(pFirst->m_id, pFirst->m_pEntry,pFirst->m_minDist);
                    delete pFirst->m_pEntry;
                    delete pFirst;
                    break;
                }
                case 3: {
                    ps.pop();
                    ++(m_stats.m_u64QueryResults);
                    ++count;
                    knearest = pFirst->m_minDist;
                    simpleData d(pFirst->m_id, pFirst->m_minDist);
                    v.visitData(d);
                    delete pFirst;
                    break;
                }

                default:
                    throw Tools::IllegalStateException("illegal NNEntry state" + std::to_string(pFirst->m_type));
            }
        }
        ps.clean();
    }
//    std::cout<<"knearest is"<<knearest<<std::endl;
//    std::cerr<<"iternum is "<<iternum<<"\n";
    m_stats.m_doubleExactQueryResults+=knearest;
}

void SpatialIndex::xRTree::xRTree::nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v)
{
//	if (query.getDimension() != m_dimension) throw Tools::IllegalArgumentException("nearestNeighborQuery: Shape has the wrong number of dimensions.");
	NNComparator nnc;
	nearestNeighborQuery(k, query, v, nnc);
}


void SpatialIndex::xRTree::xRTree::selfJoinQuery(const IShape& query, IVisitor& v)
{
	if (query.getDimension() != m_dimension)
		throw Tools::IllegalArgumentException("selfJoinQuery: Shape has the wrong number of dimensions.");

#ifdef HAVE_PTHREAD_H
	Tools::LockGuard lock(&m_lock);
#endif

	xMBRPtr mbr = m_xMBRPool.acquire();
	query.getMBR(*mbr);
	selfJoinQuery(m_rootID, m_rootID, *mbr, v);
}

void SpatialIndex::xRTree::xRTree::queryStrategy(IQueryStrategy& qs)
{
#ifdef HAVE_PTHREAD_H
	Tools::LockGuard lock(&m_lock);
#endif

	id_type next = m_rootID;
	bool hasNext = true;

	while (hasNext)
	{
		NodePtr n = readNode(next);
		qs.getNextEntry(*n, next, hasNext);
	}
}

void SpatialIndex::xRTree::xRTree::getIndexProperties(Tools::PropertySet& out) const
{
	Tools::Variant var;

	// dimension
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_dimension;
	out.setProperty("Dimension", var);

	// index capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_indexCapacity;
	out.setProperty("IndexCapacity", var);

	// leaf capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_leafCapacity;
	out.setProperty("LeafCapacity", var);

	// R-tree variant
	var.m_varType = Tools::VT_LONG;
	var.m_val.lVal = m_treeVariant;
	out.setProperty("TreeVariant", var);

	// fill factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_fillFactor;
	out.setProperty("FillFactor", var);

	// near minimum overlap factor
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_nearMinimumOverlapFactor;
	out.setProperty("NearMinimumOverlapFactor", var);

	// split distribution factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_splitDistributionFactor;
	out.setProperty("SplitDistributionFactor", var);

	// reinsert factor
	var.m_varType = Tools::VT_DOUBLE;
	var.m_val.dblVal = m_reinsertFactor;
	out.setProperty("ReinsertFactor", var);

	// tight MBRs
	var.m_varType = Tools::VT_BOOL;
	var.m_val.blVal = m_bTightMBRs;
	out.setProperty("EnsureTightMBRs", var);

	// index pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_indexPool.getCapacity();
	out.setProperty("IndexPoolCapacity", var);

	// leaf pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_leafPool.getCapacity();
	out.setProperty("LeafPoolCapacity", var);

	// xMBR pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_xMBRPool.getCapacity();
	out.setProperty("xMBRPoolCapacity", var);

	// xPoint pool capacity
	var.m_varType = Tools::VT_ULONG;
	var.m_val.ulVal = m_xPointPool.getCapacity();
	out.setProperty("xPointPoolCapacity", var);
}

void SpatialIndex::xRTree::xRTree::addCommand(ICommand* pCommand, CommandType ct)
{
	switch (ct)
	{
		case CT_NODEREAD:
			m_readNodeCommands.emplace_back(Tools::SmartPointer<ICommand>(pCommand));
			break;
		case CT_NODEWRITE:
			m_writeNodeCommands.emplace_back(Tools::SmartPointer<ICommand>(pCommand));
			break;
		case CT_NODEDELETE:
			m_deleteNodeCommands.emplace_back(Tools::SmartPointer<ICommand>(pCommand));
			break;
	}
}

bool SpatialIndex::xRTree::xRTree::isIndexValid()
{
	bool ret = true;
	std::stack<ValidateEntry> st;
	NodePtr root = readNode(m_rootID);

	if (root->m_level != m_stats.m_u32TreeHeight - 1)
	{
		std::cerr << "Invalid tree height." << std::endl;
		return false;
	}

	std::map<uint32_t, uint32_t> nodesInLevel;
	nodesInLevel.insert(std::pair<uint32_t, uint32_t>(root->m_level, 1));

	ValidateEntry e(root->m_nodeMBR, root);
	st.push(e);

	while (! st.empty())
	{
		e = st.top(); st.pop();

		xMBR tmpxMBR;
		tmpxMBR = m_infinitexMBR;

		for (uint32_t cDim = 0; cDim < tmpxMBR.m_dimension; ++cDim)
		{
			tmpxMBR.m_pLow[cDim] = std::numeric_limits<double>::max();
			tmpxMBR.m_pHigh[cDim] = -std::numeric_limits<double>::max();

			for (uint32_t cChild = 0; cChild < e.m_pNode->m_children; ++cChild)
			{
				tmpxMBR.m_pLow[cDim] = std::min(tmpxMBR.m_pLow[cDim], e.m_pNode->m_ptrMBR[cChild]->m_pLow[cDim]);
				tmpxMBR.m_pHigh[cDim] = std::max(tmpxMBR.m_pHigh[cDim], e.m_pNode->m_ptrMBR[cChild]->m_pHigh[cDim]);
			}
		}

		if (! (tmpxMBR == e.m_pNode->m_nodeMBR))
		{
			std::cerr << "Invalid parent information." << std::endl;
			ret = false;
		}
		else if (! (tmpxMBR == e.m_parentMBR))
		{
			std::cerr << "Error in parent." << std::endl;
			ret = false;
		}

		if (e.m_pNode->m_level != 0)
		{
			for (uint32_t cChild = 0; cChild < e.m_pNode->m_children; ++cChild)
			{
				NodePtr ptrN = readNode(e.m_pNode->m_pIdentifier[cChild]);
				ValidateEntry tmpEntry(*(e.m_pNode->m_ptrMBR[cChild]), ptrN);

				std::map<uint32_t, uint32_t>::iterator itNodes = nodesInLevel.find(tmpEntry.m_pNode->m_level);

				if (itNodes == nodesInLevel.end())
				{
					nodesInLevel.insert(std::pair<uint32_t, uint32_t>(tmpEntry.m_pNode->m_level, 1l));
				}
				else
				{
					nodesInLevel[tmpEntry.m_pNode->m_level] = nodesInLevel[tmpEntry.m_pNode->m_level] + 1;
				}

				st.push(tmpEntry);
			}
		}
	}

	uint32_t nodes = 0;
	for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
	{
		if (nodesInLevel[cLevel] != m_stats.m_nodesInLevel[cLevel])
		{
			std::cerr << "Invalid nodesInLevel information." << std::endl;
			ret = false;
		}

		nodes += m_stats.m_nodesInLevel[cLevel];
	}

	if (nodes != m_stats.m_u32Nodes)
	{
		std::cerr << "Invalid number of nodes information." << std::endl;
		ret = false;
	}

	return ret;
}

void SpatialIndex::xRTree::xRTree::getStatistics(IStatistics** out) const
{
	*out = new Statistics(m_stats);
}

void SpatialIndex::xRTree::xRTree::initNew(Tools::PropertySet& ps)
{
	Tools::Variant var;

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			var.m_val.lVal != RV_QUADRATIC &&
			var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("initNew: Property TreeVariant must be Tools::VT_LONG and of xRTreeVariant type");

		m_treeVariant = static_cast<xRTreeVariant>(var.m_val.lVal);
	}

	// fill factor
	// it cannot be larger than 50%, since linear and quadratic split algorithms
	// require assigning to both nodes the same number of entries.
	var = ps.getProperty("FillFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
	    if (var.m_varType != Tools::VT_DOUBLE)
            throw Tools::IllegalArgumentException("initNew: Property FillFactor was not of type Tools::VT_DOUBLE");

        if (var.m_val.dblVal <= 0.0)
            throw Tools::IllegalArgumentException("initNew: Property FillFactor was less than 0.0");

        if (((m_treeVariant == RV_LINEAR || m_treeVariant == RV_QUADRATIC) && var.m_val.dblVal > 0.5))
            throw Tools::IllegalArgumentException(  "initNew: Property FillFactor must be in range "
                                                    "(0.0, 0.5) for LINEAR or QUADRATIC index types");
        if ( var.m_val.dblVal >= 1.0)
            throw Tools::IllegalArgumentException(  "initNew: Property FillFactor must be in range "
                                                    "(0.0, 1.0) for RSTAR index type");
		m_fillFactor = var.m_val.dblVal;
	}

	// index capacity
	var = ps.getProperty("IndexCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG || var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("initNew: Property IndexCapacity must be Tools::VT_ULONG and >= 4");

		m_indexCapacity = var.m_val.ulVal;
	}

	// leaf capacity
	var = ps.getProperty("LeafCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG )//|| var.m_val.ulVal < 4)
			throw Tools::IllegalArgumentException("initNew: Property LeafCapacity must be Tools::VT_ULONG and >= 4");

		m_leafCapacity = var.m_val.ulVal;
	}

	// near minimum overlap factor
	var = ps.getProperty("NearMinimumOverlapFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_ULONG ||
			var.m_val.ulVal < 1 ||
			var.m_val.ulVal > m_indexCapacity ||
			var.m_val.ulVal > m_leafCapacity)
			throw Tools::IllegalArgumentException("initNew: Property NearMinimumOverlapFactor must be Tools::VT_ULONG and less than both index and leaf capacities");

		m_nearMinimumOverlapFactor = var.m_val.ulVal;
	}

	// split distribution factor
	var = ps.getProperty("SplitDistributionFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_DOUBLE ||
			var.m_val.dblVal <= 0.0 ||
			var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initNew: Property SplitDistributionFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_splitDistributionFactor = var.m_val.dblVal;
	}

	// reinsert factor
	var = ps.getProperty("ReinsertFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_DOUBLE ||
			var.m_val.dblVal <= 0.0 ||
			var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initNew: Property ReinsertFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_reinsertFactor = var.m_val.dblVal;
	}

	// dimension
	var = ps.getProperty("Dimension");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property Dimension must be Tools::VT_ULONG");
		if (var.m_val.ulVal <= 1)
			throw Tools::IllegalArgumentException("initNew: Property Dimension must be greater than 1");

		m_dimension = var.m_val.ulVal;
	}

	// tight MBRs
	var = ps.getProperty("EnsureTightMBRs");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_BOOL)
			throw Tools::IllegalArgumentException("initNew: Property EnsureTightMBRs must be Tools::VT_BOOL");

		m_bTightMBRs = var.m_val.blVal;
	}

	// index pool capacity
	var = ps.getProperty("IndexPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property IndexPoolCapacity must be Tools::VT_ULONG");

		m_indexPool.setCapacity(var.m_val.ulVal);
	}

	// leaf pool capacity
	var = ps.getProperty("LeafPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property LeafPoolCapacity must be Tools::VT_ULONG");

		m_leafPool.setCapacity(var.m_val.ulVal);
	}

	// xMBR pool capacity
	var = ps.getProperty("xMBRPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property xMBRPoolCapacity must be Tools::VT_ULONG");

		m_xMBRPool.setCapacity(var.m_val.ulVal);
	}

	// xPoint pool capacity
	var = ps.getProperty("xPointPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG)
			throw Tools::IllegalArgumentException("initNew: Property xPointPoolCapacity must be Tools::VT_ULONG");

		m_xPointPool.setCapacity(var.m_val.ulVal);
	}

	m_infinitexMBR.makeInfinite(m_dimension);

	m_stats.m_u32TreeHeight = 1;
	m_stats.m_nodesInLevel.emplace_back(0);

	Leaf root(this, -1);
	m_rootID = writeNode(&root);

	storeHeader();
}

void SpatialIndex::xRTree::xRTree::initOld(Tools::PropertySet& ps)
{
	loadHeader();

	// only some of the properties may be changed.
	// the rest are just ignored.

	Tools::Variant var;

	// tree variant
	var = ps.getProperty("TreeVariant");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_LONG ||
			(var.m_val.lVal != RV_LINEAR &&
			 var.m_val.lVal != RV_QUADRATIC &&
			 var.m_val.lVal != RV_RSTAR))
			throw Tools::IllegalArgumentException("initOld: Property TreeVariant must be Tools::VT_LONG and of xRTreeVariant type");

		m_treeVariant = static_cast<xRTreeVariant>(var.m_val.lVal);
	}

	// near minimum overlap factor
	var = ps.getProperty("NearMinimumOverlapFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (
			var.m_varType != Tools::VT_ULONG ||
			var.m_val.ulVal < 1 ||
			var.m_val.ulVal > m_indexCapacity ||
			var.m_val.ulVal > m_leafCapacity)
			throw Tools::IllegalArgumentException("initOld: Property NearMinimumOverlapFactor must be Tools::VT_ULONG and less than both index and leaf capacities");

		m_nearMinimumOverlapFactor = var.m_val.ulVal;
	}

	// split distribution factor
	var = ps.getProperty("SplitDistributionFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_DOUBLE || var.m_val.dblVal <= 0.0 || var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initOld: Property SplitDistributionFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_splitDistributionFactor = var.m_val.dblVal;
	}

	// reinsert factor
	var = ps.getProperty("ReinsertFactor");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_DOUBLE || var.m_val.dblVal <= 0.0 || var.m_val.dblVal >= 1.0)
			throw Tools::IllegalArgumentException("initOld: Property ReinsertFactor must be Tools::VT_DOUBLE and in (0.0, 1.0)");

		m_reinsertFactor = var.m_val.dblVal;
	}

	// tight MBRs
	var = ps.getProperty("EnsureTightMBRs");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_BOOL) throw Tools::IllegalArgumentException("initOld: Property EnsureTightMBRs must be Tools::VT_BOOL");

		m_bTightMBRs = var.m_val.blVal;
	}

	// index pool capacity
	var = ps.getProperty("IndexPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property IndexPoolCapacity must be Tools::VT_ULONG");

		m_indexPool.setCapacity(var.m_val.ulVal);
	}

	// leaf pool capacity
	var = ps.getProperty("LeafPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property LeafPoolCapacity must be Tools::VT_ULONG");

		m_leafPool.setCapacity(var.m_val.ulVal);
	}

	// xMBR pool capacity
	var = ps.getProperty("xMBRPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property xMBRPoolCapacity must be Tools::VT_ULONG");

		m_xMBRPool.setCapacity(var.m_val.ulVal);
	}

	// xPoint pool capacity
	var = ps.getProperty("xPointPoolCapacity");
	if (var.m_varType != Tools::VT_EMPTY)
	{
		if (var.m_varType != Tools::VT_ULONG) throw Tools::IllegalArgumentException("initOld: Property xPointPoolCapacity must be Tools::VT_ULONG");

		m_xPointPool.setCapacity(var.m_val.ulVal);
	}

	m_infinitexMBR.makeInfinite(m_dimension);
}

void SpatialIndex::xRTree::xRTree::storeHeader()
{
	const uint32_t headerSize =
		sizeof(id_type) +						// m_rootID
		sizeof(xRTreeVariant) +					// m_treeVariant
		sizeof(double) +						// m_fillFactor
		sizeof(uint32_t) +						// m_indexCapacity
		sizeof(uint32_t) +						// m_leafCapacity
		sizeof(uint32_t) +						// m_nearMinimumOverlapFactor
		sizeof(double) +						// m_splitDistributionFactor
		sizeof(double) +						// m_reinsertFactor
		sizeof(uint32_t) +						// m_dimension
		sizeof(char) +							// m_bTightMBRs
		sizeof(uint32_t) +						// m_stats.m_nodes
		sizeof(uint64_t) +						// m_stats.m_data
		sizeof(uint32_t) +						// m_stats.m_treeHeight
		m_stats.m_u32TreeHeight * sizeof(uint32_t);	// m_stats.m_nodesInLevel

	uint8_t* header = new uint8_t[headerSize];
	uint8_t* ptr = header;

	memcpy(ptr, &m_rootID, sizeof(id_type));
	ptr += sizeof(id_type);
	memcpy(ptr, &m_treeVariant, sizeof(xRTreeVariant));
	ptr += sizeof(xRTreeVariant);
	memcpy(ptr, &m_fillFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_indexCapacity, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_leafCapacity, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_nearMinimumOverlapFactor, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_splitDistributionFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_reinsertFactor, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, &m_dimension, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	char c = (char) m_bTightMBRs;
	memcpy(ptr, &c, sizeof(char));
	ptr += sizeof(char);
	memcpy(ptr, &(m_stats.m_u32Nodes), sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &(m_stats.m_u64Data), sizeof(uint64_t));
	ptr += sizeof(uint64_t);
	memcpy(ptr, &(m_stats.m_u32TreeHeight), sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
	{
		memcpy(ptr, &(m_stats.m_nodesInLevel[cLevel]), sizeof(uint32_t));
		ptr += sizeof(uint32_t);
	}

	m_pStorageManager->storeByteArray(m_headerID, headerSize, header);

	delete[] header;
}

void SpatialIndex::xRTree::xRTree::loadHeader()
{
	uint32_t headerSize;
	uint8_t* header = 0;
	m_pStorageManager->loadByteArray(m_headerID, headerSize, &header);

	uint8_t* ptr = header;

	memcpy(&m_rootID, ptr, sizeof(id_type));
	ptr += sizeof(id_type);
	memcpy(&m_treeVariant, ptr, sizeof(xRTreeVariant));
	ptr += sizeof(xRTreeVariant);
	memcpy(&m_fillFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_indexCapacity, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_leafCapacity, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_nearMinimumOverlapFactor, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_splitDistributionFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_reinsertFactor, ptr, sizeof(double));
	ptr += sizeof(double);
	memcpy(&m_dimension, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	char c;
	memcpy(&c, ptr, sizeof(char));
	m_bTightMBRs = (c != 0);
	ptr += sizeof(char);
	memcpy(&(m_stats.m_u32Nodes), ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&(m_stats.m_u64Data), ptr, sizeof(uint64_t));
	ptr += sizeof(uint64_t);
	memcpy(&(m_stats.m_u32TreeHeight), ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	for (uint32_t cLevel = 0; cLevel < m_stats.m_u32TreeHeight; ++cLevel)
	{
		uint32_t cNodes;
		memcpy(&cNodes, ptr, sizeof(uint32_t));
		ptr += sizeof(uint32_t);
		m_stats.m_nodesInLevel.emplace_back(cNodes);
	}

	delete[] header;
}

void SpatialIndex::xRTree::xRTree::insertData_impl(uint32_t dataLength, uint8_t* pData, xMBR& mbr, id_type id)
{
	assert(mbr.getDimension() == m_dimension);

	std::stack<id_type> pathBuffer;
	uint8_t* overflowTable = 0;

	try
	{
		NodePtr root = readNode(m_rootID);

		overflowTable = new uint8_t[root->m_level];
		memset(overflowTable, 0, root->m_level);

		NodePtr l = root->chooseSubtree(mbr, 0, pathBuffer);
		if (l.get() == root.get())
		{
			assert(root.unique());
			root.relinquish();
		}
		l->insertData(dataLength, pData, mbr, id, pathBuffer, overflowTable);

		delete[] overflowTable;
		++(m_stats.m_u64Data);
	}
	catch (...)
	{
		delete[] overflowTable;
		throw;
	}
}

void SpatialIndex::xRTree::xRTree::insertData_impl(uint32_t dataLength, uint8_t* pData, xMBR& mbr, id_type id, uint32_t level, uint8_t* overflowTable)
{
	assert(mbr.getDimension() == m_dimension);

	std::stack<id_type> pathBuffer;
	NodePtr root = readNode(m_rootID);
	NodePtr n = root->chooseSubtree(mbr, level, pathBuffer);

	assert(n->m_level == level);

	if (n.get() == root.get())
	{
		assert(root.unique());
		root.relinquish();
	}
	n->insertData(dataLength, pData, mbr, id, pathBuffer, overflowTable);
}

bool SpatialIndex::xRTree::xRTree::deleteData_impl(const xMBR& mbr, id_type id)
{
	assert(mbr.m_dimension == m_dimension);

	std::stack<id_type> pathBuffer;
	NodePtr root = readNode(m_rootID);
	NodePtr l = root->findLeaf(mbr, id, pathBuffer);
	if (l.get() == root.get())
	{
		assert(root.unique());
		root.relinquish();
	}

	if (l.get() != 0)
	{
		Leaf* pL = static_cast<Leaf*>(l.get());
		pL->deleteData(id, pathBuffer);
		--(m_stats.m_u64Data);
		return true;
	}

	return false;
}

SpatialIndex::id_type SpatialIndex::xRTree::xRTree::writeNode(Node* n)
{
	uint8_t* buffer;
	uint32_t dataLength;
	n->storeToByteArray(&buffer, dataLength);

	id_type page;
	if (n->m_identifier < 0) page = StorageManager::NewPage;
	else page = n->m_identifier;

	try
	{
		m_pStorageManager->storeByteArray(page, dataLength, buffer);
		delete[] buffer;
	}
	catch (InvalidPageException& e)
	{
		delete[] buffer;
		std::cerr << e.what() << std::endl;
		throw;
	}

	if (n->m_identifier < 0)
	{
		n->m_identifier = page;
		++(m_stats.m_u32Nodes);

#ifndef NDEBUG
		try
		{
			m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel.at(n->m_level) + 1;
		}
		catch(...)
		{
			throw Tools::IllegalStateException("writeNode: writing past the end of m_nodesInLevel.");
		}
#else
		m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] + 1;
#endif
	}

	++(m_stats.m_u64Writes);

	for (size_t cIndex = 0; cIndex < m_writeNodeCommands.size(); ++cIndex)
	{
		m_writeNodeCommands[cIndex]->execute(*n);
	}
	return page;
}

SpatialIndex::xRTree::NodePtr SpatialIndex::xRTree::xRTree::readNode(id_type page)
{
    m_ts->m_indexIO++;
	uint32_t dataLength;
	uint8_t* buffer;

	try
	{
		m_pStorageManager->loadByteArray(page, dataLength, &buffer);
	}
	catch (InvalidPageException& e)
	{
		std::cerr << e.what() << std::endl;
		throw;
	}

	try
	{
		uint32_t level;
		memcpy(&level, buffer, sizeof(uint32_t));

		NodePtr n;

		if (level>0) n = m_indexPool.acquire();
		else if (level==0) n = m_leafPool.acquire();
		else throw Tools::IllegalStateException("readNode: failed reading the correct node type information");

		if (n.get() == nullptr)
		{
			if (level>0) n = NodePtr(new Index(this, -1, 0), &m_indexPool);
			else if (level == 0) n = NodePtr(new Leaf(this, -1), &m_leafPool);
		}

//		n->m_pTree = this;
		n->m_identifier = page;
		n->loadFromByteArray(buffer);

		++(m_stats.m_u64Reads);

		for (size_t cIndex = 0; cIndex < m_readNodeCommands.size(); ++cIndex)
		{
			m_readNodeCommands[cIndex]->execute(*n);
		}

		delete[] buffer;
		return n;
	}
	catch (...)
	{
		delete[] buffer;
		throw;
	}
}

void SpatialIndex::xRTree::xRTree::deleteNode(Node* n)
{
	try
	{
		m_pStorageManager->deleteByteArray(n->m_identifier);
	}
	catch (InvalidPageException& e)
	{
		std::cerr << e.what() << std::endl;
		throw;
	}

	--(m_stats.m_u32Nodes);
	m_stats.m_nodesInLevel[n->m_level] = m_stats.m_nodesInLevel[n->m_level] - 1;

	for (size_t cIndex = 0; cIndex < m_deleteNodeCommands.size(); ++cIndex)
	{
		m_deleteNodeCommands[cIndex]->execute(*n);
	}
}

void SpatialIndex::xRTree::xRTree::rangeQuery(RangeQueryType type, const IShape& query, IVisitor& v)
{
#ifdef HAVE_PTHREAD_H
    Tools::LockGuard lock(&m_lock);
#endif
    const xMBR *querybr= dynamic_cast<const xMBR*>(&query);
    const Cylinder *querycy= dynamic_cast<const Cylinder*>(&query);
    std::set<id_type > results;
    std::multimap<id_type, storeEntry> pending;
    storeEntry storee;
    if(querybr!= nullptr) {
        /* deprecated */
        bool isSlice;
        if (querybr->m_pLow[querybr->m_dimension - 1] == querybr->m_pHigh[querybr->m_dimension - 1])
            isSlice = true;
        else
            isSlice = false;
        std::stack<NodePtr> st;
        NodePtr root = readNode(m_rootID);

        if (root->m_children > 0 && root->m_nodeMBR.intersectsShape(query)) st.push(root);
        while (!st.empty()) {

            NodePtr n = st.top();
            st.pop();
            if (n->m_level == 0) {
                v.visitNode(*n);
                for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                    bool b;
                    if(results.count(m_ts->getTrajId(n->m_pIdentifier[cChild]))>0) continue;
                    if (m_bUsingMBR == true) {
                        if (type == ContainmentQuery) b = n->m_ptrMBR[cChild]->containsShape(query);
                        else b = n->m_ptrMBR[cChild]->intersectsShape(query);
                    } else {
                        if (type == ContainmentQuery) b = n->m_ptrxMBC[cChild]->containsShape(query);
                        else b = n->m_ptrxMBC[cChild]->intersectsShape(query);
                    }
                    if (b) {
                        sb+=1;
                        simpleData data = simpleData(n->m_pIdentifier[cChild], 0);
                        ++(m_stats.m_u64QueryResults);
                        if (m_DataType == TrajectoryType) {
                            //check if the timed slice is included in query
                            xMBR spatialbr(querybr->m_pLow, querybr->m_pHigh, 2);
                            xMBR timedbr;
                            if (isSlice) {
                                if (m_bUsingMBR)
                                    timedbr = xMBR((n->m_ptrMBR[cChild])->m_pLow, (n->m_ptrMBR[cChild])->m_pHigh,
                                                     (n->m_ptrMBR[cChild])->m_dimension - 1);
                                else
                                    (n->m_ptrxMBC[cChild])->getMBRAtTime(querybr->m_pLow[2], timedbr);
                                if (spatialbr.containsxMBR(timedbr)) {
                                    sbb+=1;
                                    m_stats.m_doubleExactQueryResults += 1;
                                    results.insert(m_ts->getTrajId(data.m_id));
                                    v.visitData(data);
                                } else {
                                        m_ts->m_trajIO += std::ceil(n->m_dataLen[cChild] / 4096.0);
                                        uint32_t len = n->m_dataLen[cChild] + n->m_pageOff[cChild];
                                        uint8_t *load = new uint8_t[len];
                                        m_ts->loadByteArray(n->m_pageNum[cChild], len, &load);
                                        uint8_t *ldata = load + n->m_pageOff[cChild];
                                        Trajectory partTraj;
                                        partTraj.loadFromByteArray(ldata);
                                        delete[](load);
                                        if (partTraj.intersectsxMBR(*querybr)) {
                                            m_stats.m_doubleExactQueryResults += 1;
                                            results.insert(m_ts->getTrajId(data.m_id));
                                            v.visitData(data);
                                        }
                                }
                            } else {
                                //time-period range query
                                bool bb;
                                if (m_bUsingMBR) {
                                    timedbr = xMBR((n->m_ptrMBR[cChild])->m_pLow, (n->m_ptrMBR[cChild])->m_pHigh,
                                                     (n->m_ptrMBR[cChild])->m_dimension - 1);
                                    bb = spatialbr.containsxMBR(timedbr);
                                } else {
                                    bb = n->m_ptrxMBC[cChild]->prevalidate(*querybr);
                                }
                                if (bb) {
                                    sbb+=1;
                                    m_stats.m_doubleExactQueryResults += 1;
                                    results.insert(m_ts->getTrajId(data.m_id));
                                    v.visitData(data);
                                } else {
                                        m_ts->m_trajIO += std::ceil(n->m_dataLen[cChild] / 4096.0);
                                        uint32_t len = n->m_dataLen[cChild] + n->m_pageOff[cChild];
                                        uint8_t *load = new uint8_t[len];
                                        m_ts->loadByteArray(n->m_pageNum[cChild], len, &load);
                                        uint8_t *ldata = load + n->m_pageOff[cChild];
                                        Trajectory partTraj;
                                        partTraj.loadFromByteArray(ldata);
                                        delete[](load);
                                        if (partTraj.intersectsxMBR(*querybr)) {
                                            m_stats.m_doubleExactQueryResults += 1;
                                            results.insert(m_ts->getTrajId(data.m_id));
                                            v.visitData(data);
                                        }
                                }
                            }
                        } else {
                            v.visitData(data);
                        }
                    }
                }
            } else {
                v.visitNode(*n);
                for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                    if (query.intersectsShape(*(n->m_ptrMBR[cChild]))) {
                        st.push(readNode(n->m_pIdentifier[cChild]));
                    }
                }
            }
        }
    }
    else if(querycy!= nullptr){
        bool isSlice=(querycy->m_startTime==querycy->m_endTime);
        std::stack<NodePtr> st;
        NodePtr root = readNode(m_rootID);
        if (root->m_children > 0 && root->m_nodeMBR.intersectsShape(query)) st.push(root);
        while (!st.empty()) {
            NodePtr n = st.top();
            st.pop();
            if (n->m_level == 0) {
                v.visitNode(*n);
                for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                    if(results.count(m_ts->getTrajId(n->m_pIdentifier[cChild]))>0) continue;
                    int b;
                    if (m_bUsingMBR) {
                        b = querycy->checkRel(*(n->m_ptrMBR[cChild]));
                    } else {
                        b = querycy->checkRel(*(n->m_ptrxMBC[cChild]));
                    }
                    if (b>0) {
                        sb += 1;
                        simpleData data = simpleData(n->m_pIdentifier[cChild], 0);
                        if (b==2) {
                            sbb += 1;
                            m_stats.m_doubleExactQueryResults += 1;
                            results.insert(m_ts->getTrajId(data.m_id));
                            pending.erase(m_ts->getTrajId(data.m_id));
                            v.visitData(data);
                            ++(m_stats.m_u64QueryResults);
                        } else {
                            storee.m_page = n->m_pageNum[cChild];
                            storee.m_off = n->m_pageOff[cChild];
                            storee.m_len = n->m_dataLen[cChild];
                            id_type id = m_ts->getTrajId(data.m_id);
                            if(bUsingSBBD) {
                                pending.insert(make_pair(id, storee));
                            }else{
                                uint32_t len = storee.m_off + storee.m_len;
                                m_ts->m_trajIO += std::ceil(len / 4096.0);
                                uint8_t *load;
                                m_ts->loadByteArray(storee.m_page, len, &load);
                                uint8_t *ldata = load + storee.m_off;
                                Trajectory partTraj;
                                partTraj.loadFromByteArray(ldata);
                                delete[](load);
                                if (partTraj.intersectsCylinder(*querycy)) {
                                    m_stats.m_doubleExactQueryResults += 1;
                                    results.insert(id);
                                    simpleData data = simpleData(id, 0);
                                    v.visitData(data);
                                    ++(m_stats.m_u64QueryResults);
                                }
                            }
                        }
                    }
                }
            } else {
                v.visitNode(*n);
                for (uint32_t cChild = 0; cChild < n->m_children; ++cChild) {
                    if (query.intersectsShape(*(n->m_ptrMBR[cChild]))) {
                        st.push(readNode(n->m_pIdentifier[cChild]));
                    }
                }
            }
        }
        if(bUsingSBBD) {
            while (!pending.empty()) {
                auto iter = pending.begin();
                storee = iter->second;
                id_type id = iter->first;
                uint32_t len = storee.m_off + storee.m_len;
                m_ts->m_trajIO += std::ceil(len / 4096.0);
                uint8_t *load;
                m_ts->loadByteArray(storee.m_page, len, &load);
                uint8_t *ldata = load + storee.m_off;
                Trajectory partTraj;
                partTraj.loadFromByteArray(ldata);
                delete[](load);
                if (partTraj.intersectsCylinder(*querycy)) {
                    m_stats.m_doubleExactQueryResults += 1;
                    results.insert(id);
                    pending.erase(id);
                    simpleData data = simpleData(id, 0);
                    v.visitData(data);
                    ++(m_stats.m_u64QueryResults);
                } else {
                    pending.erase(iter);
                }
            }
        }
    }
}

void SpatialIndex::xRTree::xRTree::selfJoinQuery(id_type id1, id_type id2, const xMBR& r, IVisitor& vis)
{
	throw;
}

void SpatialIndex::xRTree::xRTree::visitSubTree(NodePtr subTree, IVisitor& v)
{
	std::stack<NodePtr> st;
	st.push(subTree);

	while (! st.empty())
	{
		NodePtr n = st.top(); st.pop();
		v.visitNode(*n);

		if(n->m_level == 0)
		{
			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				Data data = Data(0, 0, *(n->m_ptrxMBC[cChild]), n->m_pIdentifier[cChild]);
				v.visitData(data);
				++(m_stats.m_u64QueryResults);
			}
		}
		else
		{
			for (uint32_t cChild = 0; cChild < n->m_children; ++cChild)
			{
				st.push(readNode(n->m_pIdentifier[cChild]));
			}
		}
	}
}

std::ostream& SpatialIndex::xRTree::operator<<(std::ostream& os, const xRTree& t)
{
	os	<< "Dimension: " << t.m_dimension << std::endl
		<< "Fill factor: " << t.m_fillFactor << std::endl
		<< "Index capacity: " << t.m_indexCapacity << std::endl
		<< "Leaf capacity: " << t.m_leafCapacity << std::endl
		<< "Tight MBRs: " << ((t.m_bTightMBRs) ? "enabled" : "disabled") << std::endl;

	if (t.m_treeVariant == RV_RSTAR)
	{
		os	<< "Near minimum overlap factor: " << t.m_nearMinimumOverlapFactor << std::endl
			<< "Reinsert factor: " << t.m_reinsertFactor << std::endl
			<< "Split distribution factor: " << t.m_splitDistributionFactor << std::endl;
	}

	if (t.m_stats.getNumberOfNodesInLevel(0) > 0)
		os	<< "Utilization: " << 100 * t.m_stats.getNumberOfData() / (t.m_stats.getNumberOfNodesInLevel(0) * t.m_leafCapacity) << "%" << std::endl
			<< t.m_stats;

	#ifndef NDEBUG
	os	<< "Leaf pool hits: " << t.m_leafPool.m_hits << std::endl
		<< "Leaf pool misses: " << t.m_leafPool.m_misses << std::endl
		<< "Index pool hits: " << t.m_indexPool.m_hits << std::endl
		<< "Index pool misses: " << t.m_indexPool.m_misses << std::endl
		<< "xMBR pool hits: " << t.m_xMBRPool.m_hits << std::endl
		<< "xMBR pool misses: " << t.m_xMBRPool.m_misses << std::endl
        << "xPoint pool hits: " << t.m_xPointPool.m_hits << std::endl
        << "xPoint pool misses: " << t.m_xPointPool.m_misses << std::endl;
#endif
    return os;
}
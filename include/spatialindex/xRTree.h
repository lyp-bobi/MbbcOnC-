/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2004, Marios Hadjieleftheriou
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

namespace SpatialIndex
{
	namespace xRTreeNsp
	{
#define PageSizeDefault 4096
		SIDX_DLL enum xRTreeVariant
		{
			RV_LINEAR = 0x0,
			RV_QUADRATIC,
			RV_RSTAR
		};

		SIDX_DLL enum BulkLoadMethod
		{
			BLM_STR = 0x0
		};

		SIDX_DLL enum PersistenObjectIdentifier
		{
			PersistentIndex = 0x1,
			PersistentLeaf = 0x2
		};

		SIDX_DLL enum RangeQueryType
		{
			ContainmentQuery = 0x1,
			IntersectionQuery = 0x2
		};
        class xRTree;
        xRTree* createNewxRTree(IStorageManager* store,long indexfan, long leaffan);


        SIDX_DLL xRTree* buildMBRRTreeWP(IStorageManager* store, const CUTFUNC &f);
        SIDX_DLL xRTree* buildMBCRTreeWP(IStorageManager* store, const CUTFUNC &f);
        SIDX_DLL xRTree* buildTBTreeWP(IStorageManager* store);
        SIDX_DLL xRTree* buildSTRTreeWP(IStorageManager* store);

        SIDX_DLL xRTree* buildMBRRTreeWoP(IStorageManager* store, const CUTFUNC &f);
        SIDX_DLL xRTree* buildMBCRTreeWoP(IStorageManager* store, const CUTFUNC &f);
        SIDX_DLL xRTree* buildTBTreeWoP(IStorageManager* store);
        SIDX_DLL xRTree* buildSTRTreeWoP(IStorageManager* store);
//        xRTree* buildMBCRTreeWP(xStore* store,function<void(xTrajectory&,list<xSBB>&)> cut);
//        xRTree* buildTBTree2(xStore* store,function<void(xTrajectory&,list<xSBB>&)> cut);
//        xRTree* buildTBTree(xStore* store);
//        xRTree* build3DRTree(xStore* store);
//        xRTree* buildSTRTree(xStore* store);
	}
}

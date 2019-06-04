//
// Created by chuang on 4/8/19.
//


#pragma once


namespace SpatialIndex
{
    namespace LERTree
    {
        class Leaf : public Node
        {
        public:
            virtual ~Leaf();

        protected:
            Leaf(LERTree* pTree, id_type id);

            //virtual NodePtr Node.hchooseSubtree(const Region& mbc, uint32_t level, std::stack<id_type>& pathBuffer);
//            virtual NodePtr findLeaf(const Region& mbc, id_type id, std::stack<id_type>& pathBuffer);

            //virtual void split(uint32_t dataLength, uint8_t* pData, Region& mbc, id_type id, NodePtr& left, NodePtr& right);

            //virtual void deleteData(id_type id, std::stack<id_type>& pathBuffer);

            friend class LERTree;
            friend class BulkLoader;
        }; // Leaf
    }
}

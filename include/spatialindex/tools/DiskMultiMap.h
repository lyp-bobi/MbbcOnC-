//
// Created by Chuang on 2020/3/1.
//

#ifndef SPATIALINDEX_DISKMULTIMAP_H
#define SPATIALINDEX_DISKMULTIMAP_H

#include <list>
#include <vector>
#include <cstring>
#include <algorithm>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <type_traits>
#include <cstdint> // if Offset is int32_t instead of ios::streamoff
using namespace std;

template<typename T> struct False : false_type {};
template<typename T>
struct is_string {
    enum { value = is_same<string, typename remove_cv<T>::type>::value };
};

SIDX_DLL class BinaryFile {
public:
    typedef int32_t Offset;

    ~BinaryFile() {
        m_stream.close();
    }

    bool openExisting(const std::string& filename) {
        if (m_stream.is_open())
            return false;
        m_stream.open(filename, ios::in|ios::out|ios::binary);
        return m_stream.good();
    }

    bool createNew(const std::string& filename) {
        if (m_stream.is_open())
            return false;
        m_stream.open(filename, ios::in|ios::out|ios::binary|ios::trunc);
        return m_stream.good();
    }

    void close() {
        if (m_stream.is_open())
            m_stream.close();
    }

    template<typename T>
    bool write(const T& data, Offset toOffset) {
        static_assert(!is_pointer<T>::value && !is_member_pointer<T>::value,
                      "BinaryFile::write can not be used to write a pointer");
        static_assert(!is_string<T>::value,
                      "BinaryFile::write can not be used to write a std::string");
        static_assert(is_trivially_copyable<T>::value ||
                      is_string<T>::value,  // suppress msg for string
                      "BinaryFile::write can not be used to write a non-trivially copyable class");

        return write(reinterpret_cast<const char*>(&data), sizeof(data), toOffset);
    }

    bool write(const char* data, size_t length, Offset toOffset) {
        return m_stream.seekp(toOffset, ios::beg)  &&
               m_stream.write(data, length);
    }

    template<typename T>
    bool read(T& data, Offset fromOffset) {
        static_assert(!is_pointer<T>::value && !is_member_pointer<T>::value,
                      "BinaryFile::read can not be used to read a pointer");
        static_assert(!is_string<T>::value,
                      "BinaryFile::read can not be used to read a std::string");
        static_assert(is_trivially_copyable<T>::value ||
                      is_string<T>::value,  // suppress msg for string
                      "BinaryFile::read can not be used to read a non-trivially copyable class");

        return read(reinterpret_cast<char*>(&data), sizeof(data), fromOffset);
    }

    bool read(char* data, size_t length, Offset fromOffset) {
        return m_stream.seekg(fromOffset, ios::beg)  &&
               m_stream.read(data, length);
    }

    template<typename T>
    bool read(const T&, Offset) {
        static_assert(False<T>::value,
                      "BinaryFile::read can not be used to read into a const or a temporary object");
        return false;
    }

    template<typename T>
    bool read(T*&&, Offset) {
        static_assert(False<T>::value,
                      "BinaryFile::read can not be used to read a pointer");
        return false;
    }

    Offset fileLength() {
        if (!m_stream.is_open())
            return -1;
        ios::streamoff currPos = m_stream.tellg();
        m_stream.seekg(0, ios::end);
        ios::streamoff length = m_stream.tellg();
        m_stream.seekg(currPos, ios::beg);
        return static_cast<Offset>(length);
    }

    bool isOpen() const {
        return m_stream.is_open();
    }

private:
    fstream m_stream; // fstreams are not copyable, so BinaryFiles won't be copyable.
};

typedef BinaryFile::Offset Location;
SIDX_DLL class DiskMultiMap {
public:
    struct MultiMapTuple {
        std::string key, value, context;
    };

    class Iterator {
    public:
        Iterator(){}
        Iterator(const std::list<MultiMapTuple>& nodeList): nodes(nodeList){}
        bool isValid() const{
            return nodes.size();
        }
        Iterator& operator++(){
            if (isValid()) nodes.pop_front();
            return *this;
        }
        MultiMapTuple operator*(){
            return isValid() ? nodes.front() : MultiMapTuple();
        }
    private:
        std::list<MultiMapTuple> nodes;
    };

    DiskMultiMap(){}
    ~DiskMultiMap(){close();}
    bool createNew(const std::string& filename, unsigned int numBuckets){
        bf.close();
        if (!bf.createNew(filename))
            return false;

        header.hashTableStart = sizeof(DiskMultiMapHeader);
        header.nodeDataStart = header.hashTableStart + numBuckets * sizeof(DiskMultiMapBucket);
        header.numBuckets = numBuckets;
        syncHeader();

        DiskMultiMapBucket b;
        for (int i = 0; i < numBuckets; i++) {
            b.location = bf.fileLength();
            if (!bf.write(b, b.location))
                return false;
        }

        return true;
    }
    bool openExisting(const std::string& filename){
        bf.close();
        return bf.openExisting(filename) && bf.read(header, 0);
    }
    void close(){
        bf.close();
    }
    bool insert(const std::string& key, const std::string& value, const std::string& context){
        BucketNode node;
        if (createNode(key, value, context, node)) {
            DiskMultiMapBucket bucket = bucketAt(hash(key));
            node.location = nextLocationToAddAt();
            node.nextNode = LIST_END;

            if (bucket.numNodes++ == 0) {
                bucket.firstNode = node.location;
                bucket.lastNode = node.location;
            } else {
                BucketNode lastNode = nodeAt(bucket.lastNode);
                lastNode.nextNode = node.location;
                bucket.lastNode = node.location;
                updateNode(lastNode);
            }

            updateNode(node);
            updateBucket(bucket);
            return true;
        }

        return false;
    }
    Iterator search(const std::string& key){
        DiskMultiMapBucket bucket = bucketAt(hash(key));
        if (bucket.numNodes == 0)
            return Iterator();

        std::list<BucketNode> bucketNodes(1, nodeAt(bucket.firstNode));
        for (int i = 1; i < bucket.numNodes; i++)
            bucketNodes.push_back(nodeAt(bucketNodes.back().nextNode));

        std::list<MultiMapTuple> nodeList;
        std::for_each(bucketNodes.begin(), bucketNodes.end(), [this, &nodeList, &key](const BucketNode& node) { if (!strcmp(node.key, key.c_str())) nodeList.push_back(convertToTuple(node)); });
        return Iterator(nodeList);
    }

    int erase(const std::string& key, const std::string& value, const std::string& context){
        DiskMultiMapBucket bucket = bucketAt(hash(key));
        if (bucket.numNodes == 0) return 0;

        int itemsErased = 0;
        std::vector<BucketNode> nodes(1, nodeAt(bucket.firstNode));
        for (int i = 1; i < bucket.numNodes; i++)
            nodes.push_back(nodeAt(nodes.back().nextNode));

        for (std::vector<BucketNode>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
            if (nodeEquals(*it, key, value, context)) {
                addDeletedNode(*it);
                it = --nodes.erase(it);
                itemsErased++;
            }
        }

        if (nodes.size() > 0) {
            nodes.back().nextNode = LIST_END;
            updateNode(nodes.back());
            for (int i = 0; i < nodes.size() - 1; i++) {
                nodes[i].nextNode = nodes[i + 1].location;
                updateNode(nodes[i]);
            }
        }

        bucket.firstNode = (nodes.size() ? nodes.front().location : LIST_END);
        bucket.lastNode = (nodes.size() ? nodes.back().location : LIST_END);
        bucket.numNodes = nodes.size();

        syncHeader();
        updateBucket(bucket);
        return itemsErased;
    }
private:
    static const Location LIST_END = 2147483647;
    struct DiskMultiMapHeader {
        unsigned int numBuckets = 0;
        long numDeletedNodes = 0;
        Location hashTableStart = 0;
        Location nodeDataStart = 0;
        Location lastDeletedNode = LIST_END;
        Location firstDeletedNode = LIST_END;
    };

    struct DiskMultiMapBucket {
        long numNodes = 0;
        Location location = 0;
        Location lastNode = LIST_END;
        Location firstNode = LIST_END;
    };

    struct BucketNode {
        char key[121];
        char value[121];
        char context[121];
        Location location;
        Location nextNode;
    };

    bool createNode(const std::string& key, const std::string& value, const std::string& context, BucketNode& node) const{
        if (strlen(key.c_str()) <= 120 && strlen(value.c_str()) <= 120 && strlen(context.c_str()) <= 120) {
            strcpy(node.key, key.c_str());
            strcpy(node.value, value.c_str());
            strcpy(node.context, context.c_str());
            return true;
        }

        return false;
    }
    bool nodeEquals(const BucketNode& node, const std::string& key, const std::string& value, const std::string& context) const{
        return strcmp(node.key, key.c_str()) == 0 && strcmp(node.value, value.c_str()) == 0 && strcmp(node.context, context.c_str()) == 0;
    }

    MultiMapTuple convertToTuple(const BucketNode& node) const{
        MultiMapTuple tuple;
        tuple.key = node.key;
        tuple.value = node.value;
        tuple.context = node.context;
        return tuple;
    }

    unsigned long hash(const std::string& key) const{
        unsigned long hash = 5381;	// djb2 hash algorithm
        for (const char& c : key)
            hash = ((hash << 5) + hash) + static_cast<int>(c);
        return hash % header.numBuckets;
    }
    Location nextLocationToAddAt(){
        if (header.numDeletedNodes > 0) {
            BucketNode deletedNode = nodeAt(header.firstDeletedNode);
            header.firstDeletedNode = deletedNode.nextNode;
            if (--header.numDeletedNodes == 0)
                header.lastDeletedNode = LIST_END;
            syncHeader();
            return deletedNode.location;
        }

        return bf.fileLength();
    }
    void addDeletedNode(BucketNode node){
        if (header.numDeletedNodes++ > 0) {
            BucketNode lastNode = nodeAt(header.lastDeletedNode);
            lastNode.nextNode = node.location;
            updateNode(lastNode);
        } else header.firstDeletedNode = node.location;

        header.lastDeletedNode = node.location;
        node.nextNode = LIST_END;
        updateNode(node);
        syncHeader();
    }

    DiskMultiMapBucket bucketAt(unsigned long index){
        DiskMultiMapBucket bucket;
        bf.read(bucket, header.hashTableStart + index * sizeof(DiskMultiMapBucket));
        return bucket;
    }
    BucketNode nodeAt(Location offset){
        BucketNode node;
        bf.read(node, offset);
        return node;
    }
    bool updateBucket(DiskMultiMapBucket bucket){ return bf.write(bucket, bucket.location); }
    bool updateNode(BucketNode node){ return bf.write(node, node.location); }
    bool syncHeader(){ return bf.write(header, 0); }

    DiskMultiMapHeader header;
    BinaryFile bf;

};

#endif //SPATIALINDEX_DISKMULTIMAP_H

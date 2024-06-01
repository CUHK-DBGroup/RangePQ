#ifndef TESTTREE_WBTCHUNK_HPP
#define TESTTREE_WBTCHUNK_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

const int MAX_INT=2e9;

class WBTChunk {
    struct Element{
        int key,id,label;
        Element(int _key, int _id, int _label):key(_key),id(_id),label(_label){};
        bool operator<(const Element &b) const{
            if(key!=b.key) return key<b.key;
            else return id<b.id;
        }

    };

    struct Chunk{
        std::pair<int,int> L,R;
        int size;
        std::unordered_map<int, std::vector<Element> > ivfIndex;

        Chunk(const std::vector<Element> &elements){
            size = elements.size();
            L = R = {elements[0].key,elements[0].id};
            for(auto e:elements)
            {
                L = min(L,{e.key,e.id});
                R = max(R, {e.key,e.id});
                int label = e.label;
                if(ivfIndex.count(label) == 0) ivfIndex[label] = std::vector<Element>();
                ivfIndex[label].push_back(e);
            }
            for(auto &vec :ivfIndex) sort(vec.second.begin(),vec.second.end());
        }

        void insertElement(int key, int id, int label){
            Element element(key,id,label);
            auto &elements = ivfIndex[label];
            elements.insert(upper_bound(elements.begin(),elements.end(), element),element);
            size++;
            L = std::min(L,{key,id});
            R = std::max(R,{key,id});
        }

        Chunk *splitChunk(int chunkSize){
            std::vector<Element> newElements, tmpElements;

            for(const auto& vec:ivfIndex){
                for(auto e:vec.second){
                    tmpElements.push_back(e);
                }
            }

            sort(tmpElements.begin(),tmpElements.end());

            auto th = tmpElements[chunkSize - 1];
            R = {th.key, th.id};

            std::vector<int> del;
            for(auto &vec:ivfIndex){
                auto &List = vec.second;
                while(List.size() > 0 && th < List[List.size()-1]){
                    newElements.push_back(List[List.size()-1]);
                    List.pop_back();
                }
                if(List.size()==0) del.push_back(vec.first);
            }
            for(int id:del)
                ivfIndex.erase(id);
            size = chunkSize;
            return new Chunk(newElements);
        }

        void erase(int key, int id, int label){
            // if(key == 7548) std::cout<<id<<" "<<label<<std::endl;
            auto &elements = ivfIndex[label];
            for(auto it = elements.begin();it!=elements.end();it++)
                if((*it).key==key && (*it).id == id) {
                    size--;
                    elements.erase(it);
                    break;
                }
            if(ivfIndex[label].size()==0) ivfIndex.erase(label);
        }


        void getAnsInChunk(std::unordered_set<int> &ans, int L, int R){
            Element searchElement(L,0,0);
            for(auto pr:ivfIndex) {
                auto &elements = pr.second;
                // printf("%d:",pr.first);
                // for(auto x:elements) printf("%d:%d ",x.id,x.key);
                // std::cout<<std::endl;
                auto it = std::lower_bound(elements.begin(), elements.end(), searchElement);
                if(it!=elements.end()) {
                    if (it->key <= R) {
                        ans.insert(it->label);
                    }
                }
            }
        }

        void getIdInChunk(std::vector<int> &ans, int L, int R, int label){
            if(ivfIndex.count(label) == 0) return;
            auto &elements = ivfIndex[label];
            Element searchElement(L,0,0);
            auto it = std::lower_bound(elements.begin(), elements.end(), searchElement);
            for(;it != elements.end(); ++it){
                if(it->key <= R){
                    ans.push_back(it->id);
                }
                else break;
            }
        }

        void getIvf(std::unordered_map<int,int> &ivfSet){
            ivfSet.clear();
            for(auto p:ivfIndex){
                ivfSet.insert({p.first,p.second.size()});
            }
        }
    };

    struct WBTNode{
        std::pair<int,int> L,R;
        Chunk *chunk;
        int size;
        WBTNode *son[2];
        std::unordered_map<int,int> ivfSet;

        WBTNode(const std::vector<Element> &elements){
            chunk = new Chunk(elements);
            size = 1;
            son[0] = son[1] = nullptr;
            L = chunk->L;
            R = chunk->R;
            chunk->getIvf(ivfSet);
        }

        WBTNode(Chunk *_chunk){
            chunk = _chunk;
            size = 1;
            son[0] = son[1] = nullptr;
            L = chunk->L;
            R = chunk->R;
            chunk->getIvf(ivfSet);
        }

        ~WBTNode(){
            delete chunk;
            ivfSet.clear();
        }
    };

private:
    double alpha = 0.29;
    int activeChunks = 0;
    WBTNode *root;
    std::vector<Element> elementsBuf;
    std::vector<Chunk*> chunksBuf;

    bool insertToChunk(Chunk *chunk, int key, int id, int label){
        chunk->insertElement(key, id, label);
        if(chunkSize == chunk->size * 2) return true; else return false;
    }

    bool eraseInChunk (Chunk *chunk, int key, int id, int label){
        chunk->erase(key, id, label);
        if(chunkSize/2 -1  == chunk->size) return true; else return false;
    }

    void addElementsToIvf(WBTNode *node, int label, int value){
        if(node->ivfSet.count(label) > 0) node->ivfSet[label] += value;
        else node->ivfSet[label] = value;
        // node->size += value;
    }

    void addElementToNode(WBTNode *node, int key,int id, int label){
        activeChunks += insertToChunk(node->chunk,key,id,label);

        node->L = std::min(node->L,{key,id});
        node->R = std::max(node->R,{key,id});
        addElementsToIvf(node, label, 1);
    }

    void eraseElementFromIvf(WBTNode *node, int label){
        node->ivfSet[label]--;
        if(node->ivfSet[label] == 0) node->ivfSet.erase(label);
    }

    void eraseElementFromNode(WBTNode *node, int key,int id, int label){
        eraseElementFromIvf(node, label);
        activeChunks-= eraseInChunk(node->chunk,key, id, label);
    }

    static WBTNode *newNode(Chunk *chunk){
        return new WBTNode(chunk);
    }

    void updateRange(WBTNode *node){
        node->L = node->chunk->L;
        node->R = node->chunk->R;
        if(node->son[0] != nullptr) node->L = node->son[0]->L;
        if(node->son[1] != nullptr) node->R = node->son[1]->R;
    }

    void updateSize(WBTNode *node){
        node->size = 1;
        if(node->son[0] != nullptr) node->size += node->son[0]->size;
        if(node->son[1] != nullptr) node->size += node->son[1]->size;
    }

    void update(WBTNode *node){
        node->chunk->getIvf(node->ivfSet);

        if(node->son[0] != nullptr) {
            for (auto it = node->son[0]->ivfSet.begin(); it != node->son[0]->ivfSet.end(); ++it) {
                addElementsToIvf(node, it->first, it->second);
            }
        }

        if(node->son[1] != nullptr) {
            for (auto it = node->son[1]->ivfSet.begin(); it != node->son[1]->ivfSet.end(); ++it) {
                addElementsToIvf(node, it->first, it->second);
            }
        }
        updateRange(node);
        updateSize(node);
    }

    WBTNode *rotate(WBTNode *k1, bool d) {  // d 表示将哪个儿子旋转到 x 的位置
        WBTNode *k2 = k1->son[d];
        k1->son[d] = k2->son[!d];
        k2->son[!d] = k1;
        update(k1);
        update(k2);
        return k2;
    }

    WBTNode *maintain(WBTNode *node) {
        int d;
        if (node->size == 1) return node;
        if (node->son[0] != nullptr && node->son[0]->size < node->size * alpha) d = 1;
        else if (node->son[1] != nullptr && node->son[1]->size < node->size * alpha) d = 0;
        else return node;
        if (node->son[d] == nullptr) return node;
        if (node->son[d]->son[!d] != nullptr &&node->son[d]->son[!d]->size * (1 - alpha) >= node->son[d]->size * (1 - 2 * alpha))
            node->son[d] =  rotate(node->son[d], d ^ 1);
        return rotate(node, d);
    }

    void insertNode(WBTNode* &node, WBTNode *newNode){
        if (node == nullptr) {
            node = newNode;
            return;
        }
        for(auto it = newNode->ivfSet.begin();it!=newNode->ivfSet.end();it++)
            addElementsToIvf(node,it->first,it->second);
        bool d = newNode->R > node->R;
        node->L = std::min(node->L, newNode->L);
        node->R = std::max(node->R, newNode->R);
        insertNode(node->son[d], newNode);
        node = maintain(node);
        updateSize(node);
    }

    void splitNode(WBTNode* &node){
        std::cout<<"split"<<std::endl;
        Chunk *newChunk = node->chunk->splitChunk(chunkSize);
        insertNode(node->son[1], new WBTNode(newChunk));
        updateSize(node);
    }

    void insertElement(WBTNode* &node, std::pair<int,int> keyId, int label) {
        bool d = keyId > node->chunk->R;
        // printf("%d %d %d\n",keyId.first,node->chunk->L.first, node->chunk->R.first);
        // std::cout<<std::endl;
        if((node->chunk->L <= keyId && keyId <= node->chunk->R)||(node->son[d] == nullptr)) {
            // printf("In:%d %d %d\n",keyId.first,node->chunk->L.first, node->chunk->R.first);
            // std::cout<<std::endl;
            addElementToNode(node, keyId.first,keyId.second,label);
            if(node->chunk->size == 2*chunkSize)
                splitNode(node);
            return;
        }
        addElementsToIvf(node, label, 1);
        node->L = std::min(node->L, keyId);
        node->R = std::max(node->R, keyId);
        insertElement(node->son[d], keyId ,label);
        node = maintain(node);
        updateSize(node);
    }

    void eraseElement(WBTNode* &node, std::pair<int,int> keyId, int label){
        if(node->chunk->L <= keyId && keyId <= node->chunk->R){
            eraseElementFromNode(node, keyId.first,keyId.second, label);
            return;
        }
        int d = keyId > node->chunk->R;
        eraseElementFromIvf(node, label);
        eraseElement(node->son[d], keyId, label);
    }

    void askChunk(std::unordered_set<int> &ans,WBTNode *node, std::pair<int,int> L, std::pair<int,int> R){
        if(L < node->L && node->R < R){
            for(auto it = node->ivfSet.begin(); it != node->ivfSet.end(); ++it){
                ans.insert(it->first);
            }
            return;
        }
        if(L < node->chunk->L && node->chunk->R < R)
            for(auto it = node->chunk->ivfIndex.begin(); it != node->chunk->ivfIndex.end(); ++it){
                ans.insert(it->first);
            }
        if(node->son[0] != nullptr && L < node->son[0]->R) askChunk(ans,node->son[0], L, R);
        if(node->son[1] != nullptr && R > node->son[1]->L) askChunk(ans,node->son[1], L, R);
    }

    Chunk *findChunk(WBTNode *node, std::pair<int,int> keyId){
        if(node->chunk->L <= keyId && keyId <= node->chunk->R)
            return node->chunk;
        int d = keyId > node->chunk->R;
        if(node->son[d] == nullptr) return nullptr;
        return findChunk(node->son[d], keyId);
    }

    WBTNode *build(int L, int R){
        if(L > R) return nullptr;
        int mid = (L + R) / 2;
        WBTNode *node = newNode(chunksBuf[mid]);
        node->son[0] = build(L, mid - 1);
        node->son[1] = build(mid + 1, R);
        update(node);
        return node;
    }

    void clearNode(WBTNode *node, bool rebuild){
        if(node == nullptr) return;
        clearNode(node->son[0],rebuild);
        if(rebuild)
            for(auto it = node->chunk->ivfIndex.begin(); it != node->chunk->ivfIndex.end(); ++it){
                auto vec = it->second;
                for(auto it2 = vec.begin(); it2 !=vec.end(); ++it2){
                    elementsBuf.push_back(*it2);
                }
            }
        clearNode(node->son[1],rebuild);
        delete node;
    }

    void getIdInTree(std::vector<int> &ans,WBTNode *node, std::pair<int,int> L, std::pair<int,int> R, int labelId, int num){
        if(L < node->chunk->L && node->chunk->R < R && node->chunk->ivfIndex.count(labelId)>0){
            for(auto it = node->chunk->ivfIndex[labelId].begin(); it != node->chunk->ivfIndex[labelId].end(); ++it){
                ans.push_back(it->id);
            }
        }
        if(node->son[0] != nullptr && L < node->son[0]->R
            && ans.size()<num && node->son[0]->ivfSet.count(labelId)>0)
            getIdInTree(ans,node->son[0], L, R, labelId, num);
        if(node->son[1] != nullptr && R > node->son[1]->L
            && ans.size()<num && node->son[1]->ivfSet.count(labelId)>0)
            getIdInTree(ans,node->son[1], L, R, labelId, num);
    }

public:
    int chunkSize;

    void insert(int key, int id, int label){
        insertElement(root, {key, id}, label);
    }

    void rebuildTree(){
        clearNode(root, true);
        buildTree();
    }

    void erase(int key, int id, int label){
        // std::cout<<"erase"<<key<<" "<<id<<" "<<std::flush;
        eraseElement(root, {key, id}, label);
        if(activeChunks * 2 < root->size){
            std::cout<<"rebuild"<<activeChunks<<" "<<root->size<<std::endl;
            rebuildTree();
        }
    }

    void ask(std::unordered_set<int> &ans,int L, int R){
        auto l = std::max(root->L, {L,0});
        auto r = std::min(root->R, {R,MAX_INT});
        if(l > r) return;
        askChunk(ans,root, l, r);
        auto LChunk = findChunk(root,l);
        auto RChunk = findChunk(root,r);
        if(LChunk == nullptr && RChunk == nullptr) return;
        if(LChunk == RChunk){
            // std::cout<<"111"<<std::endl;
            LChunk->getAnsInChunk(ans,L,R);
        }
        else{
            if(LChunk != nullptr) LChunk->getAnsInChunk(ans,L,R);
            if(RChunk != nullptr) RChunk->getAnsInChunk(ans,L,R);
        }
    }

    void getId(const std::vector<int> &labelSet, std::vector<int> &ans, int L, int R, int LNum){
        int i = 0;
        auto l = std::max(root->L, {L,0});
        auto r = std::min(root->R, {R,MAX_INT});
        if(l > r) return;
        auto LChunk = findChunk(root,l);
        auto RChunk = findChunk(root,r);
        if(LChunk == RChunk) RChunk = nullptr;
        while(ans.size()<LNum){
            getIdInTree(ans, root,l,r,labelSet[i], LNum);
            if(ans.size()<LNum && LChunk != nullptr) LChunk->getIdInChunk(ans, L, R, labelSet[i]);
            if(ans.size()<LNum && RChunk != nullptr) RChunk->getIdInChunk(ans, L, R, labelSet[i]);
            i++;
            if(i == labelSet.size()) break;
        }
        if(ans.size()>LNum) ans.resize(LNum);
    }

    void clear(){
        clearNode(root,false);
    }

    void buildTree(){
        std::vector<Element> elements;
        for(int i = 0; i < elementsBuf.size(); i++){
            elements.push_back(elementsBuf[i]);
            if(elements.size() == chunkSize){
                chunksBuf.push_back(new Chunk(elements));
                elements.clear();
                activeChunks ++;
            }
        }
        if(elements.size()>0) {
            chunksBuf.push_back(new Chunk(elements));
            elements.clear();
        }
        root = build(0, chunksBuf.size()-1);
        chunksBuf.clear();
    }

    WBTChunk(){
        root = nullptr;
    }

    WBTChunk(const std::vector<int> &keyList, const std::vector<int> &idList, const std::vector<int> &labelList, int _chunkSize){
        for(int i = 0; i < keyList.size(); i++){
            elementsBuf.emplace_back(keyList[i],idList[i], labelList[i]);
        }
        sort(elementsBuf.begin(),elementsBuf.end());
        // std::cout<<std::endl;
        chunkSize = _chunkSize;
        buildTree();
    }

    ~WBTChunk(){
        clear();
    }
};


#endif //TESTTREE_WBTCHUNK_HPP

#ifndef TESTTREE_WBTROTATE_HPP
#define TESTTREE_WBTROTATE_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

class WBTRotate {
    struct WBTNode{
        int id, label;
        int value, L, R;
        int size;
        WBTNode *son[2];
        std::unordered_map<int,int> ivfSet;
        WBTNode(int _key, int _value, int _label):id(_key), value(_value), L(_value), R(_value), label(_label){
            size = 1;
            son[0] = son[1] = nullptr;
            ivfSet[label] = 1;
        }

        ~WBTNode(){
            ivfSet.clear();
        }
    };
private:

    bool cmp(WBTNode *n1, WBTNode *n2){
        if(n1->value == n2->value) return n1->id > n2->id;
        else return n1->value > n2->value;
    }

    bool cmp(int id, int value, WBTNode *n2){
        if(value == n2->value) return id > n2->id;
        else return value > n2->value;
    }

    double alpha = 0.29;
    WBTNode *root;

    void addElementToNode(WBTNode *node, int key, int value){
        if(node->ivfSet.count(key) > 0) node->ivfSet[key] += value;
        else node->ivfSet[key] = value;
        node->size += value;
    }

    void eraseElementFromNode(WBTNode *node, int key){
        node->ivfSet[key]--;
        node->size--;
        if(node->ivfSet[key] == 0) node->ivfSet.erase(key);
    }

    WBTNode *newNode(int id, int value, int label){
        return new WBTNode(id, value, label);
    }

    void updateRange(WBTNode *node){
        node->L = node->R = node->value;
        if(node->son[0] != nullptr) node->L = node->son[0]->L;
        if(node->son[1] != nullptr) node->R = node->son[1]->R;
    }

    void update(WBTNode *node){
        node->ivfSet.clear();

        node->L = node->R = node->value;
        addElementToNode(node, node->label, 1);

        if(node->son[0] != nullptr) {
            for (auto it = node->son[0]->ivfSet.begin(); it != node->son[0]->ivfSet.end(); ++it) {
                addElementToNode(node, it->first, it->second);
            }
            node->L = node->son[0]->L;
        }

        if(node->son[1] != nullptr) {
            for (auto it = node->son[1]->ivfSet.begin(); it != node->son[1]->ivfSet.end(); ++it) {
                addElementToNode(node, it->first, it->second);
            }
            node->R = node->son[1]->R;
        }

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

    void insert(WBTNode* &node, WBTNode *nNode) {
        if (node == nullptr) {
            node = nNode;
            return;
        }
        addElementToNode(node, nNode->label, 1);
        node->L = std::min(node->L, nNode->value);
        node->R = std::max(node->R, nNode->value);
        bool d = cmp(nNode, node);
        insert(node->son[d], nNode);
        node = maintain(node);
    }

    void ask(std::unordered_set<int> &ans,WBTNode *node, int L, int R){
        if(L <= node->L && node->R <= R){
            for(auto it = node->ivfSet.begin(); it != node->ivfSet.end(); ++it){
                ans.insert(it->first);
            }
            return;
        }
        if(L <= node->value && node->value <= R)
            ans.insert(node->label);
        if(node->son[0] != nullptr && L <= node->son[0]->R) ask(ans,node->son[0], L, R);
        if(node->son[1] != nullptr && R >= node->son[1]->L) ask(ans,node->son[1], L, R);
    }

    WBTNode *build(int L, int R, const std::vector<std::pair<std::pair<int,int>, int > > &w){
        if(L > R) return nullptr;
        int mid = (L + R) / 2;
        WBTNode *node = newNode(w[mid].first.first, w[mid].first.second,w[mid].second);
        node->son[0] = build(L, mid - 1, w);
        node->son[1] = build(mid + 1, R, w);
        update(node);
        return node;
    }

    WBTNode* find(WBTNode *node, int id, int value){
        if(node == nullptr) return nullptr;
        if(node->value == value) return node;
        return find(node->son[cmp(id,value, node)], id, value);
    }

    WBTNode* findMax(WBTNode *node){
        if(node->son[1] == nullptr) return node;
        else return findMax(node->son[1]);
    }

    void eraseKey(WBTNode* &node, int id, int value, int label, bool eraseNode){
        eraseElementFromNode(node, label);
        int d = cmp(id,value, node);
        if(node->id == id){
            if(eraseNode) node = node->son[node->son[0] == nullptr];
            return;
        }
        eraseKey(node->son[d], id, value, label, eraseNode);
        updateRange(node);
        node = maintain(node);
    }

    void clearNode(WBTNode *node){
        if(node == nullptr) return;
        clearNode(node->son[0]);
        clearNode(node->son[1]);
        delete node;
    }

public:
    void insert(int id, int value, int label){
        WBTNode *nNode = newNode(id, value, label);
        insert(root, nNode);
    }

    void ask(std::unordered_set<int> &ans,int L, int R){
        L = std::max(root->L, L);
        R = std::min(root->R, R);
        if(L > R) return;
        ask(ans,root,L,R);
    }

    WBTRotate(){
        root = nullptr;
    }

    bool erase(int id, int value){
        WBTNode *delNode = find(root, id, value);
        if(delNode == nullptr) return false;
        int label = delNode->label;
        if(delNode->son[0] == nullptr || delNode->son[1] == nullptr){
            eraseKey(root, id, value, label, true);
            delete delNode;
        }
        else{
            WBTNode *swapNode = findMax(delNode->son[0]);
            eraseKey(delNode->son[0], swapNode->id,swapNode->value, swapNode->label, true);
            eraseKey(root, id, value, delNode->label, false);
            delNode->value = swapNode->value;
            delNode->label = swapNode->label;
            delNode->id = swapNode->id;
            updateRange(delNode);
            delete swapNode;
        }
        return true;
    }

    void clear(){
        clearNode(root);
    }

    WBTRotate(const std::vector<std::pair<std::pair<int,int>, int > > &w){
        auto W = w;
        std::sort(W.begin(),W.end());
        int n = W.size();
        root = build(0, n - 1, W);
    }

    ~WBTRotate(){
        clear();
    }
};


#endif //TESTTREE_WBTROTATE_HPP

#ifndef TESTTREE_BST_HPP
#define TESTTREE_BST_HPP

#include <unordered_set>
#include <vector>
#include <algorithm>

class BST {
    struct WBTNode{
        int id;
        int value, L, R;
        int size;
        WBTNode *son[2];
        WBTNode(int _key, int _value):id(_key), value(_value), L(_value), R(_value){
            size = 1;
            son[0] = son[1] = nullptr;
        }

        ~WBTNode(){
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

    WBTNode *newNode(int id, int value){
        return new WBTNode(id, value);
    }

    void update(WBTNode *node){
        node->L = node->R = node->value;
        node->size = 1;
        if(node->son[0] != nullptr) {
            node->L = node->son[0]->L;
            node->size += node->son[0]->size;
        }
        if(node->son[1] != nullptr) {
            node->R = node->son[1]->R;
            node->size += node->son[1]->size;
        }
    }


    WBTNode *rotate(WBTNode *k1, bool d) {
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
        bool d = cmp(nNode, node);
        insert(node->son[d], nNode);
        node = maintain(node);
        update(node);
    }

    void ask(std::unordered_set<size_t> &ans,WBTNode *node, int L, int R){
        if(L <= node->value && node->value <= R){
            ans.insert(node->id);
            // return;
        }
        if(node->son[0] != nullptr && L <= node->son[0]->R) ask(ans,node->son[0], L, R);
        if(node->son[1] != nullptr &&R >= node->son[1]->L) ask(ans,node->son[1], L, R);
    }

    WBTNode *build(int L, int R, const std::vector<std::pair<int,int> > &w){
        if(L > R) return nullptr;
        int mid = (L + R) / 2;
        WBTNode *node = newNode(w[mid].first, w[mid].second);
        node->son[0] = build(L, mid - 1, w);
        node->son[1] = build(mid + 1, R, w);
        update(node);
        return node;
    }

    WBTNode* find(WBTNode *node, int id, int value){
        if(node == nullptr) return nullptr;
        if(node->value == value && node->id == id) return node;
        return find(node->son[cmp(id,value, node)], id, value);
    }

    WBTNode* findMax(WBTNode *node){
        if(node->son[1] == nullptr) return node;
        else return findMax(node->son[1]);
    }

    void eraseKey(WBTNode* &node, int id, int value, bool eraseNode){
        int d = cmp(id,value, node);
        if(node->id == id){
            if(eraseNode) node = node->son[node->son[0] == nullptr];
            return;
        }
        eraseKey(node->son[d], id, value, eraseNode);
        update(node);
        // if(id==289)puts("Maintain begin");
        node = maintain(node);
        // if(id==289)puts("Maintain end");
    }

    void clearNode(WBTNode *node){
        if(node == nullptr) return;
        clearNode(node->son[0]);
        clearNode(node->son[1]);
        delete node;
    }

public:
    void insert(int id, int value){
        WBTNode *nNode = newNode(id, value);
        insert(root, nNode);
    }

    void ask(std::unordered_set<size_t> &ans,int L, int R){
        L = std::max(root->L, L);
        R = std::min(root->R, R);
        if(L > R) return;
        ask(ans,root,L,R);
    }

    BST(){
        root = nullptr;
    }

    bool erase(int id, int value){
        WBTNode *delNode = find(root, id, value);
        if(delNode == nullptr) return false;
        if(delNode->son[0] == nullptr || delNode->son[1] == nullptr){
            eraseKey(root, id, value, true);
            delete delNode;
        }
        else{
            WBTNode *swapNode = findMax(delNode->son[0]);
            eraseKey(delNode->son[0], swapNode->id,swapNode->value, true);
            eraseKey(root, id, value, false);
            delNode->value = swapNode->value;
            delNode->id = swapNode->id;
            update(delNode);
            delete swapNode;
        }
        return true;
    }

    void clear(){
        clearNode(root);
    }

    BST(const std::vector<std::pair<int,int> > &w){
        auto W = w;
        std::sort(W.begin(),W.end());
        int n = W.size();
        root = build(0, n - 1, W);
    }

    ~BST(){
        clear();
    }
};


#endif //TESTTREE_BST_HPP

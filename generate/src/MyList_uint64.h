#include <vector>
#include <numeric>
#include <random>
#include <utility>
#include <algorithm>
#include <set>
#include <map>
#include <iostream>
#include <cstring>

using std::vector;
using std::pair;
using std::set;
using std::map;
using std::string;

struct LIST_node_64 {
    u_int64_t val;
    LIST_node_64 *pre, *nxt;
    LIST_node_64() { pre = nxt = NULL; }
    LIST_node_64(u_int64_t x) {
        val = x;
        pre = nxt = NULL;
    }
};

class MyList_64 {

public:
    MyList_64() {
        head = tail = NULL, List_size = 0;
    }
    ~MyList_64() { }

    void insert(u_int64_t x) {
        LIST_node_64 *t = new LIST_node_64;
        t->val = x, t->pre = tail, t->nxt = NULL;
        if(head == NULL)
            head = tail = t;
        else
            tail->nxt = t, tail = t;
        List_size ++;
    }

    void delete_node(LIST_node_64 *t) {

        LIST_node_64 *pre = t->pre, *nxt = t->nxt;
        if(List_size == 1)
            head = tail = NULL;
        else if(t == head)
            head = t->nxt;
        else if(t == tail)
            pre->nxt = NULL, tail = pre;
        else    
            pre->nxt = nxt, nxt->pre = pre;

        t->pre = t->nxt = NULL;
        delete t;
        List_size --;
    }

    int size() {
        return List_size;
    }

    u_int64_t get_head_val() {
        return head->val;
    }

    void Free() {
        for(LIST_node_64 *p = head, *q; p != NULL; p = q) {
            q = p->nxt; delete p;
        }
        head = tail = NULL, List_size = 0;
    }

    LIST_node_64* get_pre(LIST_node_64 *t) {
        return t->pre;
    }
    LIST_node_64* get_nxt(LIST_node_64 *t) {
        return t->nxt;
    }

    LIST_node_64 *head, *tail;
private:
    int List_size;
};
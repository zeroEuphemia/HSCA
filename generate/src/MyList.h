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

struct LIST_node {
    int val;
    LIST_node *pre, *nxt;
    LIST_node() { pre = nxt = NULL; }
    LIST_node(int x) {
        val = x;
        pre = nxt = NULL;
    }
};

class MyList {

public:
    MyList() {
        head = tail = NULL, List_size = 0;
    }
    ~MyList() { }

    void insert(int x) {
        LIST_node *t = new LIST_node;
        t->val = x, t->pre = tail, t->nxt = NULL;
        if(head == NULL)
            head = tail = t;
        else
            tail->nxt = t, tail = t;
        List_size ++;
    }

    void delete_node(LIST_node *t) {

        LIST_node *pre = t->pre, *nxt = t->nxt;
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

    int get_head_val() {
        return head->val;
    }

    void Free() {
        for(LIST_node *p = head, *q; p != NULL; p = q) {
            q = p->nxt; delete p;
        }
        head = tail = NULL, List_size = 0;
    }

    LIST_node* get_pre(LIST_node *t) {
        return t->pre;
    }
    LIST_node* get_nxt(LIST_node *t) {
        return t->nxt;
    }

    LIST_node *head, *tail;
private:
    int List_size;
};
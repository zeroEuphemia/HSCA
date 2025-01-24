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

class MyBitSet {
public:
    MyBitSet(): siz(0), cap(0) {
        bitmap = (std::nullptr_t) nullptr;
    }
    MyBitSet(long long _siz): siz(_siz) {
        cap = (_siz + 63) >> 6;
        bitmap = new uint64_t[cap]{0};
    }
    MyBitSet(const MyBitSet& bs): siz(bs.siz), cap(bs.cap) {
        bitmap = new uint64_t[bs.cap]{0};
        memcpy(bitmap, bs.bitmap, sizeof(uint64_t) * cap);
    }
    ~MyBitSet(){  }

    void Free() {
        delete [] bitmap;
    }

    MyBitSet& operator = (const MyBitSet &bs) {
        if (this == &bs) return *this;
        delete [] bitmap;
        siz = bs.siz;
        cap = bs.cap;
        bitmap = new uint64_t[cap]{0};
        memcpy(bitmap, bs.bitmap, sizeof(uint64_t) * cap);
        return *this;
    }

    MyBitSet& operator |= (const MyBitSet &bs) {
        int len = cap < bs.cap ? cap : bs.cap;
        for(int i = 0; i < len; i ++)
            bitmap[i] |= bs.bitmap[i];
        return *this;
    }

    MyBitSet& operator %= (const MyBitSet &bs) { // |= ~bs
        int len = cap < bs.cap ? cap : bs.cap;
        for(int i = 0; i < len; i ++)
            bitmap[i] |= ~bs.bitmap[i];
        return *this;
    }

    bool set(int idx) {
        uint64_t& t1 = bitmap[idx >> 6];
        uint64_t cg = (1ull << (idx & 63));
        if ((t1 & cg) == 0ull){
            t1 |= cg;
            return true;
        }
        return false;
    }

    bool unset(int idx) {
        uint64_t& t1 = bitmap[idx >> 6];
        uint64_t cg = (1ull << (idx & 63));
        if ((t1 & cg) != 0ull){
            t1 ^= cg;
            return true;
        }
        return false;
    }

    bool get(int idx) const {
        return (bitmap[idx >> 6] >> (idx & 63)) & 1;
    }

    int getone() const {
        int cnt = 0;
        for(int i = 0; i < cap; i ++)
            for(uint64_t x = bitmap[i]; x; x >>= 1)
                cnt += (x & 1);
        return cnt;
    }

    void print() {
        for(long long i = 0; i < siz; i ++)
            std::cout << get(i) << " ";
        std::cout << "\n";
    }

private:
    uint64_t* bitmap;
    int cap;
    long long siz;
};
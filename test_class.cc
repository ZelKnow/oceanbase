#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
#include "assert.h"

#define ATOMIC_LOAD_ACQ(x) __atomic_load_n((x), __ATOMIC_ACQUIRE)
#define ATOMIC_STORE_REL(x, v) ({ __atomic_store_n((x), (v), __ATOMIC_RELEASE);})
#define LA_ATOMIC_ID 0
#define UNUSED(v) ((void)(v))
#define ATOMIC_BCAS(val, cmpv, newv) __sync_bool_compare_and_swap((val), (cmpv), (newv))


class Version {
public:
	uint64_t data_;
    enum {
        latch_bit = (1ULL << 63),
        inserting_bit = (1ULL << 62),
        splitting_bit = (1ULL << 61),
        vinsert_lowbit = (1ULL << 56),
        vsplit_lowbit = (1ULL << 0),
        split_mask = ~(~0ULL<<56)<<0,
        insert_mask = ~(~0ULL<<5)<<56,
        ulatch_mask = ~(latch_bit | inserting_bit | splitting_bit),
        dirty_mask = inserting_bit | splitting_bit
    };
    Version(): data_(0) {}
    Version get_stable_version() {
        Version stable_version;
        do {
            stable_version = *this;
        } while(stable_version.data_ & dirty_mask);
        //__atomic_thread_fence(__ATOMIC_ACQUIRE);
        return stable_version;
    }
	uint64_t get_vsplit() {
        return (data_ & split_mask) >> 0;
    }
	uint64_t get_vinsert() {
        return (data_ & insert_mask) >> 56;
    }
    bool is_inserting() {
        return data_ & inserting_bit;
    }
    bool is_splitting() {
        return data_ & splitting_bit;
    }
    bool is_latched() {
        return data_ & latch_bit;
    }
    void increase_vinsert() {
        data_ = (data_ & (~insert_mask)) + (((data_ & insert_mask) + vinsert_lowbit) & insert_mask);
    }
    void increase_vsplit() {
        data_ = (data_ & (~split_mask)) + (((data_ & split_mask) + vsplit_lowbit) & split_mask);
    }
    bool has_splited(Version &snapshot_version) {
        //__atomic_thread_fence(__ATOMIC_ACQUIRE);
        if(is_splitting()) {
            return true;
        }
        if(get_vsplit() != snapshot_version.get_vsplit()) {
            return true;
        }
        return false;
    }
    bool has_inserted(Version &snapshot_version) {
        //__atomic_thread_fence(__ATOMIC_ACQUIRE);
        if(is_inserting()) {
            return true;
        }
        if(get_vinsert() != snapshot_version.get_vinsert()) {
            return true;
        }
        return false;
    }
    void set_inserting() {
        data_ = data_ | (uint64_t(1) << 62);
        //__atomic_thread_fence(__ATOMIC_RELEASE);
    }
    void set_splitting() {
        data_ = data_ | (uint64_t(1) << 61);
        //__atomic_thread_fence(__ATOMIC_RELEASE);
    }
    void latch() {
        while(true) {
            uint64_t old_data = data_ & (~latch_bit);
            uint64_t new_data = old_data | latch_bit;
            if (ATOMIC_BCAS(&data_, old_data, new_data)) {
                return;
            }
        }
    }
    void unlatch() {
        if (is_inserting()) {
            increase_vinsert();
        }
        if (is_splitting()) {
            increase_vsplit();
        }
        //__atomic_thread_fence(__ATOMIC_RELEASE);
        data_ = data_ & ulatch_mask;
    }
};

class Path {
};

class Permutation {
private:
	uint64_t data_;
public:
	uint8_t get_count();
	uint8_t at();
	void insert(uint8_t pos, uint8_t value);
};

#define NODE_KEY_COUNT 15
class BtreeNode {
private:
	Version version_;
	Permutation permutation_;
public:
    void set_inserting() {
        // precondition: latched.
        version_.set_inserting();
    }
    void set_splitting() {
        // precondition: latched.
        version_.set_splitting();
    }
    void latch() {
        version_.latch();
    }
    void unlatch() {
        version_.unlatch();
    }
};

template<typename Keytype>
class LeafNode: public BtreeNode {
    Keytype k;
    BtreeNode* father;
public:
    void test() {
        std::cout<<k<<std::endl;
    }
    void set(Keytype kk) {
        k = kk;
    }
};

class MockNode {
public:
    uint64_t a;
    uint64_t b;
    uint64_t c;
    uint64_t id;
};

void test_version() {
    Version v;
    for(int i=1;i<1000;i++) {
        v.increase_vinsert();
        assert(v.get_vinsert() == i % 32);
        v.increase_vsplit();
        assert(v.get_vsplit() == i);
        assert(v.is_inserting() == false);
        assert(v.is_splitting() == false);
        assert(v.is_latched() == false);
    }
    v.set_inserting();
    assert(v.is_inserting());
    v.set_splitting();
    assert(v.is_splitting());
    v.latch();
    assert(v.is_latched());
    int vinsert = v.get_vinsert();
    int vsplit = v.get_vsplit();
    v.unlatch();
    assert(v.is_inserting() == false);
    assert(v.is_splitting() == false);
    assert(v.is_latched() == false);
    assert(v.get_vinsert() == (vinsert+1)%32);
    assert(v.get_vsplit() == vsplit+1);
}

#define WRITE_THREAD_COUNT 4
#define READ_THREAD_COUNT 8
#define REPEAT_COUNT 10000000

void test_sync() {
    Version v;
    MockNode mock_node;
    mock_node.a = mock_node.b = mock_node.c = 0;
    mock_node.id = -1;
    std::thread write_threads[WRITE_THREAD_COUNT];
    std::thread read_threads[READ_THREAD_COUNT];
    for(int i=0;i<WRITE_THREAD_COUNT;i++) {
        write_threads[i] = std::thread([&](int i) {
            for (int i=0;i<REPEAT_COUNT;i++) {
                v.latch();
                assert(mock_node.id == -1);
                v.set_splitting();
                mock_node.id = i;
                mock_node.a++;
                mock_node.b++;
                mock_node.c++;
                mock_node.id = -1;
                v.unlatch();
            }
        }, i);
    }
    for(int i=0;i<READ_THREAD_COUNT;i++) {
        read_threads[i] = std::thread([&](int i) {
            int success_count = 0;
            for (int i=0;i<REPEAT_COUNT;i++) {
                uint64_t a, b, c, id;
                Version snapshot_version;
                do {
                    snapshot_version = v.get_stable_version();
                    a = mock_node.a;
                    b = mock_node.b;
                    c = mock_node.c;
                    id = mock_node.id;
                } while(v.has_splited(snapshot_version));
                assert(a == b && b == c && id == -1);
            }
        }, i);
    }
    for(int i=0;i<WRITE_THREAD_COUNT;i++) {
        write_threads[i].join();
    }
    for(int i=0;i<READ_THREAD_COUNT;i++) {
        read_threads[i].join();
    }
    assert(mock_node.a == mock_node.b && mock_node.b == mock_node.c);
    assert(mock_node.a == WRITE_THREAD_COUNT*REPEAT_COUNT);
    assert(mock_node.id == -1);
}

int main() {
    while(true) {
        test_version();
        test_sync();
    }
}
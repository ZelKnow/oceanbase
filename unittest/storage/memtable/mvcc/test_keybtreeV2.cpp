/**
 * Copyright (c) 2021 OceanBase
 * OceanBase CE is licensed under Mulan PubL v2.
 * You can use this software according to the terms and conditions of the Mulan PubL v2.
 * You may obtain a copy of Mulan PubL v2 at:
 *          http://license.coscl.org.cn/MulanPubL-2.0
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
 * EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
 * MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 * See the Mulan PubL v2 for more details.
 */

#include "storage/memtable/mvcc/ob_keybtree.h"
#include "common/object/ob_object.h"
#include "common/rowkey/ob_store_rowkey.h"
#include "lib/allocator/ob_malloc.h"
#include "lib/random/ob_random.h"
#include "storage/memtable/ob_memtable_key.h"
#include "storage/memtable/mvcc/ob_mvcc_row.h"
#include "storage/memtable/mvcc/ob_query_engine.h"
#include "lib/oblog/ob_log_module.h"
#include <gtest/gtest.h>
#include <thread>
#include <vector>
#include <unistd.h>

namespace oceanbase
{
namespace unittest
{

//#define IS_EQ(x, y) ASSERT_EQ(x, y)
//#define IS_EQ(x, y) EXPECT_EQ(x, y)
#define DUMP_BTREE \
{ \
  FILE *file = fopen("dump_btree.txt", "w+"); \
  btree.dump(file); \
  fclose(file); \
}

#define IS_EQ(x, y) if ((x) != (y)) { abort(); }

#define judge(key, val) \
{ \
  if ((int64_t)(val) >> 3 != get_v(key)) { \
    abort(); \
  } \
}

const char *attr = ObModIds::TEST;

class FakeAllocator : public ObIAllocator
{
public:
  void *alloc(int64_t size) override { return ob_malloc(size, attr); }
  void* alloc(const int64_t size, const ObMemAttr &attr) override
  {
    UNUSED(attr);
    return alloc(size);
  }
  void free(void *ptr) override { ob_free(ptr); }
  static FakeAllocator*get_instance()
  {
    static FakeAllocator allocator;
    return &allocator;
  }
};

typedef keybtree::ObKeyBtree<memtable::ObStoreRowkeyWrapper, memtable::ObMvccRow*> ObKeyBtree;
typedef keybtree::BtreeIterator<memtable::ObStoreRowkeyWrapper, memtable::ObMvccRow*> BtreeIterator;
typedef keybtree::BtreeNodeAllocator<memtable::ObStoreRowkeyWrapper, memtable::ObMvccRow*> BtreeNodeAllocator;
typedef memtable::ObQueryEngine::IteratorAlloc<BtreeIterator> BtreeIteratorAlloc;
typedef memtable::ObQueryEngine::Iterator<BtreeIterator> BIterator;
typedef memtable::ObStoreRowkeyWrapper ObKey;

memtable::ObStoreRowkeyWrapper build_int_key(int64_t key) {
  ObObj *obj = new ObObj(int64_t(key));
  ObStoreRowkey *rk = new ObStoreRowkey();
  rk->assign(obj, 1);
  return memtable::ObStoreRowkeyWrapper(rk);
}

memtable::ObStoreRowkeyWrapper build_char_key(int cnt, int len) {
  auto alloc = FakeAllocator::get_instance();
  void *block;
  block = alloc->alloc(sizeof(ObObj));
  ObObj *obj = new(block) ObObj(ObObjType::ObCharType);
  block = alloc->alloc(sizeof(ObStoreRowkey));
  ObStoreRowkey *rk = new(block) ObStoreRowkey();
  char *aa = (char *)alloc->alloc(len);
  if (cnt < 0) {
    for (int i=0;i<len;i++) {
      if (cnt == -1) {
        aa[i] = '0' + ObRandom::rand(0, 9);
      } else {
        aa[i] = '9';
      }
    }
  } else {
    sprintf(aa, "%0*d", len, cnt);
  }
  obj->set_string(ObObjType::ObCharType, aa, len);
  obj->set_collation_type(CS_TYPE_UTF8MB4_BIN);
  rk->assign(obj, 1);
  return memtable::ObStoreRowkeyWrapper(rk);
}

const int THREAD_COUNT = 32;
const int SINGLE_THREAD_INSERT_COUNT = 300000;
const int KEY_LEN = 16;

void test_concurrent_insert(bool prefetch) {
  std::thread insert_threads[THREAD_COUNT];
  ObKeyBtree *btree;
  BtreeNodeAllocator allocator(*FakeAllocator::get_instance());
  btree = new ObKeyBtree(allocator);
  btree->set_prefetch(prefetch);

  auto alloc = FakeAllocator::get_instance();
  int64_t start_time = rdtsc();
  for (int i=0;i<THREAD_COUNT;i++) {
    insert_threads[i] = std::thread([&](int j) {
      for (int k=0;k<SINGLE_THREAD_INSERT_COUNT;k++) {
        int insert_key = k*THREAD_COUNT+j;
        memtable::ObStoreRowkeyWrapper key = build_char_key(-1, KEY_LEN);
        memtable::ObMvccRow *value = (memtable::ObMvccRow *)(uint64_t(insert_key+1)<<3);
        btree->insert(key, value);
      }
    }, i);
  }
  for (int i=0;i<THREAD_COUNT;i++) {
    insert_threads[i].join();
  }
  int64_t end_time = rdtsc();
  std::cout<<"Btree Insert: "<<end_time-start_time<<std::endl;
  auto start_key = build_char_key(0, KEY_LEN);
  auto end_key = build_char_key(-2, KEY_LEN);
  BtreeIteratorAlloc btree_iter_alloc_;
  BIterator *btree_iter = nullptr;
  btree_iter = btree_iter_alloc_.alloc();
  btree->set_key_range(btree_iter->get_read_handle(), start_key, false, end_key, false, 0);
  memtable::ObStoreRowkeyWrapper prev;
  memtable::ObStoreRowkeyWrapper key;
  memtable::ObMvccRow *value;
  int i = 0;
  int ret = OB_SUCCESS;
  while(OB_SUCC(btree_iter->get_read_handle().get_next(key, value))) {
    i+=1;
    if(i == 1) {
      prev = key;
      continue;
    }
    ASSERT_EQ(prev.get_ptr()->compare(*key.get_ptr()), -1);
    prev = key;
  }
  std::cout<<i<<std::endl;
}

TEST(TestBtreePrefetchInsert, insert_test) {
  test_concurrent_insert(true);
}

TEST(TestBtreeInsert, insert_test) {
  test_concurrent_insert(false);
}

const int TOTAL_INSERT_COUNT = 2000000;
const int SINGLE_THREAD_SCAN_COUNT = 3000;

void test_scan(bool prefetch) {
  std::thread scan_threads[THREAD_COUNT];
  ObKeyBtree *btree;
  BtreeNodeAllocator allocator(*FakeAllocator::get_instance());
  btree = new ObKeyBtree(allocator);
  btree->set_prefetch(prefetch);

  auto alloc = FakeAllocator::get_instance();
  for (int k=0;k<TOTAL_INSERT_COUNT;k++) {
    memtable::ObStoreRowkeyWrapper key = build_char_key(k, KEY_LEN);
    memtable::ObMvccRow *value = nullptr;
    btree->insert(key, value);
  }
  std::cout<<"inserted"<<std::endl;
  int64_t start_time = rdtsc();
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i] = std::thread([&](int j) {
      for (int k=0;k<SINGLE_THREAD_SCAN_COUNT;k++) {
        int insert_key = ObRandom::rand(0, TOTAL_INSERT_COUNT-1);
        //std::cout<<insert_key<<std::endl;
        auto start_key = build_char_key(insert_key, KEY_LEN);
        auto end_key = build_char_key(insert_key+10000, KEY_LEN);
        BtreeIteratorAlloc btree_iter_alloc_;
        BIterator *btree_iter = nullptr;
        btree_iter = btree_iter_alloc_.alloc();
        btree->set_key_range(btree_iter->get_read_handle(), start_key, false, end_key, false, 0);
        memtable::ObStoreRowkeyWrapper key;
        memtable::ObMvccRow *value;
        int ret = OB_SUCCESS;
        while(OB_SUCC(btree_iter->get_read_handle().get_next(key, value))) {}
        //ASSERT_EQ(btree_iter->get_read_handle().get_range(), 0);
      }
    }, i);
  }
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i].join();
  }
  int64_t end_time = rdtsc();
  std::cout<<"Btree Scan: "<<end_time-start_time<<std::endl;

}

TEST(TestBtreePrefetchScan, scan_test) {
  test_scan(true);
}

TEST(TestBtreeScan, scan_test) {
  test_scan(false);
}

const int PRE_INSERT_COUNT = 1500000;

void test_scan_insert(bool prefetch) {
  std::thread scan_threads[THREAD_COUNT];
  ObKeyBtree *btree;
  BtreeNodeAllocator allocator(*FakeAllocator::get_instance());
  btree = new ObKeyBtree(allocator);
  btree->set_prefetch(prefetch);
  auto alloc = FakeAllocator::get_instance();
  for (int k=0;k<TOTAL_INSERT_COUNT;k++) {
    memtable::ObStoreRowkeyWrapper key = build_char_key(-1, KEY_LEN);
    memtable::ObMvccRow *value = nullptr;
    btree->insert(key, value);
  }
  std::cout<<"inserted"<<std::endl;
  int64_t start_time = rdtsc();
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i] = std::thread([&](int j) {
      for (int k=0;k<SINGLE_THREAD_SCAN_COUNT;k++) {
        int insert_key;
        for (int a=0;a<10;a++) {
          insert_key = 10*(k*THREAD_COUNT+j) + a + PRE_INSERT_COUNT;
          memtable::ObStoreRowkeyWrapper key = build_char_key(-1, KEY_LEN);
          memtable::ObMvccRow *value = nullptr;
          btree->insert(key, value);
        }
        insert_key = ObRandom::rand(0, insert_key);
        auto start_key = build_char_key(-1, KEY_LEN);
        auto end_key = build_char_key(-1, KEY_LEN);
        if (end_key.get_ptr()->compare(*start_key.get_ptr()) == -1) {
          auto tmp_key = end_key;
          end_key = start_key;
          start_key = tmp_key;
        }
        //auto start_key = build_char_key(-1, KEY_LEN);
        //auto end_key = build_char_key(-2, KEY_LEN);
        BtreeIteratorAlloc btree_iter_alloc_;
        BIterator *btree_iter = nullptr;
        btree_iter = btree_iter_alloc_.alloc();
        btree->set_key_range(btree_iter->get_read_handle(), start_key, false, end_key, false, 0);
        memtable::ObStoreRowkeyWrapper key;
        memtable::ObMvccRow *value;
        int ret = OB_SUCCESS;
        while(OB_SUCC(btree_iter->get_read_handle().get_next(key, value))) {}
      }
    }, i);
  }
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i].join();
  }
  int64_t end_time = rdtsc();
  std::cout<<"Btree Scan + Insert: "<<end_time-start_time<<std::endl;
  //sleep(10);
}
TEST(TestBtreePrefetchHybrid, Hybrid_test) {
  test_scan_insert(true);
}

TEST(TestBtreeHybrid, Hybrid_test) {
  test_scan_insert(false);
}

}
}

int main(int argc, char **argv)
{
  //oceanbase::unittest::BIND_CPU(pthread_self());
  oceanbase::common::ObLogger::get_logger().set_file_name("test_keybtreeV2.log", true);
  oceanbase::common::ObLogger::get_logger().set_log_level("INFO");
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

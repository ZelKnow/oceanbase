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

#include "storage/memtable/mvcc/ob_keybtreeV2.h"
#include "storage/memtable/mvcc/ob_keybtree.h"
#include "storage/memtable/mvcc/ob_query_engine.h"
#include "lib/allocator/ob_malloc.h"
#include "lib/oblog/ob_log.h"
#include "lib/random/ob_random.h"
#include "common/object/ob_object.h"
#include <gtest/gtest.h>
#include <thread>
#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_set>

namespace oceanbase {
namespace unittest {
using namespace oceanbase::keybtreeV2;

#define DUMP_BTREE                              \
  {                                             \
    FILE *file = fopen("dump_btree.txt", "w+"); \
    btree.dump(file);                           \
    fclose(file);                               \
  }

class FakeAllocator : public ObIAllocator {
public:
  const char *attr = ObModIds::TEST;
  void *alloc(int64_t size) override
  {
    return ob_malloc(size, attr);
  }
  void *alloc(const int64_t size, const ObMemAttr &attr) override
  {
    UNUSED(attr);
    return alloc(size);
  }
  void free(void *ptr) override
  {
    ob_free(ptr);
  }
  static FakeAllocator *get_instance()
  {
    static FakeAllocator allocator;
    return &allocator;
  }
};

class FakeKey {
public:
  FakeKey() : obj_(nullptr)
  {}
  FakeKey(ObObj *obj) : obj_(obj)
  {}
  void set_int(int64_t data)
  {
    obj_->set_int(data);
  }
  int compare(FakeKey other, int &cmp) const
  {
    return obj_->compare(*other.obj_, cmp);
  }
  int64_t to_string(char *buf, const int64_t limit) const
  {
    return obj_->to_string(buf, limit);
  }

  ObObj *get_ptr() const
  {
    return obj_;
  }
  ObObj *obj_;
};

FakeKey build_int_key(int64_t key)
{
  auto alloc = FakeAllocator::get_instance();
  void *block = alloc->alloc(sizeof(ObObj));
  ObObj *obj = new (block) ObObj(key);
  return FakeKey(obj);
}

void free_key(FakeKey &key)
{
  auto alloc = FakeAllocator::get_instance();
  alloc->free((void *)key.obj_);
}

/*ObStoreRowkeyWrapper build_char_key(int cnt, int len) {
  auto alloc = FakeAllocator::get_instance();
  void *block = alloc->alloc(sizeof(ObObj));
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
  return ObStoreRowkeyWrapper(rk);
}*/

using NewBtree = ObKeyBtree<FakeKey, int64_t*>;
using NewBtreeNodeAllocator = BtreeNodeAllocator<FakeKey, int64_t*>;
using NewBtreeIterator = BtreeIterator<FakeKey, int64_t*>;
using OldBtreeNodeAllocator = keybtree::BtreeNodeAllocator<FakeKey, int64_t*>;
using OldBtree = keybtree::ObKeyBtree<FakeKey, int64_t*>;
using OldBtreeIterator = keybtree::BtreeIterator<FakeKey, int64_t*>;

void test_concurrent_insert(bool is_new) {
  constexpr uint64_t KEY_NUM = 3000000;
  constexpr uint64_t THREAD_COUNT = 10;
  constexpr uint64_t PER_THREAD_INSERT_COUNT = KEY_NUM / THREAD_COUNT;

  NewBtree *new_btree;
  OldBtree *old_btree;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  NewBtreeNodeAllocator new_node_allocator(*allocator);
  OldBtreeNodeAllocator old_node_allocator(*allocator);

  if(is_new) {
    new_btree = new NewBtree(new_node_allocator);
    ASSERT_EQ(new_btree->init(), OB_SUCCESS);
  } else {
    old_btree = new OldBtree(old_node_allocator);
    ASSERT_EQ(old_btree->init(), OB_SUCCESS);
  }

  std::thread threads[THREAD_COUNT];
  std::vector<std::vector<FakeKey>> data(THREAD_COUNT, std::vector<FakeKey>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = build_int_key(THREAD_COUNT * j + i);
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  // concurrent insert
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id] = std::thread(
      [&](int i) {
        for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
          int64_t *val = reinterpret_cast<int64_t*>(data[i][j].get_ptr());
          if(is_new) {
            new_btree->insert(data[i][j], val);
          } else {
            old_btree->insert(data[i][j], val);
          }
        }
      },
      thread_id);
  }

  uint64_t start_time = rdtsc();
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id].join();
  }
  if(is_new) {
    std::cout<<"New Btree Concurrent Insert: "<<rdtsc()-start_time<<std::endl;
  } else {
    std::cout<<"Old Btree Concurrent Insert: "<<rdtsc()-start_time<<std::endl;
  }

  for(int i=0;i<data.size();i++) {
    for(int j=0;j<data[i].size();j++) {
      allocator->free(data[i][j].get_ptr());
    }
  }
}

TEST(TestNewBtreeConcurrentInsert, smoke_test)
{
  test_concurrent_insert(true);
}

TEST(TestOldBtreeConcurrentInsert, smoke_test)
{
  test_concurrent_insert(false);
}

void test_concurrent_insert_scan(bool is_new) {
  constexpr int KEY_NUM = 3000000;
  constexpr int THREAD_COUNT = 10;
  constexpr int PER_THREAD_INSERT_COUNT = KEY_NUM / THREAD_COUNT;
  constexpr int INSERT_SCAN_RATIO = 20;

  NewBtree *new_btree;
  OldBtree *old_btree;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  NewBtreeNodeAllocator new_node_allocator(*allocator);
  OldBtreeNodeAllocator old_node_allocator(*allocator);

  if(is_new) {
    new_btree = new NewBtree(new_node_allocator);
    ASSERT_EQ(new_btree->init(), OB_SUCCESS);
  } else {
    old_btree = new OldBtree(old_node_allocator);
    ASSERT_EQ(old_btree->init(), OB_SUCCESS);
  }

  // constructing insert keys
  std::vector<std::vector<FakeKey>> data(THREAD_COUNT, std::vector<FakeKey>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = build_int_key(THREAD_COUNT * j + i);
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  std::thread threads[THREAD_COUNT];
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id] = std::thread(
        [&](int i) {
          FakeKey start_key = build_int_key(0);
          FakeKey end_key = build_int_key(0);
          bool is_backward;
          for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
            int64_t *val = reinterpret_cast<int64_t*>(data[i][j].get_ptr());
            FakeKey key;
            if(is_new) {
              new_btree->insert(data[i][j], val);
            } else {
              old_btree->insert(data[i][j], val);
            }
            if (j % INSERT_SCAN_RATIO != 0) {
              continue;
            } else if (j % (2*INSERT_SCAN_RATIO) == 0) {
              int64_t start_int = ObRandom::rand(0, KEY_NUM);
              int64_t end_int = start_int + 2000;
              start_key.set_int(start_int);
              end_key.set_int(end_int);
              is_backward = false;
            } else {
              int64_t start_int = ObRandom::rand(0, KEY_NUM);
              int64_t end_int = start_int + 2000;
              start_key.set_int(end_int);
              end_key.set_int(start_int);
              is_backward = true;
            }
            if(is_new) {
              NewBtreeIterator iter;
              iter.init(*new_btree);
              iter.set_key_range(start_key, true, end_key, true);
              int ret = OB_SUCCESS;
              while (OB_SUCC(iter.get_next(key, val))) {}
              ASSERT_EQ(ret, OB_ITER_END);
            } else {
              OldBtreeIterator iter;
              iter.init(*old_btree);
              iter.set_key_range(start_key, true, end_key, true);
              int ret = OB_SUCCESS;
              while(OB_SUCC(iter.get_next(key, val))) {}
              ASSERT_EQ(ret, OB_ITER_END);
            }
          }
        },
        thread_id);
  }
  uint64_t start_time = rdtsc();
  for(int i=0;i<THREAD_COUNT;i++) {
    threads[i].join();
  }
  if(is_new) {
    std::cout<<"New Btree Concurrent Insert And Scan: "<<rdtsc()-start_time<<std::endl;
  } else {
    std::cout<<"Old Btree Concurrent Insert And Scan: "<<rdtsc()-start_time<<std::endl;
  }

  for(int i=0;i<data.size();i++) {
    for(int j=0;j<data[i].size();j++) {
      allocator->free(data[i][j].get_ptr());
    }
  }
}

TEST(TestNewConcurrentInsertAndScan, smoke_test) {
  test_concurrent_insert_scan(true);
}

TEST(TestOldConcurrentInsertAndScan, smoke_test) {
  test_concurrent_insert_scan(false);
}

}  // namespace unittest
}  // namespace oceanbase

int main(int argc, char **argv)
{
  // oceanbase::unittest::BIND_CPU(pthread_self());
  oceanbase::common::ObLogger::get_logger().set_file_name("test_keybtreeV2.log", true);
  oceanbase::common::ObLogger::get_logger().set_log_level("INFO");
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
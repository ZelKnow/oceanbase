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
#include <atomic>
#include <unordered_set>

#define KEYTYPE_INT

namespace oceanbase {
namespace unittest {
using namespace oceanbase::keybtree;

#define DUMP_BTREE                              \
  {                                             \
    FILE *file = fopen("dump_btree.txt", "w+"); \
    btree.dump(file);                           \
    fclose(file);                               \
  }

class FakeAllocator : public ObIAllocator {
public:
  FakeAllocator(): ObIAllocator() {}
  const char *attr = ObModIds::TEST;
  void *alloc(int64_t size) override
  {
    return ob_malloc(size, attr);
  }
  void *alloc_key(int64_t size)
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

int my_compare(const char *s, const int slen, const char *t, const int tlen, int &cmp, int &prefix)
{
  int start = prefix;
  int min_len = min(slen, tlen);
  for(int i=start;i<min_len;i++) {
    if (*(s+i) != *(t+i)) {
      cmp = int(*(s+i)) - int(*(t+i));
      prefix = start + i;
      return OB_SUCCESS;
    }
  }
  if (slen == tlen) {
    cmp = 0;
    prefix = slen;
  } else if (slen > tlen) {
    cmp = 1;
    prefix = min_len;
  } else {
    cmp = -1;
    prefix = min_len;
  }
  return OB_SUCCESS;
}

class FakeKey {
public:
  FakeKey() : obj_(nullptr)
  {}
  FakeKey(ObObj *obj) : obj_(obj)
  {}
  void set(int64_t data, int len) {
    #ifdef KEYTYPE_INT
    set_int(data);
    #else
    set_char(data, len);
    #endif
  }
  void set_int(int64_t data)
  {
    obj_->set_int(data);
  }
  void set_char(int64_t data, int len)
  {
    char *aa = const_cast<char*>(obj_->get_string_ptr());
    sprintf(aa, "%0*ld", len, data);
  }
  int compare(FakeKey other, int &cmp) const
  {
    return obj_->compare(*other.obj_, cmp);
  }
  int compare(FakeKey other, int &cmp, int &prefix) const 
  {
    #ifdef KEYTYPE_INT
    prefix = 0;
    return compare(other, cmp);
    #endif
    int ret = OB_SUCCESS;
    prefix = 0;
    if(obj_->get_collation_type() != other.obj_->get_collation_type()) {
      ret = OB_ERR_UNEXPECTED;
    } else if(obj_->get_collation_type() != CS_TYPE_UTF8MB4_BIN) {
      ret = compare(other, cmp);
    } else {
      ret = my_compare(obj_->v_.string_, obj_->val_len_, other.obj_->v_.string_, other.obj_->val_len_, cmp, prefix);
    }
    return ret;
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
  void *block = alloc->alloc_key(sizeof(ObObj));
  EXPECT_TRUE(OB_NOT_NULL(block));
  ObObj *obj = new (block) ObObj(key);
  return FakeKey(obj);
}

void free_key(FakeKey &key)
{
  auto alloc = FakeAllocator::get_instance();
  if (key.obj_->get_type() == ObObjType::ObCharType) {
    alloc->free((void *)key.obj_->get_string_ptr());
  }
  alloc->free((void *)key.obj_);
}

FakeKey build_char_key(int64_t key, int len) {
  auto alloc = FakeAllocator::get_instance();
  void *block = alloc->alloc_key(sizeof(ObObj));
  EXPECT_TRUE(OB_NOT_NULL(block));
  ObObj *obj = new(block) ObObj(ObObjType::ObCharType);
  char *aa = (char *)alloc->alloc_key(len+1);
  sprintf(aa, "%0*ld", len, key);
  obj->set_string(ObObjType::ObCharType, aa, len);
  obj->set_collation_type(CS_TYPE_UTF8MB4_BIN);
  return FakeKey(obj);
}

FakeKey build_key(int64_t key, int len) {
  #ifdef KEYTYPE_INT
  return build_int_key(key);
  #else
  return build_char_key(key, len);
  #endif
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

using Btree = ObKeyBtree<FakeKey, int64_t*>;
using BtreeNodeAllocator = BtreeNodeAllocator<FakeKey, int64_t*>;
using BtreeIterator = BtreeIterator<FakeKey, int64_t*>;
uint64_t THREAD_COUNT;
uint64_t KEY_NUM;
uint64_t INSERT_SCAN_RATIO = 10;

void test_concurrent_insert(bool is_new) {
  constexpr int KEY_LEN = 10;
  uint64_t PER_THREAD_INSERT_COUNT = KEY_NUM / THREAD_COUNT;
  Btree *new_btree;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator new_node_allocator(*allocator);

  new_btree = new Btree(new_node_allocator);
  ASSERT_EQ(new_btree->init(), OB_SUCCESS);

  std::vector<std::thread> threads(THREAD_COUNT);
  std::vector<std::vector<FakeKey>> data(THREAD_COUNT, std::vector<FakeKey>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = build_key(THREAD_COUNT * j + i, KEY_LEN);
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  volatile uint64_t start_time = rdtsc();

  // concurrent insert
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id] = std::thread(
      [&](int i) {
        for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
          int64_t *val = reinterpret_cast<int64_t*>(data[i][j].get_ptr());
        new_btree->insert(data[i][j], val);
        }
      },
      thread_id);
  }

  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id].join();
  }

  volatile uint64_t end_time = rdtsc();

    std::cout<<"New Btree Concurrent Insert: "<<end_time-start_time<<std::endl;
    new_btree->destroy();

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

void test_concurrent_insert_scan(bool is_new) {
  uint64_t PER_THREAD_INSERT_COUNT = KEY_NUM / THREAD_COUNT;
  constexpr int KEY_LEN = 10;
  Btree *new_btree;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator new_node_allocator(*allocator);

    new_btree = new Btree(new_node_allocator);
    ASSERT_EQ(new_btree->init(), OB_SUCCESS);

  // constructing insert keys
  std::vector<std::vector<FakeKey>> data(THREAD_COUNT, std::vector<FakeKey>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = build_key(THREAD_COUNT * j + i, KEY_LEN);
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  volatile uint64_t start_time = rdtsc();

  std::vector<std::thread> threads(THREAD_COUNT);
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id] = std::thread(
        [&](int i) {
          FakeKey start_key = build_key(0, KEY_LEN);
          FakeKey end_key = build_key(0, KEY_LEN);
          bool is_backward;
          for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
            int64_t *val = reinterpret_cast<int64_t*>(data[i][j].get_ptr());
            FakeKey key;
              new_btree->insert(data[i][j], val);
            if (j % INSERT_SCAN_RATIO != 0) {
              continue;
            } else if (j % (2*INSERT_SCAN_RATIO) == 0) {
              int64_t start_int = ObRandom::rand(0, KEY_NUM);
              int64_t end_int = start_int + 4000;
              start_key.set(start_int, KEY_LEN);
              end_key.set(end_int, KEY_LEN);
              is_backward = false;
            } else {
              int64_t start_int = ObRandom::rand(0, KEY_NUM);
              int64_t end_int = start_int + 4000;
              start_key.set(end_int, KEY_LEN);
              end_key.set(start_int, KEY_LEN);
              is_backward = true;
            }
              BtreeIterator iter;
              iter.init(*new_btree);
              iter.set_key_range(start_key, true, end_key, true, 0);
              int ret = OB_SUCCESS;
              while (OB_SUCC(iter.get_next(key, val))) {}
              ASSERT_EQ(ret, OB_ITER_END);
          }
        },
        thread_id);
  }
  for(int i=0;i<THREAD_COUNT;i++) {
    threads[i].join();
  }

  volatile uint64_t end_time = rdtsc();

    std::cout<<"New Btree Concurrent Insert And Scan: "<<end_time-start_time<<std::endl;
    new_btree->destroy();

  for(int i=0;i<data.size();i++) {
    for(int j=0;j<data[i].size();j++) {
      allocator->free(data[i][j].get_ptr());
    }
  }
}

TEST(TestNewConcurrentInsertAndScan, smoke_test) {
  test_concurrent_insert_scan(true);
}

}  // namespace unittest
}  // namespace oceanbase

int main(int argc, char **argv)
{
  // oceanbase::unittest::BIND_CPU(pthread_self());

  oceanbase::common::ObLogger::get_logger().set_file_name("test_keybtreeV2.log", true);
  oceanbase::common::ObLogger::get_logger().set_log_level("INFO");
  ::testing::InitGoogleTest(&argc, argv);
  oceanbase::unittest::THREAD_COUNT = atoi(argv[1]);
  oceanbase::unittest::KEY_NUM = atoi(argv[2]);
  return RUN_ALL_TESTS();
}
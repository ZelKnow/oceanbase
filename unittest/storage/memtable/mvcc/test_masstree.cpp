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
#include "storage/memtable/mvcc/masstree-beta/masstree_api.h"
#include "storage/memtable/mvcc/ARTSynchronized/OptimisticLockCoupling/Tree.h"
#include "storage/memtable/mvcc/ARTSynchronized/ROWEX/Tree.h"
//#include "storage/memtable/mvcc/ARTOLC/Tree.h"
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

typedef memtable::ObStoreRowkeyWrapper MasstreeKey;
typedef ART_ROWEX::Tree ART_TREE;
const char *attr = ObModIds::TEST;

void init_key(MasstreeKey *ptr, int64_t key)
{
  ptr->get_rowkey()->get_rowkey().get_obj_ptr()[0].set_int(key);
}

int alloc_key(MasstreeKey *&ret_key, int64_t key)
{
  int ret = OB_SUCCESS;
  ObObj *obj_ptr = nullptr;
  ObStoreRowkey *storerowkey = nullptr;
  if (OB_ISNULL(obj_ptr = (ObObj *)ob_malloc(sizeof(ObObj), attr)) || OB_ISNULL(new(obj_ptr)ObObj(key))) {
    ret = OB_ALLOCATE_MEMORY_FAILED;
  } else if (OB_ISNULL(storerowkey = (ObStoreRowkey *)ob_malloc(sizeof(ObStoreRowkey), attr))) {
    ret = OB_ALLOCATE_MEMORY_FAILED;
  } else if (OB_ISNULL(ret_key = (MasstreeKey *)ob_malloc(sizeof(MasstreeKey), attr)) || OB_ISNULL(new(ret_key)MasstreeKey(storerowkey))) {
    ret = OB_ALLOCATE_MEMORY_FAILED;
  }
  return ret;
}

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
typedef memtable::ObQueryEngine::IteratorAlloc<masstree::MasstreeIterator> MasstreeIteratorAlloc;
typedef memtable::ObQueryEngine::Iterator<masstree::MasstreeIterator> MIterator;
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
const int SINGLE_THREAD_INSERT_COUNT = 200000;
const int KEY_LEN = 48;

enum TREETYPE {
  BTREE = 0,
  MASSTREE,
  ART,
};

void loadKey(TID tid, Key &key) {
    // Store the key of the tuple into the key vector
    // Implementation is database specific
    //char aa[KEY_LEN];
    sprintf((char *)key.stackKey, "%0*d", KEY_LEN, int(tid));
    key.len = KEY_LEN;
    key.data = key.stackKey;
}

void test_concurrent_insert(TREETYPE tree_type) {
  std::thread insert_threads[THREAD_COUNT];
  ObKeyBtree *btree;
  masstree::MassTreeIndex *masstree;
  ART_TREE *art;
  if (tree_type == MASSTREE) {
    masstree = new masstree::MassTreeIndex;
  } else if (tree_type == BTREE) {
    BtreeNodeAllocator allocator(*FakeAllocator::get_instance());
    btree = new ObKeyBtree(allocator);
  } else if (tree_type == ART) {
    art = new ART_TREE(loadKey);
  }
  auto alloc = FakeAllocator::get_instance();
  int64_t start_time = rdtsc();
  for (int i=0;i<THREAD_COUNT;i++) {
    insert_threads[i] = std::thread([&](int j) {
      threadinfo *ti;
      ART::ThreadInfo *t;
      if (tree_type == MASSTREE) {
        ti = threadinfo::make(threadinfo::TI_PROCESS, -1);
      } else if (tree_type == ART) {
        t = new ThreadInfo(art->getThreadInfo());
      }
      for (int k=0;k<SINGLE_THREAD_INSERT_COUNT;k++) {
        int insert_key = k*THREAD_COUNT+j;
        memtable::ObStoreRowkeyWrapper key = build_char_key(insert_key, KEY_LEN);
        memtable::ObMvccRow *value;
        if (tree_type == MASSTREE) {
          masstree->insert(key, value, ti);
        } else if (tree_type == BTREE) {
          btree->insert(key, value);
        } else if (tree_type == ART) {
          Key art_key;
          art_key.set(key.get_ptr()->get_string_ptr(), key.get_ptr()->get_string_len());
          art->insert(art_key, insert_key, *t);
        }
      }
    }, i);
  }
  for (int i=0;i<THREAD_COUNT;i++) {
    insert_threads[i].join();
  }
  int64_t end_time = rdtsc();
  /*
  if (tree_type == ART) {
    auto start_key = build_char_key(0, KEY_LEN);
    auto end_key = build_char_key(1000, KEY_LEN);
    Key art_start_key, art_end_key;
    art_start_key.set(start_key.get_ptr()->get_string_ptr(), start_key.get_ptr()->get_string_len());
    art_end_key.set(end_key.get_ptr()->get_string_ptr(), end_key.get_ptr()->get_string_len());
    TID results[1000];
    size_t result_count = 0;
    auto t = art->getThreadInfo();
    art->lookupRange(art_start_key, art_end_key, art_start_key, results, 1000, result_count, t);
    for (int i=0;i<1000;i++) {
      ASSERT_EQ(results[i], i);
    }
  }*/
  if (tree_type == MASSTREE) {
    std::cout<<"Masstree Insert: "<<end_time-start_time<<std::endl;
  } else if (tree_type == BTREE) {
    std::cout<<"Btree Insert: "<<end_time-start_time<<std::endl;
  } else if (tree_type == ART) {
    std::cout<<"ART Insert: "<<end_time-start_time<<std::endl;
  }
}

TEST(TestMasstreeInsert, insert_test) {
  test_concurrent_insert(MASSTREE);
}

TEST(TestBtreeInsert, insert_test) {
  test_concurrent_insert(BTREE);
}

TEST(TestARTInsert, insert_test) {
  test_concurrent_insert(ART);
}

const int TOTAL_INSERT_COUNT = 2000000;
const int SINGLE_THREAD_SCAN_COUNT = 10000;


void test_scan(TREETYPE tree_type) {
  std::thread scan_threads[THREAD_COUNT];
  ObKeyBtree *btree;
  masstree::MassTreeIndex *masstree;
  ART_TREE *art;
  if (tree_type == MASSTREE) {
    masstree = new masstree::MassTreeIndex;
  } else if (tree_type == BTREE) {
    BtreeNodeAllocator allocator(*FakeAllocator::get_instance());
    btree = new ObKeyBtree(allocator);
  } else if (tree_type == ART) {
    art = new ART_TREE(loadKey);
  }
  auto alloc = FakeAllocator::get_instance();
  threadinfo *ti;
  ART::ThreadInfo *t;
  if (tree_type == MASSTREE) {
    ti = threadinfo::make(threadinfo::TI_MAIN, -1);
  } else if (tree_type == ART) {
    t = new ThreadInfo(art->getThreadInfo());
  }
  for (int k=0;k<TOTAL_INSERT_COUNT;k++) {
    memtable::ObStoreRowkeyWrapper key = build_char_key(k, KEY_LEN);
    memtable::ObMvccRow *value = nullptr;
    if (tree_type == MASSTREE) {
      masstree->insert(key, value, ti);
    } else if (tree_type == BTREE) {
      btree->insert(key, value);
    } else if (tree_type == ART) {
      Key art_key;
      art_key.set(key.get_ptr()->get_string_ptr(), key.get_ptr()->get_string_len());
      art->insert(art_key, k, *t);    
    }
  }
  std::cout<<"inserted"<<std::endl;
  int64_t start_time = rdtsc();
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i] = std::thread([&](int j) {
      threadinfo *ti;
      ART::ThreadInfo *t;
      
      if (tree_type == MASSTREE) {
        ti = threadinfo::make(threadinfo::TI_PROCESS, -1);
      } else if (tree_type == ART) {
        t = new ThreadInfo(art->getThreadInfo());
      }
      for (int k=0;k<SINGLE_THREAD_SCAN_COUNT;k++) {
        int insert_key = ObRandom::rand(0, TOTAL_INSERT_COUNT-1);
        //std::cout<<insert_key<<std::endl;
        auto start_key = build_char_key(insert_key, KEY_LEN);
        auto end_key = build_char_key(TOTAL_INSERT_COUNT, KEY_LEN);
        if (tree_type == MASSTREE) {
          MasstreeIteratorAlloc masstree_iter_alloc_;
          MIterator *masstree_iter = nullptr;
          masstree_iter = masstree_iter_alloc_.alloc();
          masstree_iter->get_read_handle().set_key_range(start_key, false, end_key, false, 1000, *masstree, ti);
          //std::cout<<masstree_iter->get_read_handle().kv_queue_.size()<<" "<<masstree_iter->get_read_handle().get_version()<<std::endl;
          //ASSERT_EQ(masstree_iter->get_read_handle().get_version(), 0);
        } else if (tree_type == BTREE) {
          BtreeIteratorAlloc btree_iter_alloc_;
          BIterator *btree_iter = nullptr;
          btree_iter = btree_iter_alloc_.alloc();
          btree->set_key_range(btree_iter->get_read_handle(), start_key, false, end_key, false, 0);
          memtable::ObStoreRowkeyWrapper key;
          memtable::ObMvccRow *value;
          int ret = OB_SUCCESS;
          while(OB_SUCC(btree_iter->get_read_handle().get_next(key, value))) {}
          //ASSERT_EQ(btree_iter->get_read_handle().get_range(), 0);
        } else if (tree_type == ART) {
          Key art_start_key, art_end_key;
          art_start_key.set(start_key.get_ptr()->get_string_ptr(), start_key.get_ptr()->get_string_len());
          art_end_key.set(end_key.get_ptr()->get_string_ptr(), end_key.get_ptr()->get_string_len());
          TID results[1000];
          size_t result_count = 0;
          art->lookupRange(art_start_key, art_end_key, art_end_key, results, 1000, result_count, *t);
        }
      }
    }, i);
  }
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i].join();
  }
  int64_t end_time = rdtsc();
  if (tree_type == MASSTREE) {
    std::cout<<"Masstree Scan: "<<end_time-start_time<<std::endl;
  } else if (tree_type == BTREE) {
    std::cout<<"Btree Scan: "<<end_time-start_time<<std::endl;
  } else if (tree_type == ART) {
    std::cout<<"ART Scan: "<<end_time-start_time<<std::endl;
  }
}

TEST(TestMasstreeScan, scan_test) {
  test_scan(MASSTREE);
}

TEST(TestBtreeScan, scan_test) {
  test_scan(BTREE);
}

TEST(TestARTScan, scan_test) {
  test_scan(ART);
}

const int PRE_INSERT_COUNT = 1000000;

void test_scan_insert(TREETYPE tree_type) {
  std::thread scan_threads[THREAD_COUNT];
  ObKeyBtree *btree;
  masstree::MassTreeIndex *masstree;
  ART_TREE *art;
  if (tree_type == MASSTREE) {
    masstree = new masstree::MassTreeIndex;
  } else if (tree_type == BTREE) {
    BtreeNodeAllocator allocator(*FakeAllocator::get_instance());
    btree = new ObKeyBtree(allocator);
  } else if (tree_type == ART) {
    art = new ART_TREE(loadKey);
  }
  auto alloc = FakeAllocator::get_instance();
  threadinfo *ti;
  ART::ThreadInfo *t;
  if (tree_type == MASSTREE) {
    ti = threadinfo::make(threadinfo::TI_MAIN, -1);
  } else if (tree_type == ART) {
    t = new ThreadInfo(art->getThreadInfo());
  }
  for (int k=0;k<TOTAL_INSERT_COUNT;k++) {
    memtable::ObStoreRowkeyWrapper key = build_char_key(-1, KEY_LEN);
    memtable::ObMvccRow *value = nullptr;
    if (tree_type == MASSTREE) {
      masstree->insert(key, value, ti);
    } else if (tree_type == BTREE) {
      btree->insert(key, value);
    } else if (tree_type == ART) {
      Key art_key;
      art_key.set(key.get_ptr()->get_string_ptr(), key.get_ptr()->get_string_len());
      art->insert(art_key, k, *t);    
    }
  }
  std::cout<<"inserted"<<std::endl;
  int64_t start_time = rdtsc();
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i] = std::thread([&](int j) {
      threadinfo *ti;
      ART::ThreadInfo *t;
      if (tree_type == MASSTREE) {
        ti = threadinfo::make(threadinfo::TI_PROCESS, -1);
      } else if (tree_type == ART) {
        t = new ThreadInfo(art->getThreadInfo());
      }
      for (int k=0;k<SINGLE_THREAD_SCAN_COUNT;k++) {
        int insert_key;
        for (int a=0;a<5;a++) {
          insert_key = 5*(k*THREAD_COUNT+j) + a + PRE_INSERT_COUNT;
          memtable::ObStoreRowkeyWrapper key = build_char_key(insert_key, KEY_LEN);
          memtable::ObMvccRow *value = nullptr;
          if (tree_type == MASSTREE) {
            masstree->insert(key, value, ti);
          } else if (tree_type == BTREE) {
            btree->insert(key, value);
          } else if (tree_type == ART) {
            Key art_key;
            art_key.set(key.get_ptr()->get_string_ptr(), key.get_ptr()->get_string_len());
            art->insert(art_key, insert_key, *t);
          }
        }
        auto start_key = build_char_key(ObRandom::rand(0, insert_key), KEY_LEN);
        auto end_key = build_char_key(TOTAL_INSERT_COUNT + insert_key, KEY_LEN);
        //auto start_key = build_char_key(-1, KEY_LEN);
        //auto end_key = build_char_key(-2, KEY_LEN);
        if (tree_type == MASSTREE) {
          MasstreeIteratorAlloc masstree_iter_alloc_;
          MIterator *masstree_iter = nullptr;
          masstree_iter = masstree_iter_alloc_.alloc();
          masstree_iter->get_read_handle().set_key_range(start_key, false, end_key, false, 1000, *masstree, ti);
        } else if (tree_type == BTREE) {
          BtreeIteratorAlloc btree_iter_alloc_;
          BIterator *btree_iter = nullptr;
          btree_iter = btree_iter_alloc_.alloc();
          btree->set_key_range(btree_iter->get_read_handle(), start_key, false, end_key, false, 0);
          memtable::ObStoreRowkeyWrapper key;
          memtable::ObMvccRow *value;
          int ret = OB_SUCCESS;
          while(OB_SUCC(btree_iter->get_read_handle().get_next(key, value))) {}
        } else if (tree_type == ART) {
          Key art_start_key, art_end_key;
          art_start_key.set(start_key.get_ptr()->get_string_ptr(), start_key.get_ptr()->get_string_len());
          art_end_key.set(end_key.get_ptr()->get_string_ptr(), end_key.get_ptr()->get_string_len());
          TID results[1000];
          size_t result_count = 0;
          art->lookupRange(art_start_key, art_end_key, art_end_key, results, 1000, result_count, *t);
        }
      }
    }, i);
  }
  for (int i=0;i<THREAD_COUNT;i++) {
    scan_threads[i].join();
  }
  int64_t end_time = rdtsc();
  if (tree_type == MASSTREE) {
    std::cout<<"Masstree Scan + Insert: "<<end_time-start_time<<std::endl;
  } else if (tree_type == BTREE) {
    std::cout<<"Btree Scan + Insert: "<<end_time-start_time<<std::endl;
  } else if (tree_type == ART) {
    std::cout<<"Art Scan + Insert: "<<end_time-start_time<<std::endl;
  }
  //sleep(10);
}

TEST(TestMasstreeHybrid, Hybrid_test) {
  test_scan_insert(MASSTREE);
}

TEST(TestBtreeHybrid, Hybrid_test) {
  test_scan_insert(BTREE);
}

TEST(TestARTHybrid, Hybrid_test) {
  test_scan_insert(ART);
}

}
}

int main(int argc, char **argv)
{
  //oceanbase::unittest::BIND_CPU(pthread_self());
  oceanbase::common::ObLogger::get_logger().set_file_name("test_masstree.log", true);
  oceanbase::common::ObLogger::get_logger().set_log_level("INFO");
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

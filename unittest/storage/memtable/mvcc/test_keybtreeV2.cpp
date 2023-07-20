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

void generate_internal_keys(int *a, int n)
{
  // n = 15, example: 1,3,5,7,9,11,13,15,2,4,6,8,10,12,14
  int i = 0;
  for (int j = 1; j <= n; j += 2) {
    a[i++] = j;
  }
  for (int j = 2; j <= n; j += 2) {
    a[i++] = j;
  }
}

void generate_leafnode_keys(int *a, int n)
{
  // n = 15, example: 2,4,6,8,10,12,14,1,3,5,7,9,11,13,15
  int i = 0;
  for (int j = 2; j <= n; j += 2) {
    a[i++] = j;
  }
  for (int j = 1; j <= n; j += 2) {
    a[i++] = j;
  }
}

TEST(TestVersionSingleThread, smoke_test)
{
  Version v;
  for (int i = 1; i < 1000; i++) {
    v.increase_vinsert();
    ASSERT_EQ(v.get_vinsert(), i % 32);
    v.increase_vsplit();
    ASSERT_EQ(v.get_vsplit(), i);
    ASSERT_FALSE(v.is_inserting());
    ASSERT_FALSE(v.is_splitting());
    ASSERT_FALSE(v.is_latched());
  }
  v.set_inserting();
  ASSERT_TRUE(v.is_inserting());
  v.set_splitting();
  ASSERT_TRUE(v.is_splitting());
  v.latch();
  ASSERT_TRUE(v.is_latched());
  int vinsert = v.get_vinsert();
  int vsplit = v.get_vsplit();
  v.unlatch();
  ASSERT_FALSE(v.is_inserting());
  ASSERT_FALSE(v.is_splitting());
  ASSERT_FALSE(v.is_latched());
  ASSERT_EQ(v.get_vinsert(), (vinsert + 1) % 32);
  ASSERT_EQ(v.get_vsplit(), vsplit + 1);
}

TEST(TestVersionMultiThread, smoke_test)
{
  constexpr int BCOUNT = 16;
  class MockNode {
  public:
    int64_t a;
    int64_t b[BCOUNT];
    int64_t id;
  };

  constexpr int WRITE_THREAD_COUNT = 4;
  constexpr int READ_THREAD_COUNT = 8;
  constexpr int REPEAT_COUNT = 1000000;

  Version v;
  MockNode mock_node;
  mock_node.a = 0;
  for (int i = 0; i < BCOUNT; i++) {
    mock_node.b[i] = 0;
  }
  mock_node.id = -1;
  std::thread write_threads[WRITE_THREAD_COUNT];
  std::thread read_threads[READ_THREAD_COUNT];
  for (int i = 0; i < WRITE_THREAD_COUNT; i++) {
    write_threads[i] = std::thread(
        [&](int j) {
          for (int i = 0; i < REPEAT_COUNT; i++) {
            v.latch();
            ASSERT_EQ(mock_node.id, -1);
            v.set_splitting();
            mock_node.id = j;
            mock_node.a++;
            for (int i = 0; i < BCOUNT; i++) {
              mock_node.b[i]++;
            }
            mock_node.id = -1;
            v.unlatch();
            if (i % (WRITE_THREAD_COUNT * WRITE_THREAD_COUNT) == j * WRITE_THREAD_COUNT) {
              usleep(1);
            }
          }
        },
        i);
  }
  for (int i = 0; i < READ_THREAD_COUNT; i++) {
    read_threads[i] = std::thread(
        [&](int i) {
          int success_count = 0;
          for (int i = 0; i < REPEAT_COUNT; i++) {
            int64_t a, id;
            int64_t b[BCOUNT];
            Version snapshot_version;
            do {
              snapshot_version = v.get_stable_snapshot();
              a = mock_node.a;
              for (int i = 0; i < BCOUNT; i++) {
                b[i] = mock_node.b[i];
              }
              id = mock_node.id;
            } while (v.has_splitted(snapshot_version));
            ASSERT_EQ(id, -1);
            for (int i = 0; i < BCOUNT; i++) {
              ASSERT_EQ(b[i], a);
            }
            ASSERT_EQ(a, snapshot_version.get_vsplit());
          }
        },
        i);
  }
  for (int i = 0; i < WRITE_THREAD_COUNT; i++) {
    write_threads[i].join();
  }
  for (int i = 0; i < READ_THREAD_COUNT; i++) {
    read_threads[i].join();
  }
  for (int i = 0; i < BCOUNT; i++) {
    ASSERT_EQ(mock_node.b[i], mock_node.a);
  }
  ASSERT_EQ(mock_node.a, WRITE_THREAD_COUNT * REPEAT_COUNT);
  ASSERT_EQ(mock_node.id, -1);
}

TEST(TestPermutation, smoke_test)
{
  Permutation permutation;
  constexpr uint8_t n = Permutation::max_size();
  int data[n];
  generate_leafnode_keys(data, n);
  for (int len = 1; len <= n; len++) {
    int cur = data[len - 1];
    std::sort(data, data + len);
    int pos = 0;
    for (; pos < n; pos++) {
      if (data[pos] == cur) {
        break;
      }
    }
    permutation.insert(pos, cur);
    for (int i = 0; i < len; i++) {
      ASSERT_EQ(data[i], permutation.at(i));
    }
    ASSERT_EQ(len, permutation.size());
  }
  ASSERT_TRUE(permutation.is_full());

  permutation.set_size(n / 2);
  ASSERT_EQ(n / 2, permutation.size());
  ASSERT_FALSE(permutation.is_full());

  for (int i = 0; i < n / 2; i++) {
    ASSERT_EQ(data[i], permutation.at(i));
  }
}

using InternalNode = InternalNode<FakeKey, int64_t>;
using BtreeNode = BtreeNode<FakeKey, int64_t>;
using LeafNode = LeafNode<FakeKey, int64_t>;
using ObKeyBtree = ObKeyBtree<FakeKey, int64_t>;
using BtreeIterator = BtreeIterator<FakeKey, int64_t>;

TEST(TestInternalNode, smoke_test)
{
  constexpr int n = InternalNode::max_size();
  int data[n];

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  InternalNode *internal_node;
  ASSERT_EQ(node_allocator.make_internal(internal_node), OB_SUCCESS);

  std::vector<FakeKey> keys;

  FakeKey search_key = build_int_key(0);
  keys.push_back(search_key);

  generate_internal_keys(data, n);
  for (int i = 0; i < n; i++) {
    data[i] *= 2;
  }

  internal_node->set_leftmost_child(1);

  // test insert and search
  for (int len = 1; len <= n; len++) {
    int cur = data[len - 1];
    FakeKey key = build_int_key(cur);
    keys.push_back(key);
    internal_node->insert(key, cur);
    for (int i = 0; i < len * 2; i++) {
      int64_t val;
      search_key.set_int(data[i / 2] + i % 2);
      ASSERT_EQ(internal_node->search(search_key, val), OB_SUCCESS);
      ASSERT_EQ(data[i / 2] / 2 * 2, val);
    }
  }
  search_key.set_int(0);
  int64_t val;
  ASSERT_EQ(internal_node->search(search_key, val), OB_SUCCESS);
  ASSERT_EQ(val, 1);

  // test split
  std::sort(data, data + n);
  InternalNode *new_internal;
  ASSERT_EQ(node_allocator.make_internal(new_internal), OB_SUCCESS);
  new_internal->reset();
  FakeKey overflow_key = build_int_key(data[n - 1] + 1);
  FakeKey fence_key;
  keys.push_back(overflow_key);

  ASSERT_EQ(internal_node->split_and_insert(new_internal, overflow_key, data[n - 1] + 1, fence_key), OB_SUCCESS);
  ASSERT_EQ(fence_key.obj_->get_int(), data[n / 2]);

  // test left node
  ASSERT_EQ(internal_node->size(), n / 2);
  for (int i = 0; i < n / 2; i++) {
    search_key.set_int(data[i]);
    ASSERT_EQ(internal_node->search(search_key, val), OB_SUCCESS);
    ASSERT_EQ(data[i], val);
  }
  search_key.set_int(0);
  ASSERT_EQ(internal_node->search(search_key, val), OB_SUCCESS);
  ASSERT_EQ(val, 1);  // leftmost child

  // test right node
  ASSERT_EQ(new_internal->size(), (n + 1) / 2);
  for (int i = n / 2 + 1; i < n; i++) {
    search_key.set_int(data[i]);
    ASSERT_EQ(new_internal->search(search_key, val), OB_SUCCESS);
    ASSERT_EQ(data[i], val);
  }
  search_key.set_int(data[n - 1] + 1);
  ASSERT_EQ(new_internal->search(search_key, val), OB_SUCCESS);
  ASSERT_EQ(val, data[n - 1] + 1);
  ASSERT_EQ(new_internal->search(fence_key, val), OB_SUCCESS);
  ASSERT_EQ(val, fence_key.obj_->get_int());  // leftmost child

  for (int i = 0; i < keys.size(); i++) {
    allocator->free(keys[i].obj_);
  }
}

TEST(TestLeafNode, smoke_test)
{
  constexpr int n = LeafNode::max_size();
  int data[n];

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  LeafNode *leaf_node;
  ASSERT_EQ(node_allocator.make_leaf(leaf_node), OB_SUCCESS);

  std::vector<FakeKey> keys;

  FakeKey search_key = build_int_key(0);
  keys.push_back(search_key);

  generate_leafnode_keys(data, n);
  for (int i = 0; i < n; i++) {
    data[i] *= 2;
  }

  leaf_node->reset();

  // test insert and search
  for (int len = 1; len <= n; len++) {
    int cur = data[len - 1];
    FakeKey key = build_int_key(cur);
    keys.push_back(key);
    leaf_node->insert(key, cur);
    for (int i = 0; i < len * 2; i++) {
      int64_t val;
      search_key.set_int(data[i / 2] + i % 2);
      if (i % 2 == 0) {
        ASSERT_EQ(leaf_node->search(search_key, val), OB_SUCCESS);
        ASSERT_EQ(data[i / 2], val);
      } else {
        ASSERT_EQ(leaf_node->search(search_key, val), OB_ENTRY_NOT_EXIST);
      }
    }
  }
  search_key.set_int(0);
  int64_t val;
  ASSERT_EQ(leaf_node->search(search_key, val), OB_ENTRY_NOT_EXIST);

  // test split
  std::sort(data, data + n);
  LeafNode *new_leaf;
  ASSERT_EQ(node_allocator.make_leaf(new_leaf), OB_SUCCESS);
  new_leaf->reset();
  FakeKey overflow_key = build_int_key(data[n - 1] + 1);
  FakeKey fence_key;
  keys.push_back(overflow_key);

  ASSERT_EQ(leaf_node->split_and_insert(new_leaf, overflow_key, data[n - 1] + 1, fence_key), OB_SUCCESS);
  ASSERT_EQ(fence_key.obj_->get_int(), data[(n + 1) / 2]);

  // test left node
  ASSERT_EQ(leaf_node->size(), (n + 1) / 2);
  for (int i = 0; i < (n + 1) / 2; i++) {
    search_key.set_int(data[i]);
    ASSERT_EQ(leaf_node->search(search_key, val), OB_SUCCESS);
    ASSERT_EQ(data[i], val);
  }

  // test right node
  ASSERT_EQ(new_leaf->size(), n / 2 + 1);
  for (int i = (n + 1) / 2; i < n; i++) {
    search_key.set_int(data[i]);
    ASSERT_EQ(new_leaf->search(search_key, val), OB_SUCCESS);
    ASSERT_EQ(data[i], val);
  }
  search_key.set_int(data[n - 1] + 1);
  ASSERT_EQ(new_leaf->search(search_key, val), OB_SUCCESS);
  ASSERT_EQ(val, data[n - 1] + 1);

  for (int i = 0; i < keys.size(); i++) {
    allocator->free(keys[i].obj_);
  }
}

void judge_node_scan(LeafNode *leaf, FakeKey start_key, FakeKey end_key, bool include_start_key, bool include_end_key,
    bool is_backward, bool real_is_end, std::vector<int64_t> answer)
{
  KVQueue<FakeKey, int64_t> q;
  bool is_end;
  BtreeKV<FakeKey, int64_t> kv;
  int i = 0;

  ASSERT_EQ(leaf->scan(start_key, end_key, include_start_key, include_end_key, is_backward, q, is_end), OB_SUCCESS);
  ASSERT_EQ(is_end, real_is_end);
  ASSERT_EQ(answer.size(), q.size());
  while (q.size() > 0) {
    ASSERT_EQ(q.pop(kv), OB_SUCCESS);
    ASSERT_EQ(answer[i], kv.val_);
    i++;
  }
}

TEST(TestLeafNodeScan, smoke_test)
{
  constexpr int n = 7;
  int data[n] = {1, 3, 5, 7, 9, 11, 13};

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  LeafNode *leaf_node;
  ASSERT_EQ(node_allocator.make_leaf(leaf_node), OB_SUCCESS);

  std::vector<FakeKey> keys;

  FakeKey start_key = build_int_key(0);
  keys.push_back(start_key);

  FakeKey end_key = build_int_key(0);
  keys.push_back(end_key);

  for (int len = 1; len <= n; len++) {
    int cur = data[len - 1];
    FakeKey key = build_int_key(cur);
    keys.push_back(key);
    leaf_node->insert(key, cur);
  }

  std::vector<int64_t> ans;

  // [2, 10] -> 3,5,7,9
  ans = std::vector<int64_t>({3, 5, 7, 9});
  start_key.set_int(2);
  end_key.set_int(10);
  judge_node_scan(leaf_node, start_key, end_key, true, true, false, true, ans);

  // [5, 10] -> 5,7,9
  ans = std::vector<int64_t>({5, 7, 9});
  start_key.set_int(5);
  end_key.set_int(10);
  judge_node_scan(leaf_node, start_key, end_key, true, true, false, true, ans);

  // [2, 11] -> 3,5,7,9,11
  ans = std::vector<int64_t>({3, 5, 7, 9, 11});
  start_key.set_int(2);
  end_key.set_int(11);
  judge_node_scan(leaf_node, start_key, end_key, true, true, false, true, ans);

  // (3, 11) -> 5,7,9
  ans = std::vector<int64_t>({5, 7, 9});
  start_key.set_int(3);
  end_key.set_int(11);
  judge_node_scan(leaf_node, start_key, end_key, false, false, false, true, ans);

  // (2, 7) -> 3,5
  ans = std::vector<int64_t>({3, 5});
  start_key.set_int(2);
  end_key.set_int(7);
  judge_node_scan(leaf_node, start_key, end_key, false, false, false, true, ans);

  // (3, 8) -> 5,7
  ans = std::vector<int64_t>({5, 7});
  start_key.set_int(3);
  end_key.set_int(8);
  judge_node_scan(leaf_node, start_key, end_key, false, false, false, true, ans);

  // [3, 8) -> 3,5,7
  ans = std::vector<int64_t>({3, 5, 7});
  start_key.set_int(3);
  end_key.set_int(8);
  judge_node_scan(leaf_node, start_key, end_key, true, false, false, true, ans);

  // (3, 9] -> 5,7,9
  ans = std::vector<int64_t>({5, 7, 9});
  start_key.set_int(3);
  end_key.set_int(9);
  judge_node_scan(leaf_node, start_key, end_key, false, true, false, true, ans);

  // [9, 13] -> 9,11,13
  ans = std::vector<int64_t>({9, 11, 13});
  start_key.set_int(9);
  end_key.set_int(13);
  judge_node_scan(leaf_node, start_key, end_key, true, true, false, true, ans);

  // [0, 14] -> 1,3,5,7,9,11,13
  ans = std::vector<int64_t>({1, 3, 5, 7, 9, 11, 13});
  start_key.set_int(0);
  end_key.set_int(14);
  judge_node_scan(leaf_node, start_key, end_key, true, true, false, false, ans);

  // [9, 13) -> 9,11
  ans = std::vector<int64_t>({9, 11});
  start_key.set_int(9);
  end_key.set_int(13);
  judge_node_scan(leaf_node, start_key, end_key, true, false, false, true, ans);

  // [9, 14) -> 9,11,13
  ans = std::vector<int64_t>({9, 11, 13});
  start_key.set_int(9);
  end_key.set_int(14);
  judge_node_scan(leaf_node, start_key, end_key, true, false, false, false, ans);

  // (13, 15) -> 
  ans = std::vector<int64_t>({});
  start_key.set_int(13);
  end_key.set_int(15);
  judge_node_scan(leaf_node, start_key, end_key, false, false, false, false, ans);

  // [100, 200] -> 
  ans = std::vector<int64_t>({});
  start_key.set_int(100);
  end_key.set_int(200);
  judge_node_scan(leaf_node, start_key, end_key, true, true, false, false, ans);

  // [10,2] -> 9,7,5,3
  ans = std::vector<int64_t>({9, 7, 5, 3});
  start_key.set_int(10);
  end_key.set_int(2);
  judge_node_scan(leaf_node, start_key, end_key, true, true, true, true, ans);

  // [11,1] -> 11,9,7,5,3,1
  ans = std::vector<int64_t>({11, 9, 7, 5, 3, 1});
  start_key.set_int(11);
  end_key.set_int(1);
  judge_node_scan(leaf_node, start_key, end_key, true, true, true, true, ans);

  // (11,1) -> 9,7,5,3
  ans = std::vector<int64_t>({9, 7, 5, 3});
  start_key.set_int(11);
  end_key.set_int(1);
  judge_node_scan(leaf_node, start_key, end_key, false, false, true, true, ans);

  // (11,1] -> 9,7,5,3,1
  ans = std::vector<int64_t>({9, 7, 5, 3, 1});
  start_key.set_int(11);
  end_key.set_int(1);
  judge_node_scan(leaf_node, start_key, end_key, false, true, true, true, ans);

  // [13,0] -> 13,11,9,7,5,3,1
  ans = std::vector<int64_t>({13, 11, 9, 7, 5, 3, 1});
  start_key.set_int(13);
  end_key.set_int(0);
  judge_node_scan(leaf_node, start_key, end_key, true, true, true, false, ans);

  // [13,0) -> 13,11,9,7,5,3,1
  ans = std::vector<int64_t>({13, 11, 9, 7, 5, 3, 1});
  start_key.set_int(13);
  end_key.set_int(0);
  judge_node_scan(leaf_node, start_key, end_key, true, false, true, false, ans);
  
  // [-100,-50] -> 
  ans = std::vector<int64_t>({});
  start_key.set_int(-100);
  end_key.set_int(-50);
  judge_node_scan(leaf_node, start_key, end_key, true, true, true, false, ans);

  for (int i = 0; i < keys.size(); i++) {
    allocator->free(keys[i].obj_);
  }
}

void judge_tree_scan(ObKeyBtree *btree, FakeKey start_key, FakeKey end_key, bool include_start_key,
    bool include_end_key, bool is_backward, std::vector<int64_t> &answer)
{
  FakeKey key;
  int64_t val;
  int i = 0;

  BtreeIterator iter(btree, start_key, end_key, include_start_key, include_end_key, is_backward);
  iter.init();
  while (iter.iter_next(key, val) == OB_SUCCESS) {
    ASSERT_EQ(val, answer[i]);
    i++;
  }
  if(i != answer.size()) {
    std::cout<<start_key.get_ptr()->get_int()<<" "<<end_key.get_ptr()->get_int()<<std::endl;
    std::cout<<include_start_key<<" "<<include_end_key<<" "<<is_backward<<std::endl;
    for(int i=0;i<answer.size();i++) {
      std::cout<<answer[i]<<" ";
    }
    std::cout<<std::endl;
  }
  ASSERT_EQ(i, answer.size());
}

void free_btree(ObKeyBtree &btree)
{
  FakeKey start_key = build_int_key(INT64_MIN);
  FakeKey end_key = build_int_key(INT64_MAX);
  FakeKey key;
  int64_t val;
  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeIterator iter(&btree, start_key, end_key, true, true, false);
  iter.init();
  while (iter.iter_next(key, val) == OB_SUCCESS) {
    allocator->free(key.get_ptr());
  }
}

TEST(TestBtree, smoke_test)
{
  constexpr int KEY_NUM = 100000;
  std::vector<int64_t> data(KEY_NUM);

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  ObKeyBtree btree(node_allocator);

  ASSERT_EQ(btree.init(), OB_SUCCESS);

  FakeKey search_key = build_int_key(0);

  for (int i = 0; i < KEY_NUM; i++) {
    data[i] = i * 2;
  }
  std::random_shuffle(data.begin(), data.end());

  // test insert and search
  for (int len = 1; len <= KEY_NUM; len++) {
    int cur = data[len - 1];
    FakeKey key = build_int_key(cur);
    btree.insert(key, cur);
    if (len % (KEY_NUM / 49) == 0) {
      for (int i = 0; i < len * 2; i++) {
        int64_t val;
        search_key.set_int(data[i / 2] + i % 2);
        if (i % 2 == 0) {
          ASSERT_EQ(btree.search(search_key, val), OB_SUCCESS);
          ASSERT_EQ(data[i / 2], val);
        } else {
          ASSERT_EQ(btree.search(search_key, val), OB_ENTRY_NOT_EXIST);
        }
      }
    }
  }
  std::sort(data.begin(), data.end());
  search_key.set_int(-1);
  int64_t val;
  ASSERT_EQ(btree.search(search_key, val), OB_ENTRY_NOT_EXIST);

  FakeKey start_key = build_int_key(0);
  FakeKey end_key = build_int_key(0);

  // test scan
  int REPEAT_COUNT;

  // forward include
  REPEAT_COUNT = 100;
  while (REPEAT_COUNT--) {
    int64_t start_int = ObRandom::rand(-KEY_NUM, KEY_NUM*3);
    int64_t end_int = ObRandom::rand(start_int + 1, KEY_NUM*3);
    start_key.set_int(start_int);
    end_key.set_int(end_int);
    std::vector<int64_t> ans;
    for (int i = max(0, (start_int+1)/2*2); i <= min((KEY_NUM-1)*2, end_int/2*2); i+=2) {
      ans.push_back(i);
    }
    judge_tree_scan(&btree, start_key, end_key, true, true, false, ans);
  }

  // forward exclude
  REPEAT_COUNT = 100;
  while (REPEAT_COUNT--) {
    int64_t start_int = ObRandom::rand(-KEY_NUM, KEY_NUM*3);
    int64_t end_int = ObRandom::rand(start_int + 1, KEY_NUM*3);
    start_key.set_int(start_int);
    end_key.set_int(end_int);
    std::vector<int64_t> ans;
    for (int i = max(0, (start_int+2)/2*2); i <= min((KEY_NUM-1)*2, (end_int-1)/2*2); i+=2) {
      ans.push_back(i);
    }
    judge_tree_scan(&btree, start_key, end_key, false, false, false, ans);
  }

  // backward include
  REPEAT_COUNT = 100;
  while (REPEAT_COUNT--) {
    int64_t start_int = ObRandom::rand(-KEY_NUM, KEY_NUM*3);
    int64_t end_int = ObRandom::rand(start_int + 1, KEY_NUM*3);
    start_key.set_int(start_int);
    end_key.set_int(end_int);
    std::vector<int64_t> ans;
    for (int i = min((KEY_NUM-1)*2, end_int/2*2); i >= max(0, (start_int+1)/2*2); i-=2) {
      ans.push_back(i);
    }
    judge_tree_scan(&btree, end_key, start_key, true, true, true, ans);
  }

  // backward exclude
  REPEAT_COUNT = 100;
  while (REPEAT_COUNT--) {
    int64_t start_int = ObRandom::rand(-KEY_NUM, KEY_NUM*3);
    int64_t end_int = ObRandom::rand(start_int + 1, KEY_NUM*3);
    start_key.set_int(start_int);
    end_key.set_int(end_int);
    std::vector<int64_t> ans;
    for (int i = min((KEY_NUM-1)*2, (end_int-1)/2*2); i >= max(0, (start_int+2)/2*2); i-=2) {
      ans.push_back(i);
    }
    judge_tree_scan(&btree, end_key, start_key, false, false, true, ans);
  }

  free_btree(btree);
  allocator->free(search_key.get_ptr());
  allocator->free(start_key.get_ptr());
  allocator->free(end_key.get_ptr());
}

TEST(TestEventualConsistency, smoke_test)
{
  constexpr uint64_t KEY_NUM = 3000000;
  constexpr uint64_t THREAD_COUNT = 10;
  constexpr uint64_t PER_THREAD_INSERT_COUNT = KEY_NUM / THREAD_COUNT;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  ObKeyBtree btree(node_allocator);

  ASSERT_EQ(btree.init(), OB_SUCCESS);

  std::thread threads[THREAD_COUNT];
  std::vector<std::vector<int>> data(THREAD_COUNT, std::vector<int>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = THREAD_COUNT * j + i;
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  // concurrent insert
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id] = std::thread(
        [&](int i) {
          for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
            btree.insert(build_int_key(data[i][j]), data[i][j]);
          }
        },
        thread_id);
  }

  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id].join();
  }

  // evaluate the tree
  FakeKey start_key = build_int_key(0);
  FakeKey end_key = build_int_key(KEY_NUM);
  FakeKey key;
  int64_t val;

  BtreeIterator iter(&btree, start_key, end_key, true, true, false);
  iter.init();
  int i = 0;
  while (iter.iter_next(key, val) == OB_SUCCESS) {
    ASSERT_EQ(val, i);
    i++;
  }
  ASSERT_EQ(i, KEY_NUM);

  free_btree(btree);
  allocator->free(start_key.get_ptr());
  allocator->free(end_key.get_ptr());
}

TEST(TestMonotonicReadWrite, smoke_test)
{
  constexpr int KEY_NUM = 3000000;
  constexpr int WRITE_THREAD_COUNT = 10;
  constexpr int PER_THREAD_INSERT_COUNT = KEY_NUM / WRITE_THREAD_COUNT;
  constexpr int SCAN_THREAD_COUNT = 10;
  constexpr int PER_THREAD_SCAN_COUNT = 8;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  ObKeyBtree btree(node_allocator);

  ASSERT_EQ(btree.init(), OB_SUCCESS);

  // constructing insert keys
  std::vector<std::vector<int64_t>> data(WRITE_THREAD_COUNT, std::vector<int64_t>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < WRITE_THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = WRITE_THREAD_COUNT * j + i;
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  std::thread write_threads[WRITE_THREAD_COUNT];
  for (int thread_id = 0; thread_id < WRITE_THREAD_COUNT; thread_id++) {
    write_threads[thread_id] = std::thread(
        [&](int i) {
          // insert in order
          for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
            btree.insert(build_int_key(data[i][j]), data[i][j]);
            usleep(1);
          }
        },
        thread_id);
  }

  std::thread scan_threads[SCAN_THREAD_COUNT];
  for (int thread_id = 0; thread_id < SCAN_THREAD_COUNT; thread_id++) {
    scan_threads[thread_id] = std::thread(
        [&](int thread_id) {
          FakeKey start_key = build_int_key(-1);
          FakeKey end_key = build_int_key(KEY_NUM + 1);
          int scan_count = PER_THREAD_SCAN_COUNT;
          std::unordered_set<int64_t> last_results;
          while (scan_count--) {
            std::unordered_set<int64_t> results;
            FakeKey key;
            int64_t val;
            if (thread_id % 2 == 0) {
              // scan forward
              BtreeIterator iter(&btree, start_key, end_key, true, true, false);
              iter.init();
              int64_t last = -1;
              while (iter.iter_next(key, val) == OB_SUCCESS) {
                results.insert(val);
                ASSERT_GT(val, last);
                last = val;
              }
            } else {
              // scan backward
              BtreeIterator iter(&btree, end_key, start_key, true, true, true);
              iter.init();
              int64_t last = KEY_NUM + 1;
              while (iter.iter_next(key, val) == OB_SUCCESS) {
                results.insert(val);
                ASSERT_LT(val, last);
                last = val;
              }
            }
            // test monotonic write, if a thread see a key A, then it should see all keys inserted before A
            for (int i = 0; i < WRITE_THREAD_COUNT; i++) {
              if (thread_id % 2 == 0) {
                int64_t min = KEY_NUM + 1;
                for (int j = PER_THREAD_INSERT_COUNT - 1; j >= 0; j--) {
                  ASSERT_TRUE(data[i][j] < min || results.count(data[i][j]) == 1);
                  if (results.count(data[i][j]) == 1) {
                    min = std::min(min, data[i][j]);
                  }
                }
              } else {
                int64_t max = -1;
                for (int j = PER_THREAD_INSERT_COUNT - 1; j >= 0; j--) {
                  ASSERT_TRUE(data[i][j] > max || results.count(data[i][j]) == 1);
                  if (results.count(data[i][j]) == 1) {
                    max = std::max(max, data[i][j]);
                  }
                }
              }
            }
            // test monotonic read, if a thread do two scan, then the frist scan result should be the subset of
            // the second scan result.
            for (auto i : last_results) {
              ASSERT_TRUE(results.count(i) == 1);
            }
            last_results = results;
          }
          allocator->free(start_key.get_ptr());
          allocator->free(end_key.get_ptr());
        },
        thread_id);
  }

  for (int thread_id = 0; thread_id < WRITE_THREAD_COUNT; thread_id++) {
    write_threads[thread_id].join();
  }
  for (int thread_id = 0; thread_id < SCAN_THREAD_COUNT; thread_id++) {
    scan_threads[thread_id].join();
  }

  free_btree(btree);
}

TEST(TestSequentialConsistency, smoke_test)
{
  constexpr int PER_THREAD_INSERT_COUNT = 200000;
  constexpr int READ_THREAD_COUNT = 10;

  int progress = -1;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  BtreeNodeAllocator<FakeKey, int64_t> node_allocator(*allocator);
  ObKeyBtree btree(node_allocator);
  ASSERT_EQ(btree.init(), OB_SUCCESS);

  std::thread main_thread([&] {
    for (; progress < PER_THREAD_INSERT_COUNT; progress++) {
      usleep(1);
    }
  });

  std::vector<int64_t> insert_keys(PER_THREAD_INSERT_COUNT * 2);
  for (int i = 0; i < insert_keys.size(); i++) {
    insert_keys[i] = i;
  }
  std::random_shuffle(insert_keys.begin(), insert_keys.end());

  std::thread write_threads[2];
  for (int thread_id = 0; thread_id < 2; thread_id++) {
    write_threads[thread_id] = std::thread(
        [&](int thread_id) {
          int last = -1;
          int insert_id = 0;
          while (last < PER_THREAD_INSERT_COUNT) {
            while (last >= progress) {}
            last++;
            insert_id = last * 2 + thread_id;
            btree.insert(build_int_key(insert_keys[insert_id]), insert_keys[insert_id]);
          }
        },
        thread_id);
  }

  std::vector<std::vector<bool>> read_results(READ_THREAD_COUNT, std::vector<bool>(PER_THREAD_INSERT_COUNT));

  std::thread read_threads[READ_THREAD_COUNT];
  for (int thread_id = 0; thread_id < READ_THREAD_COUNT; thread_id++) {
    read_threads[thread_id] = std::thread(
        [&](int thread_id) {
          FakeKey search_key1 = build_int_key(0);
          FakeKey search_key2 = build_int_key(0);
          int64_t val;
          for (int i = 0; i < PER_THREAD_INSERT_COUNT; i++) {
            search_key1.set_int(insert_keys[i * 2]);
            search_key2.set_int(insert_keys[i * 2 + 1]);
            if (thread_id % 2 == 0) {
              while (btree.search(search_key1, val) != OB_SUCCESS) {}
              if (btree.search(search_key2, val) == OB_ENTRY_NOT_EXIST) {
                // the order this thread sees is: search_key1 -> search_key2
                read_results[thread_id][i] = true;
              }
            } else {
              while (btree.search(search_key2, val) != OB_SUCCESS) {}
              if (btree.search(search_key1, val) == OB_ENTRY_NOT_EXIST) {
                // the order this thread sees is: search_key2 -> search_key1
                read_results[thread_id][i] = true;
              }
            }
          }
          allocator->free(search_key1.get_ptr());
          allocator->free(search_key2.get_ptr());
        },
        thread_id);
  }

  main_thread.join();
  write_threads[0].join();
  write_threads[1].join();
  for (int i = 0; i < READ_THREAD_COUNT; i++) {
    read_threads[i].join();
  }

  int count = 0;
  for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
    for (int i = 0; i < READ_THREAD_COUNT; i++) {
      read_results[i % 2][j] = read_results[i % 2][j] || read_results[i][j];
    }
    // threads shouldn't see different order
    ASSERT_FALSE(read_results[0][j] && read_results[1][j]);
  }

  free_btree(btree);
}

TEST(TestConcurrency2, smoke_test)
{
  constexpr uint64_t KEY_NUM = 3000000;
  constexpr uint64_t THREAD_COUNT = 10;
  constexpr uint64_t PER_THREAD_INSERT_COUNT = KEY_NUM / THREAD_COUNT;

  FakeAllocator *allocator = FakeAllocator::get_instance();
  keybtree::BtreeNodeAllocator<FakeKey, int64_t *> node_allocator(*allocator);
  keybtree::ObKeyBtree<FakeKey, int64_t *> btree(node_allocator);

  ASSERT_EQ(btree.init(), OB_SUCCESS);

  std::thread threads[THREAD_COUNT];
  std::vector<std::vector<int64_t>> data(THREAD_COUNT, std::vector<int64_t>(PER_THREAD_INSERT_COUNT));
  for (int i = 0; i < THREAD_COUNT; i++) {
    for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
      data[i][j] = THREAD_COUNT * j + i;
    }
    std::random_shuffle(data[i].begin(), data[i].end());
  }

  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id] = std::thread(
        [&](int i) {
          int64_t *val = nullptr;
          for (int j = 0; j < PER_THREAD_INSERT_COUNT; j++) {
            btree.insert(build_int_key(data[i][j]), val);
          }
        },
        thread_id);
  }
  uint64_t start_rdtsc = rdtsc();
  for (int thread_id = 0; thread_id < THREAD_COUNT; thread_id++) {
    threads[thread_id].join();
  }
  std::cout << rdtsc() - start_rdtsc << std::endl;
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
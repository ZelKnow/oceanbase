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

#ifndef __OB_KEYBTREEV2_H__
#define __OB_KEYBTREEV2_H__

#include "lib/metrics/ob_counter.h"
#include "lib/allocator/ob_allocator.h"
#include "lib/utility/ob_print_utils.h"
#include "lib/utility/utility.h"

namespace oceanbase {

namespace keybtreeV2 {

template <typename BtreeKey, typename BtreeVal>
struct BtreeKV {
  BtreeKey key_;  // 8byte
  BtreeVal val_;  // 8byte
};

/**
 * Version is used for BtreeNode's read/write synchronization. Readers does not need to
 * latch or modify any shared variables while reading, instead they need to double check
 * the version to make sure that no writers have modified the node while reading.
 *
 * Version layout:
 *  _______________________________________
 * | | | |           |                     |
 * | | | |  vinsert  |        vsplit       |
 * |_|_|_|___________|_____________________|
 *  1 1 1      5                56
 * the first three bits are: latch, inserting, splitting.
 *
 * Reader path:
 * 1) Get stable version snapshot. (inserting/splitting == 0, meaning no writer is modifying)
 * 2) Read the node.
 * 3) Double check version (has_inserted(), has_splitted()).
 *
 * Writer path:
 * 1) Latch. Only one writer can modify the node at the same time.
 * 2) Set the inserting/splitting bit(to tell the readers that I'm modifying the node.)
 * 3) Modify the node.
 * 4) Increase vinsert/vsplit(to tell the readers a modification has occurred)
 * 5) Unset the inserting/splitting bit and unlatch.
 */
class Version {
public:
  Version() : data_(0)
  {}
  Version get_stable_snapshot()
  {
    Version stable_version;
    do {
      stable_version.data_ = ATOMIC_LOAD_ACQ(&this->data_);
    } while (stable_version.data_ & DIRTY_MASK);
    // As a reader, we need to ensure that all our read operations on the node
    // occur after we take a snapshot of the version, so a acquire fence is needed.
    __atomic_thread_fence(__ATOMIC_ACQUIRE);
    return stable_version;
  }
  OB_INLINE uint64_t get_vsplit() const
  {
    return (data_ & SPLIT_MASK) >> 0;
  }
  OB_INLINE uint64_t get_vinsert() const
  {
    return (data_ & INSERT_MASK) >> 56;
  }
  OB_INLINE bool is_inserting() const
  {
    return data_ & INSERTING_BIT;
  }
  OB_INLINE bool is_splitting() const
  {
    return data_ & SPLITTING_BIT;
  }
  OB_INLINE bool is_latched() const
  {
    return data_ & LATCH_BIT;
  }
  OB_INLINE void increase_vinsert()
  {
    data_ = (data_ & (~INSERT_MASK)) + (((data_ & INSERT_MASK) + VINSERT_LOWBIT) & INSERT_MASK);
  }
  OB_INLINE void increase_vsplit()
  {
    data_ = (data_ & (~SPLIT_MASK)) + (((data_ & SPLIT_MASK) + VSPLIT_LOWBIT) & SPLIT_MASK);
  }
  bool has_splitted(Version &snapshot_version) const
  {
    // We need to double check the version to make sure that all our reads on node
    // are consistent, and we need to ensure all our reads occur before our double check.
    // So a acquire fence is needed here.
    __atomic_thread_fence(__ATOMIC_ACQUIRE);
    if (is_splitting()) {
      return true;
    }
    if (get_vsplit() != snapshot_version.get_vsplit()) {
      return true;
    }
    return false;
  }
  bool has_inserted(Version &snapshot_version) const
  {
    __atomic_thread_fence(__ATOMIC_ACQUIRE);
    if (is_inserting()) {
      return true;
    }
    if (get_vinsert() != snapshot_version.get_vinsert()) {
      return true;
    }
    return false;
  }
  void set_inserting()
  {
    data_ = data_ | INSERTING_BIT;
    // As a writer, we need to make sure that all our write operations to the node happen after
    // we set the inserting/splitting bit, otherwise the reader may not find themselves
    // having read an inconsistent state. So a release fence is needed.
    __atomic_thread_fence(__ATOMIC_RELEASE);
  }
  void set_splitting()
  {
    data_ = data_ | SPLITTING_BIT;
    __atomic_thread_fence(__ATOMIC_RELEASE);
  }
  void latch()
  {
    while (true) {
      uint64_t old_data = data_ & (~LATCH_BIT);
      uint64_t new_data = old_data | LATCH_BIT;
      if (ATOMIC_BCAS(&data_, old_data, new_data)) {
        return;
      }
    }
    // We don't need to set fence here because CAS do that for us.
  }
  void unlatch()
  {
    if (is_inserting()) {
      increase_vinsert();
    }
    if (is_splitting()) {
      increase_vsplit();
    }
    // We need to make sure all our writes to the node happen before we unlatch the node
    // Otherwise readers may not find themselves having read an inconsistent state. So
    // a release fence is needed here.
    // __atomic_thread_fence(__ATOMIC_SEQ_CST);
    ATOMIC_BCAS(&data_, data_, data_ & ULATCH_MASK);
  }
  void reset()
  {
    data_ = 0;
  }

private:
  uint64_t data_;
  enum {
    LATCH_BIT = (1ULL << 63),
    INSERTING_BIT = (1ULL << 62),
    SPLITTING_BIT = (1ULL << 61),
    VINSERT_LOWBIT = (1ULL << 56),
    VSPLIT_LOWBIT = (1ULL << 0),
    SPLIT_MASK = ~(~0ULL << 56) << 0,
    INSERT_MASK = ~(~0ULL << 5) << 56,
    ULATCH_MASK = ~(LATCH_BIT | INSERTING_BIT | SPLITTING_BIT),
    DIRTY_MASK = INSERTING_BIT | SPLITTING_BIT
  };
};

/**
 * Permutation represents the correct key order in BtreeNode, plus the current size of the node. Readers
 * will get the real position of the key through Permutation. The two advantages of using it are:
 * 1) To avoid key rearrangement. If a writer needs to insert a key-value pair into the node, it can just
 *    append it to the key-value array, and then modify the Permutation feild to make it visible to readers.
 *    Because Permutation is 64-bit len, we can use bit operation to achieve it efficiently.
 * 2) To avoid forcing readers to retry while doing insert in the LeafNode. The reader can snapshot the
 *    permutation to ensure that it sees consistent state during searching on leaf node, so that as long
 *    as no splits have occurred, the reader can safely assume that the value it get is correct.
 *
 * This 64-bit Permutation can be divided into 16 4-bit sub fields. The lowest 4 bits represent the size
 * of the Node, and the remaining are 15 4-bit integers that represent the actual position of the corresponding
 * key in the correct order.
 *
 * Note that the methods of Permutation DO NOT check the validity of the parameters(e.g. pos should be less than
 * max size)
 */
class Permutation {
public:
  Permutation() : data_(0)
  {}
  OB_INLINE uint8_t size() const
  {
    return data_ & COUNTER_MASK;
  }
  OB_INLINE static constexpr uint8_t max_size()
  {
    return MAX_SIZE;
  }
  OB_INLINE uint8_t at(uint8_t pos) const
  {
    return (data_ >> (ENTRY_SIZE * pos + COUNTER_SIZE)) & ENTRY_MASK;
  }
  OB_INLINE void insert(uint8_t pos, uint8_t value)
  {
    uint8_t offset = pos * ENTRY_SIZE + COUNTER_SIZE;
    data_ =
        (((data_ & (~0ULL << offset)) << ENTRY_SIZE) | (uint64_t(value) << offset) | (data_ & ((1ULL << offset) - 1))) +
        1;
  }
  // will not change the size
  OB_INLINE void set(uint8_t pos, uint8_t value)
  {
    uint8_t offset = pos * ENTRY_SIZE + COUNTER_SIZE;
    data_ =
        (((data_ & (~0ULL << (offset + ENTRY_SIZE)))) | (uint64_t(value) << offset) | (data_ & ((1ULL << offset) - 1)));
  }
  OB_INLINE void set_size(uint8_t size)
  {
    data_ = (data_ & (~COUNTER_MASK)) + size;
  }
  OB_INLINE void reset()
  {
    data_ = 0;
  }
  OB_INLINE bool is_full() const
  {
    return size() >= MAX_SIZE;
  }

private:
  uint64_t data_;
  enum {
    ENTRY_SIZE = 4,
    COUNTER_SIZE = 4,
    MAX_SIZE = (sizeof(data_) * 8 - COUNTER_SIZE) / ENTRY_SIZE,
    COUNTER_MASK = (1ULL << COUNTER_SIZE) - 1,
    ENTRY_MASK = (1ULL << ENTRY_SIZE) - 1,
  };
};

template <typename BtreeKey, typename BtreeVal>
class BtreeNode {
private:
  using BtreeKV = BtreeKV<BtreeKey, BtreeVal>;

public:
  BtreeNode() : version_(), permutation_()
  {}
  void reset()
  {
    version_.reset();
    permutation_.reset();
  }
  void dump(FILE *file);
  virtual uint8_t get_level() const = 0;
  OB_INLINE Version &get_version()
  {
    return version_;
  }
  OB_INLINE bool is_full() const
  {
    return permutation_.is_full();
  }
  OB_INLINE static constexpr int max_size()
  {
    return Permutation::max_size();
  }
  OB_INLINE int size()
  {
    return permutation_.size();
  }
  OB_INLINE BtreeKV &get_kv(const int pos, const Permutation &snapshot_permutation)
  {
    int real_pos = snapshot_permutation.at(pos);
    return kvs_[real_pos];
  }
  OB_INLINE Permutation get_permutation() const
  {
    return permutation_;
  }
  /**
   * @brief Insert \p key and \p val into the node. Doesn't check if there is enough space for insert.
   * @param key the key to be inserted.
   * @param val the value to be inserted.
   * @return OB_SUCCESS on success, others on compare fail.
   */
  int insert(BtreeKey key, BtreeVal val);
  virtual int search(const BtreeKey key, BtreeVal &val) = 0;
  /**
   * @brief Split \p this node, move right-half keys to \p new_node, and then insert \p key and \p val
   *        into \p this or \p new_node.
   * @param new_node the node that split off. Requires allocation and initialization before invoking.
   * @param key the key to be inserted.
   * @param val the value to be inserted.
   * @param[out] fence_key the fence key to be inserted into the parent node.
   * @return OB_SUCCESS on success, others on compare fail.
   */
  virtual int split_and_insert(BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key) = 0;

protected:
  /**
   * @brief Binary search the node to find the position of the first key that is
   *        greater than or equal to \p key.
   * @param key the search key
   * @param snapshot_permutation the permutation used to search
   * @return -1 if \p key is greater than all key in the node.
   */
  int search_(const BtreeKey key, const Permutation &snapshot_permutation, int &pos);
  /**
   * @brief Copy [\p start, \p end] key-value pairs to \p other.
   */
  OB_INLINE void copy_to_(BtreeNode *other, int start, int end);
  Version version_;
  Permutation permutation_;
  BtreeKV kvs_[max_size()];
  DISALLOW_COPY_AND_ASSIGN(BtreeNode);
};

// This is used to store the batch scan results. Can only store one node's data.
template <typename BtreeKey, typename BtreeVal>
class KVQueue {
private:
  using BtreeKV = BtreeKV<BtreeKey, BtreeVal>;
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  enum { capacity = BtreeNode::max_size() + 1 };

public:
  KVQueue() : push_(0), pop_(0)
  {}
  ~KVQueue()
  {}
  OB_INLINE void reset()
  {
    push_ = 0;
    pop_ = 0;
  }
  OB_INLINE int push(const BtreeKV &data)
  {
    int ret = OB_SUCCESS;
    if (push_ >= pop_ + capacity) {
      ret = OB_ARRAY_OUT_OF_RANGE;
    } else {
      items_[idx(push_++)] = data;
    }
    return ret;
  }
  OB_INLINE int pop(BtreeKV &data)
  {
    int ret = OB_SUCCESS;
    if (pop_ >= push_) {
      ret = OB_ARRAY_OUT_OF_RANGE;
    } else {
      data = items_[idx(pop_++)];
    }
    return ret;
  }
  OB_INLINE int top(BtreeKV &data)
  {
    int ret = OB_SUCCESS;
    if (pop_ >= push_) {
      ret = OB_ARRAY_OUT_OF_RANGE;
    } else {
      data = items_[idx(pop_)];
    }
    return ret;
  }
  OB_INLINE int64_t size() const
  {
    return push_ - pop_;
  }
  OB_INLINE bool empty() const
  {
    return push_ == pop_;
  }

private:
  int64_t idx(const int64_t x)
  {
    return x % capacity;
  }
  int64_t push_;
  int64_t pop_;
  BtreeKV items_[capacity];
};

template <typename BtreeKey, typename BtreeVal>
class LeafNode : public BtreeNode<BtreeKey, BtreeVal> {
private:
  using KVQueue = KVQueue<BtreeKey, BtreeVal>;
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  using BtreeKV = BtreeKV<BtreeKey, BtreeVal>;

public:
  LeafNode() : BtreeNode(), prev_(nullptr), next_(nullptr)
  {}
  /**
   * @brief Scan the leaf node, push the key-value pair in the interval [ \p start_key, \p end_key]
            to \p kv_queue.
   *
   * @param start_key the start key of the scan
   * @param end_key the end key of the scan
   * @param include_start_key true if the start key should be included in the result
                              (if the key is in the node)
   * @param include_end_key true if the start key should be included in the scan result
                            (if the key is in the node)
   * @param is_backward true if it is a backward scan
   * @param[out] kv_queue the scan result queue
   * @param[out] is_end true if there are no more keys left to scan
   * @return int OB_SUCCESS on success, others on failure
   */
  int scan(BtreeKey start_key, BtreeKey end_key, bool include_start_key, bool include_end_key, bool is_backward,
      KVQueue &kv_queue, bool &is_end);
  /**
   * @brief Scan the leaf node, push the key-value pair in the interval [MIN_KEY, \p end_key]
            to \p kv_queue.
   *
   * @param end_key the end key of the scan
   * @param include_end_key true if the start key should be included in the scan result
                            (if the key is in the node)
   * @param is_backward true if it is a backward scan
   * @param[out] kv_queue the scan result queue
   * @param[out] is_end true if there are no more keys left to scan
   * @return int OB_SUCCESS on success, others on failure
   */
  int scan(BtreeKey end_key, bool include_end_key, bool is_backward, KVQueue &kv_queue, bool &is_end);
  // Does not modify perv_ and next_ pointer.
  virtual int split_and_insert(BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key) override;
  /**
   * @brief Search \p key on the LeafNode.
   * @param key the search key.
   * @param[out] val the corresponding value.
   * @return OB_SUCCESS on success, OB_ENTRY_NOT_EXIST if key is not exist, others on compare fail.
   */
  virtual int search(const BtreeKey key, BtreeVal &val) override;
  OB_INLINE LeafNode *&get_prev()
  {
    return prev_;
  }
  OB_INLINE LeafNode *&get_next()
  {
    return next_;
  }
  OB_INLINE void set_prev(LeafNode *prev)
  {
    prev_ = prev;
  }
  OB_INLINE void set_next(LeafNode *next)
  {
    next_ = next;
  }
  virtual uint8_t get_level() const override
  {
    return 0;
  }

private:
  LeafNode *prev_, *next_;
  /**
   * @brief Scan the leaf node from \p start_pos to \p end_pos, push the key-value to \p kv_queue one by one.
   *
   * @return OB_SUCCESS on success, others on KVQueue failure.
   */
  int scan_(int start_pos, int end_pos, bool is_backward, KVQueue &kv_queue, Permutation snapshot_permutation);
  /**
   * @brief find the leftmost postion \p boundary_pos, such that \p key is greater than
   * (or equal to, depending on the \p included parameter) the key at \p boundary_pos.
   *
   * @param key the boundary key
   * @param included is \p key can be equal to the key on the \p boundary_pos
   * @param snapshot_permutation the permuation use to search
   * @param[out] boundary_pos the leftmost position that satisfies the constraint
   * @param[out] is_end true if there are no more keys left to scan, otherwise false
   * @return int OB_SUCCESS on success, others on compare failure or KVQueue failure
   */
  int find_left_boundary_(
      BtreeKey key, bool included, Permutation snapshot_permutation, int &boundary_pos, bool &is_end);
  /**
   * @brief find the rightmost postion \p boundary_pos, such that \p key is less than
   * (or equal to, depending on the \p included parameter) the key at \p boundary_pos.
   *
   * @param key the boundary key
   * @param included is \p key can be equal to the key on the \p boundary_pos
   * @param snapshot_permutation the permuation use to search
   * @param[out] boundary_pos the rightmost position that satisfies the constraint
   * @param[out] is_end true if there are no more keys left to scan, otherwise false
   * @return int OB_SUCCESS on success, others on compare failure or KVQueue failure
   */
  int find_right_boundary_(
      BtreeKey key, bool included, Permutation snapshot_permutation, int &boundary_pos, bool &is_end);
};

template <typename BtreeKey, typename BtreeVal>
class InternalNode : public BtreeNode<BtreeKey, BtreeVal> {
private:
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  using BtreeKV = BtreeKV<BtreeKey, BtreeVal>;

public:
  OB_INLINE void set_leftmost_child(BtreeVal leftmost_child)
  {
    leftmost_child_ = leftmost_child;
  }
  virtual int split_and_insert(BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key) override;
  /**
   * @brief Search \p key on InternalNode.
   * @param key the search key.
   * @param[out] val the corresponding child node of the key.
   * @return OB_SUCCESS on success, others on compare fail.
   */
  virtual int search(const BtreeKey key, BtreeVal &val) override;
  virtual uint8_t get_level() const override
  {
    return level_;
  }
  OB_INLINE void set_level(int level)
  {
    level_ = level;
  }
  OB_INLINE BtreeVal get_leftmost_child() const
  {
    return leftmost_child_;
  }

private:
  BtreeVal leftmost_child_;
  uint8_t level_;
};

template <typename BtreeKey, typename BtreeVal>
class BtreeNodeAllocator {
private:
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  using LeafNode = LeafNode<BtreeKey, BtreeVal>;
  using InternalNode = InternalNode<BtreeKey, BtreeVal>;

public:
  BtreeNodeAllocator(common::ObIAllocator &allocator) : allocator_(allocator)
  {}
  virtual ~BtreeNodeAllocator()
  {}
  int make_leaf(LeafNode *&leaf)
  {
    int ret = OB_SUCCESS;
    void *block = allocator_.alloc(sizeof(LeafNode));
    if (OB_ISNULL(block)) {
      ret = OB_ALLOCATE_MEMORY_FAILED;
    } else {
      leaf = new (block) LeafNode;
    }
    return ret;
  }
  int make_internal(InternalNode *&internal)
  {
    int ret = OB_SUCCESS;
    void *block = allocator_.alloc(sizeof(InternalNode));
    if (OB_ISNULL(block)) {
      ret = OB_ALLOCATE_MEMORY_FAILED;
    } else {
      internal = new (block) InternalNode;
    }
    return ret;
  }
  OB_INLINE void free_node(BtreeNode *node)
  {
    allocator_.free(node);
  }

private:
  common::ObIAllocator &allocator_;
};

template <typename BtreeKey, typename BtreeVal>
class Path {
private:
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;

public:
  Path() : depth_(0)
  {}
  ~Path()
  {}
  void reset()
  {
    depth_ = 0;
  }
  OB_INLINE int push(BtreeNode *node, Version version)
  {
    int ret = OB_SUCCESS;
    if (depth_ >= MAX_DEPTH) {
      ret = OB_ARRAY_OUT_OF_RANGE;
    } else {
      path_[depth_].node_ = node;
      path_[depth_].version_ = version;
      depth_++;
    }
    return ret;
  }
  OB_INLINE int pop(BtreeNode *&node, Version &version)
  {
    int ret = OB_SUCCESS;
    if (OB_UNLIKELY(depth_ <= 0)) {
      ret = OB_ARRAY_OUT_OF_RANGE;
      node = nullptr;
    } else {
      depth_--;
      node = path_[depth_].node_;
      version = path_[depth_].version_;
    }
    return ret;
  }
  OB_INLINE int top(BtreeNode *&node, Version &version)
  {
    int ret = OB_SUCCESS;
    if (depth_ <= 0) {
      ret = OB_ARRAY_OUT_OF_RANGE;
      node = nullptr;
    } else {
      node = path_[depth_ - 1].node_;
      version = path_[depth_ - 1].version_;
    }
    return ret;
  }
  OB_INLINE bool empty() const
  {
    return 0 == depth_;
  }
  OB_INLINE void resize(int64_t depth)
  {
    depth_ = std::max(depth, 0l);
  }

private:
  enum { MAX_DEPTH = 16 };
  struct Item {
    Item() : node_(nullptr)
    {}
    ~Item()
    {}
    BtreeNode *node_;
    // int pos_;
    Version version_;
  };
  int64_t depth_;
  Item path_[MAX_DEPTH];
};

template <typename BtreeKey, typename BtreeVal>
class BtreeIterator;

template <typename BtreeKey, typename BtreeVal>
class ObKeyBtree {
private:
  using BtreeNodeAllocator = BtreeNodeAllocator<BtreeKey, BtreeVal>;
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  using LeafNode = LeafNode<BtreeKey, BtreeVal>;
  using InternalNode = InternalNode<BtreeKey, BtreeVal>;
  using Path = Path<BtreeKey, BtreeVal>;
  friend class BtreeIterator<BtreeKey, BtreeVal>;

public:
  ObKeyBtree(BtreeNodeAllocator &node_allocator) : node_allocator_(node_allocator), root_(nullptr)
  {}
  ~ObKeyBtree()
  {}
  int init()
  {
    int ret = OB_SUCCESS;
    LeafNode *leaf;
    if (OB_SUCC(node_allocator_.make_leaf(leaf))) {
      root_ = leaf;
    }
    return ret;
  }
  /**
   * @brief find the value corresponding to the \p key
   * @param key the key to be searched.
   * @param[out] val the corresponding value.
   * @return OB_SUCCESS on success, OB_ENTRY_NOT_EXIST if key is not exist, others on compare fail.
   */
  int search(BtreeKey key, BtreeVal &val);
  /**
   * @brief insert key-value pair into the Btree.
   * @param key the key to be inserted.
   * @param val the value to be inserted.
   * @return OB_SUCCESS on success, OB_ALLOCATE_MEMORY_FAILED on allocation fail, others on compare fail.
   */
  int insert(BtreeKey key, BtreeVal val);
  void dump(FILE *file)
  {
    root_->dump(file);
  }
  int set_key_range(BtreeIterator<BtreeKey, BtreeVal> &iter, const BtreeKey min_key, const bool include_min_key,
      const BtreeKey max_key, const bool include_max_key);
  int64_t size() const
  {
    return size_.value();
  }
  // TODO(shouluo): ensure no threads are reading the tree
  int destroy()
  {
    int ret = OB_SUCCESS;
    free_node(root_);
    return ret;
  }

private:
  int find_node(BtreeKey &key, uint8_t level, BtreeNode *&node, Version &version, Path &path);
  int split(BtreeNode *&node, BtreeKey key, BtreeVal value, Path &path);
  int make_new_root(BtreeKey fence_key, BtreeNode *node, BtreeNode *new_node)
  {
    int ret = OB_SUCCESS;
    InternalNode *new_root;
    if (OB_SUCC(node_allocator_.make_internal(new_root))) {
      new_root->set_level(node->get_level() + 1);
      new_root->set_leftmost_child(reinterpret_cast<BtreeVal>(node));
      new_root->insert(fence_key, reinterpret_cast<BtreeVal>(new_node));
      root_ = new_root;
    }
    return ret;
  }
  OB_INLINE BtreeNode *get_root() const
  {
    return root_;
  }
  void free_node(BtreeNode *node)
  {
    if (node->get_level() >= 1) {
      for (int i = 0; i < node->size(); i++) {
        BtreeNode *child = reinterpret_cast<BtreeNode *>(node->get_kv(i, node->get_permutation()).val_);
        free_node(child);
      }
    }
    node_allocator_.free_node(node);
  }
  BtreeNodeAllocator &node_allocator_;
  BtreeNode *root_;
  common::ObSimpleCounter size_;
};

template <typename BtreeKey, typename BtreeVal>
class BtreeIterator {
private:
  using LeafNode = LeafNode<BtreeKey, BtreeVal>;
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  using ObKeyBtree = ObKeyBtree<BtreeKey, BtreeVal>;
  using Path = Path<BtreeKey, BtreeVal>;
  using BtreeKV = BtreeKV<BtreeKey, BtreeVal>;
  using KVQueue = KVQueue<BtreeKey, BtreeVal>;

public:
  BtreeIterator()
  {}
  int init(ObKeyBtree &btree);
  int set_key_range(
      const BtreeKey min_key, const bool include_min_key, const BtreeKey max_key, const bool include_max_key);
  int scan_forward();
  int scan_backward();
  int get_next(BtreeKey &key, BtreeVal &val);
  void reset()
  {
    is_end_ = true;
    kv_queue_.reset();
  }
  bool is_reverse_scan() const
  {
    return is_backward_;
  }

private:
  int first_scan();
  ObKeyBtree *tree_;
  BtreeKey start_key_;
  BtreeKey end_key_;
  bool include_start_key_;
  bool include_end_key_;
  bool is_backward_;
  LeafNode *leaf_;
  bool is_end_;
  KVQueue kv_queue_;
};

};  // namespace keybtreeV2
};  // end namespace oceanbase

#include "ob_keybtreeV2.cpp"

#endif /* __OB_KEYBTREEV2_H__ */

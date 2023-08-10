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

enum { MAX_HEIGHT = 16 };

namespace keybtreeV2 {

template <typename BtreeKey, typename BtreeVal>
struct BtreeKV {
  BtreeKey key_;  // 8byte
  BtreeVal val_;  // 8byte
  TO_STRING_KV(K_(key), KP_(val));
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
  enum class Status { NOT_CHANGED = 0, INSERTED = 1, SPLITTED = 2 };

public:
  Version() : data_(0)
  {}
  Version(uint64_t data) : data_(data)
  {}
  Version get_stable_snapshot() const
  {
    Version stable_version;
    stable_version.data_ = ATOMIC_LOAD_ACQ(&this->data_);
    while (true) {
      for (int i = 0; i < MAX_TRY_COUNT && stable_version.data_ & DIRTY_MASK; i++) {
        PAUSE();
        stable_version.data_ = this->data_;
      }
      if (stable_version.data_ & DIRTY_MASK) {
        sched_yield();
      } else {
        break;
      }
    }
    // As a reader, we need to ensure that all our read operations on the node
    // occur after we take a snapshot of the version, so a memory fence is needed.
    // TODO(shouluo): Maybe an acquire fence is enough.
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
    return stable_version;
  }
  OB_INLINE Version get_unstable_snapshot() const
  {
    return Version(data_);
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
    // So a memory fence is needed here.
    // TODO(shouluo): Maybe an acquire fence is enough.
    bool bool_ret = false;
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
    if (is_splitting()) {
      bool_ret = true;
    } else if (get_vsplit() != snapshot_version.get_vsplit()) {
      bool_ret = true;
    }
    return bool_ret;
  }
  bool has_inserted(Version &snapshot_version) const
  {
    bool bool_ret = false;
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
    if (is_inserting()) {
      bool_ret = true;
    } else if (get_vinsert() != snapshot_version.get_vinsert()) {
      bool_ret = true;
    }
    return bool_ret;
  }
  Status has_changed(Version &snapshot_version) const
  {
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
    Status s_ret = Status::NOT_CHANGED;
    Version v = get_unstable_snapshot();
    if (((v.data_ ^ snapshot_version.data_) & (~LATCH_BIT)) != 0) {
      if (v.is_splitting() || v.get_vsplit() != snapshot_version.get_vsplit()) {
        s_ret = Status::SPLITTED;
      } else if (FALSE_IT(v = get_stable_snapshot())) {
        // empty
      } else if (v.get_vsplit() != snapshot_version.get_vsplit()) {
        s_ret = Status::SPLITTED;
      } else if (v.get_vinsert() != snapshot_version.get_vinsert()) {
        s_ret = Status::INSERTED;
        snapshot_version = v;
      } else {
        // empty
      }
    }
    return s_ret;
  }
  void set_inserting()
  {
    data_ = data_ | INSERTING_BIT;
    // As a writer, we need to make sure that all our write operations to the node happen after
    // we set the inserting/splitting bit, otherwise the reader may not find themselves
    // having read an inconsistent state. So a memory fence is needed.
    // TODO(shouluo): Maybe a release fence is enough.
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
  }
  void set_splitting()
  {
    data_ = data_ | SPLITTING_BIT;
    __atomic_thread_fence(__ATOMIC_SEQ_CST);
  }
  void latch()
  {
    uint64_t old_data = 0;
    uint64_t new_data = 0;
    while (true) {
      for (int i = 0; i < MAX_TRY_COUNT; i++) {
        uint64_t old_data = data_ & (~LATCH_BIT);
        uint64_t new_data = old_data | LATCH_BIT;
        if (ATOMIC_BCAS(&data_, old_data, new_data)) {
          return;
        } else {
          PAUSE();
        }
      }
      sched_yield();
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
    // Otherwise readers may not find themselves having read an inconsistent state.
    ATOMIC_BCAS(&data_, data_, data_ & ULATCH_MASK);
  }
  void reset()
  {
    data_ = 0;
  }
  TO_STRING_KV(K_(data));

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
  enum { MAX_TRY_COUNT = 100 };
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
 * of the Node, and the remaining are 15 4-bit integers that represent the actual  * This 64-bit Permutation can be
divided into 16 4-bit sub fields. The lowest 4 bits represent the size position of the corresponding
 * key in the correct order.
 *
 *  _______________________
 * |...| 2 | 1 | 0 | size |
 * |___|___|___|___|______|
 *       4   4   4   4bit   --- 64bit
 *
 * Note that the methods of Permutation DO NOT check the validity of the parameters(e.g. pos should be less than
 * max size)
 */
class Permutation {
public:
  Permutation() : data_(0)
  {}
  Permutation(uint64_t data) : data_(data)
  {}
  Permutation(Permutation &other) : data_(other.data_)
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
  OB_INLINE uint64_t get_raw() const
  {
    return data_;
  }
  TO_STRING_KV(K_(data));

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
  BtreeNode() : level_(0), version_(), permutation_()
  {}
  void reset()
  {
    version_.reset();
    permutation_.reset();
  }
  void dump(FILE *file);
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
  OB_INLINE int size() const
  {
    return permutation_.size();
  }
  OB_INLINE BtreeKV &get_kv(const int pos, const Permutation &snapshot_permutation)
  {
    int real_pos = snapshot_permutation.at(pos);
    return kvs_[real_pos];
  }
  OB_INLINE BtreeKV &get_kv(const int pos)
  {
    return kvs_[pos];
  }
  OB_INLINE uint8_t get_level() const
  {
    return level_;
  }
  OB_INLINE void set_level(uint8_t level)
  {
    level_ = level;
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
  DEFINE_VIRTUAL_TO_STRING({
    J_KV(K_(version), K_(permutation), K_(level));
    J_COMMA();
    J_NAME("kvs_");
    J_COLON();
    (void)databuff_print_obj_array(buf, buf_len, pos, kvs_, size());
  })

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
  uint8_t level_;
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
  OB_INLINE void push(const BtreeKV &data)
  {
    items_[idx(push_++)] = data;
  }
  OB_INLINE void pop(BtreeKV &data)
  {
    data = items_[idx(pop_++)];
  }
  OB_INLINE int64_t size() const
  {
    return push_ - pop_;
  }
  OB_INLINE bool empty() const
  {
    return push_ == pop_;
  }
  DEFINE_TO_STRING({
    J_KV(K_(push), K_(pop));
    J_COMMA();
    J_NAME("items_");
    J_COLON();
    (void)databuff_print_obj_array(buf, buf_len, pos, &items_[idx(pop_)], size());
  })
private:
  int64_t idx(const int64_t x) const
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
   * @brief Scan the leaf node, push the key-value pair in the interval [ \p start_key, \p end_key ]
            to \p kv_queue.
   *
   * @param start_key the start key of the scan
   * @param end_key the end key of the scan
   * @param exclude_start_key true if the start key should be included in the result
                              (if the key is in the node)
   * @param exclude_end_key true if the start key should be included in the scan result
                            (if the key is in the node)
   * @param is_backward true if it is a backward scan
   * @param[out] kv_queue the scan result queue
   * @param[out] is_end true if there are no more keys left to scan
   * @return int OB_SUCCESS on success, others on failure
   */
  int scan(BtreeKey start_key, BtreeKey end_key, bool exclude_start_key, bool exclude_end_key, bool is_backward,
      KVQueue &kv_queue, bool &is_end);
  /**
   * @brief Scan the leaf node, push the key-value pair in the interval [MIN_KEY, \p end_key]
            to \p kv_queue.
   *
   * @param end_key the end key of the scan
   * @param exclude_end_key true if the start key should be included in the scan result
                            (if the key is in the node)
   * @param is_backward true if it is a backward scan
   * @param[out] kv_queue the scan result queue
   * @param[out] is_end true if there are no more keys left to scan
   * @return int OB_SUCCESS on success, others on failure
   */
  int scan(BtreeKey end_key, bool exclude_end_key, bool is_backward, KVQueue &kv_queue, bool &is_end);
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
  INHERIT_TO_STRING_KV("BtreeNode", BtreeNode, KP_(prev), KP_(next));

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
  OB_INLINE BtreeVal get_leftmost_child() const
  {
    return leftmost_child_;
  }
  INHERIT_TO_STRING_KV("BtreeNode", BtreeNode, K_(leftmost_child));

private:
  BtreeVal leftmost_child_;
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
  ~BtreeNodeAllocator()
  {}
  int make_leaf(LeafNode *&leaf)
  {
    return make_node_<LeafNode>(leaf);
  }
  int make_internal(InternalNode *&internal)
  {
    return make_node_<InternalNode>(internal);
  }
  OB_INLINE void free_node(BtreeNode *node)
  {
    allocator_.free(node);
  }
  void reset()
  {
    // TODO(shouluo): Modify this if we change allocator
  }

private:
  template <typename NodeType>
  int make_node_(NodeType *&node)
  {
    int ret = OB_SUCCESS;
    void *block = allocator_.alloc(sizeof(NodeType));
    if (OB_ISNULL(block)) {
      ret = OB_ALLOCATE_MEMORY_FAILED;
    } else {
      node = new (block) NodeType;
    }
    return ret;
  }
  common::ObIAllocator &allocator_;
};

template <typename T1, typename T2>
class BtreeVector {
public:
  BtreeVector() : size_(0)
  {}
  ~BtreeVector()
  {}
  void reset()
  {
    size_ = 0;
  }
  OB_INLINE void push(T1 item1, T2 item2)
  {
    data_[size_].item1_ = item1;
    data_[size_].item2_ = item2;
    ++size_;
  }
  OB_INLINE void pop(T1 &item1, T2 &item2)
  {
    --size_;
    item1 = data_[size_].item1_;
    item2 = data_[size_].item2_;
  }
  OB_INLINE bool empty() const
  {
    return 0 == size_;
  }
  OB_INLINE int size() const
  {
    return size_;
  }
  OB_INLINE void resize(int size)
  {
    size_ = std::max(size, 0);
  }
  void get(int idx, T1 &item1, T2 &item2)
  {
    item1 = data_[idx].item1_;
    item2 = data_[idx].item2_;
  }
  DEFINE_TO_STRING({
    J_KV(K_(size));
    J_COMMA();
    J_NAME("data_");
    J_COLON();
    (void)databuff_print_obj_array(buf, buf_len, pos, data_, size_);
  });

private:
  struct Pair {
    T1 item1_;
    T2 item2_;
    TO_STRING_KV(K_(item1), K_(item2));
  };
  int size_;
  Pair data_[MAX_HEIGHT];
};

template <typename BtreeKey, typename BtreeVal>
class BtreeIterator;

class NodeLatchGuard {
public:
  NodeLatchGuard() : size_(0)
  {}
  ~NodeLatchGuard()
  {
    while (size_ > 0) {
      pop();
    }
  }
  void push(Version *v)
  {
    stack_[size_++] = v;
    v->latch();
  }
  void pop()
  {
    stack_[--size_]->unlatch();
  }

private:
  int size_;
  Version *stack_[MAX_HEIGHT];
};

// TODO(shouluo): InternalNode's template
template <typename BtreeKey, typename BtreeVal>
class ObKeyBtree {
private:
  using BtreeNodeAllocator = BtreeNodeAllocator<BtreeKey, BtreeVal>;
  using BtreeNode = BtreeNode<BtreeKey, BtreeVal>;
  using LeafNode = LeafNode<BtreeKey, BtreeVal>;
  using InternalNode = InternalNode<BtreeKey, BtreeVal>;
  using Path = BtreeVector<BtreeNode *, Version>;
  using NodePairArray = BtreeVector<BtreeNode *, InternalNode *>;
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
  int set_key_range(BtreeIterator<BtreeKey, BtreeVal> &iter, const BtreeKey min_key, const bool exclude_min_key,
      const BtreeKey max_key, const bool exclude_max_key);
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
  TO_STRING_KV(KP_(root));

private:
  int find_node(BtreeKey key, uint8_t level, BtreeNode *&node, Version &version, Path &path);
  int split(BtreeNode *&node, BtreeKey key, BtreeVal value, Path &path);
  int pre_alloc_nodes(BtreeKey key, Path &path, NodePairArray &stack);
  OB_INLINE BtreeNode *get_root() const
  {
    return root_;
  }
  void free_node(BtreeNode *node)
  {
    if (node->get_level() >= 1) {
      for (int i = 0; i < node->size(); i++) {
        BtreeNode *child = reinterpret_cast<BtreeNode *>(node->get_kv(i).val_);
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
  using Path = BtreeVector<BtreeNode *, Version>;
  using BtreeKV = BtreeKV<BtreeKey, BtreeVal>;
  using KVQueue = KVQueue<BtreeKey, BtreeVal>;

public:
  BtreeIterator()
  {}
  int init(ObKeyBtree &btree);
  int set_key_range(
      const BtreeKey min_key, const bool exclude_min_key, const BtreeKey max_key, const bool exclude_max_key);
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
  TO_STRING_KV(K_(tree), K_(start_key), K_(end_key), K_(exclude_start_key), K_(exclude_end_key), K_(is_backward),
      K_(leaf), K_(is_end), K_(kv_queue));

private:
  int first_scan();
  ObKeyBtree *tree_;
  BtreeKey start_key_;
  BtreeKey end_key_;
  bool exclude_start_key_;
  bool exclude_end_key_;
  bool is_backward_;
  LeafNode *leaf_;
  bool is_end_;
  KVQueue kv_queue_;
};

};  // namespace keybtreeV2
};  // end namespace oceanbase

#include "ob_keybtreeV2.cpp"

#endif /* __OB_KEYBTREEV2_H__ */

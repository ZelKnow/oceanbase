/**
 * Copyright (c) 2021 OceanBase
 * OceanBase CE is licensed under Mulan PubL v2.
 * You can use this software according to the terms and conditions of the Mulan
 * PubL v2. You may obtain a copy of Mulan PubL v2 at:
 *          http://license.coscl.org.cn/MulanPubL-2.0
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY
 * KIND, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
 * NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE. See the
 * Mulan PubL v2 for more details.
 */

#include "ob_keybtreeV2.h"

namespace oceanbase {
namespace keybtreeV2 {

template <typename BtreeKey, typename BtreeVal>
void BtreeNode<BtreeKey, BtreeVal>::dump(FILE *file)
{
  if (get_level() == 0) {
    using LeafNode = LeafNode<BtreeKey, BtreeVal>;
    LeafNode *leaf = reinterpret_cast<LeafNode *>(this);
    fprintf(file, "%p LeafNode: next: %p, prev: %p\n", this, leaf->get_next(), leaf->get_prev());
    for (int i = 0; i < size(); i++) {
      fprintf(file, "%s,", to_cstring(this->get_kv(i, permutation_).key_.get_ptr()));
    }
    fprintf(file, "\n------------------\n");
  } else {
    using InternalNode = InternalNode<BtreeKey, BtreeVal>;
    InternalNode *internal = reinterpret_cast<InternalNode *>(this);
    fprintf(file, "%p InternalNode\n", this);
    fprintf(file, "%p,", reinterpret_cast<BtreeNode *>(internal->get_leftmost_child()));
    for (int i = 0; i < size(); i++) {
      fprintf(file, "%s,", to_cstring(this->get_kv(i, permutation_).key_.get_ptr()));
    }
    fprintf(file, "\n------------------\n");
    reinterpret_cast<BtreeNode *>(internal->get_leftmost_child())->dump(file);
    for (int i = 0; i < size(); i++) {
      reinterpret_cast<BtreeNode *>(this->get_kv(i, permutation_).val_)->dump(file);
    }
  }
}

template <typename BtreeKey, typename BtreeVal>
int BtreeNode<BtreeKey, BtreeVal>::search_(const BtreeKey key, const Permutation &snapshot_permutation, int &pos)
{
  int ret = OB_SUCCESS;
  int n = snapshot_permutation.size();
  int l = 0;
  int r = n - 1;
  int cmp = 0;
  int mid = 0;
  BtreeKey mid_key;

  while (OB_SUCC(ret) && l < r) {
    mid = (l + r + 1) / 2;
    mid_key = get_kv(mid, snapshot_permutation).key_;
    if (OB_FAIL(mid_key.compare(key, cmp))) {
      TRANS_LOG(ERROR, "compare error", K(ret), K(mid_key), K(key));
    } else if (cmp <= 0) {
      l = mid;
    } else {
      r = mid - 1;
    }
  }

  if (OB_FAIL(ret)) {
    // empty
  } else if (OB_LIKELY(l != 0)) {
    // empty
  } else if (OB_UNLIKELY(n == 0)) {
    l = -1;
  } else if (OB_FAIL(get_kv(0, snapshot_permutation).key_.compare(key, cmp))) {
    TRANS_LOG(ERROR, "compare error", K(ret), K(get_kv(0, snapshot_permutation).key_), K(key));
  } else if (cmp > 0) {
    l = -1;
  }
  pos = l;

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
void BtreeNode<BtreeKey, BtreeVal>::copy_to_(BtreeNode *other, int start, int end)
{
  for (int i = start; i <= end; i++) {
    other->kvs_[other->size()] = kvs_[permutation_.at(i)];
    other->permutation_.insert(other->size(), other->size());
  }
}

template <typename BtreeKey, typename BtreeVal>
int BtreeNode<BtreeKey, BtreeVal>::insert(BtreeKey key, BtreeVal val)
{
  int pos = 0;
  int ret = OB_SUCCESS;

  if(OB_FAIL(search_(key, permutation_, pos))) {
    TRANS_LOG(ERROR, "Search failed.", K(ret), K(key));
  } else {
    kvs_[permutation_.size()] = BtreeKV{key, val};
    __atomic_thread_fence(__ATOMIC_RELEASE);
    permutation_.insert(pos + 1, permutation_.size());
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::search(const BtreeKey key, BtreeVal &val)
{
  int ret = OB_SUCCESS;
  // snapshot the permutation, prevent seeing inconsistent state due to
  // concurrent writers modifying permutation.
  Permutation snapshot_permutation(this->permutation_);
  int pos = -1;
  int cmp = 0;

  if (OB_FAIL(this->search_(key, snapshot_permutation, pos))) {
    TRANS_LOG(ERROR, "search_ failed.", K(ret), K(key));
  } else if (OB_LIKELY(pos >= 0)) {
    BtreeKV &kv = this->get_kv(pos, snapshot_permutation);
    if (OB_FAIL(kv.key_.compare(key, cmp))) {
      TRANS_LOG(ERROR, "Compare failed", K(ret), K(key));
    } else if (cmp != 0) {
      ret = OB_ENTRY_NOT_EXIST;
    } else {
      val = kv.val_;
    }
  } else {
    ret = OB_ENTRY_NOT_EXIST;
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::split_and_insert(BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key)
{
  int ret = OB_SUCCESS;
  int n = this->size();
  int cmp = 0;
  LeafNode *new_leaf = reinterpret_cast<LeafNode *>(new_node);

  // Redistribute key-values. First copy the right half of node to the new leaf.
  this->copy_to_(new_leaf, (n + 1) / 2, n - 1);
  // Then temporarily copy the left half to the new leaf
  this->copy_to_(new_leaf, 0, (n - 1) / 2);
  // Reset this node's permutation, let size == 0.
  this->permutation_.reset();
  // Copy back the right half (originally the left half of this node) from the new leaf.
  new_leaf->copy_to_(this, n / 2, n - 1);
  // Set the size of new leaf to n/2, retain only the left half (originally the right half of this node).
  new_leaf->permutation_.set_size(n / 2);

  fence_key = new_leaf->get_kv(0, new_leaf->permutation_).key_;

  // insert the key-value
  if (OB_FAIL(fence_key.compare(key, cmp))) {
    TRANS_LOG(ERROR, "Compare failed.", K(ret), K(fence_key), K(key));
  } else if (cmp > 0) {
    ret = this->insert(key, val);
  } else {
    ret = new_leaf->insert(key, val);
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::find_left_boundary_(
    BtreeKey key, bool excluded, Permutation snapshot_permutation, int &boundary_pos, bool &is_end)
{
  int ret = OB_SUCCESS;
  int cmp = 0;
  is_end = true;

  if (OB_FAIL(key.compare(this->get_kv(0, snapshot_permutation).key_, cmp))) {
    TRANS_LOG(ERROR, "compare failed.", K(ret), K(key), K(this->get_kv(0, snapshot_permutation).key_));
  } else if (cmp < 0) {
    // Key is less than ALL keys on the node.
    boundary_pos = 0;
    is_end = false;
  } else if (OB_UNLIKELY(cmp == 0)) {
    // Key is equal to the minimum keys on the node.
    boundary_pos = excluded ? 1 : 0;
  } else if (OB_FAIL(this->search_(key, snapshot_permutation, boundary_pos))) {
    TRANS_LOG(ERROR, "search_ failed", K(ret), K(key));
  } else if (!excluded && OB_FAIL(key.compare(this->get_kv(boundary_pos, snapshot_permutation).key_, cmp))) {
    TRANS_LOG(ERROR, "compare failed", K(ret), K(key), K(this->get_kv(boundary_pos, snapshot_permutation).key_));
  } else if (excluded || cmp != 0) {
    // case 1: excluded.
    // case 2: included, but cmp != 0. 
    // In both case, because search_ is lowerbound, we need to increase boundary_pos.
    boundary_pos++;
  } else {
    // case 3: included, and cmp == 0.
    // In this case, we stay at boundary_pos.
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::find_right_boundary_(
    BtreeKey key, bool excluded, Permutation snapshot_permutation, int &boundary_pos, bool &is_end)
{
  int ret = OB_SUCCESS;
  int cmp = 0;
  is_end = true;

  if (OB_FAIL(key.compare(this->get_kv(snapshot_permutation.size() - 1, snapshot_permutation).key_, cmp))) {
    TRANS_LOG(ERROR, "compare failed.", K(ret), K(key), K(this->get_kv(snapshot_permutation.size() - 1, snapshot_permutation).key_));
  } else if (cmp > 0) {
    // Key is greater than ALL keys on the node.
    is_end = false;
    boundary_pos = snapshot_permutation.size() - 1;
  } else if (OB_UNLIKELY(cmp == 0)) {
    // Key is equal to the maximum key on the node.
    boundary_pos = excluded ? snapshot_permutation.size() - 2 : snapshot_permutation.size() - 1;
  } else if (OB_FAIL(this->search_(key, snapshot_permutation, boundary_pos))) {
    TRANS_LOG(ERROR, "search_ failed.", K(ret), K(key));
  } else if (boundary_pos == -1) {
    // Key is less than ALL keys on the node.
  } else if (excluded && OB_FAIL(key.compare(this->get_kv(boundary_pos, snapshot_permutation).key_, cmp))) {
    TRANS_LOG(ERROR, "compare failed.", K(ret), K(key), K(this->get_kv(boundary_pos, snapshot_permutation).key_));
  } else if (excluded && cmp == 0) {
    // Case 1: excluded, and cmp == 0.
    // In this case, key is equal to the key at boundary_pos, so we don't need this key
    boundary_pos--;
  } else {
    // Case 2: included
    // Case 3: excluded, but cmp != 0.
    // In both case, we stay at boundary_pos.
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::scan_(
    int start_pos, int end_pos, bool is_backward, KVQueue &kv_queue, Permutation snapshot_permutation)
{
  int ret = OB_SUCCESS;
  kv_queue.reset();
  int direction = is_backward ? -1 : 1;
  for (int i = start_pos; (end_pos - i) * direction >= 0 && OB_SUCC(ret); i += direction) {
    if(OB_FAIL(kv_queue.push(this->get_kv(i, snapshot_permutation)))) {
      TRANS_LOG(ERROR, "push to kv queue failed.", K(ret), K(start_pos), K(end_pos));
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::scan(BtreeKey start_key, BtreeKey end_key, bool exclude_start_key,
    bool exclude_end_key, bool is_backward, KVQueue &kv_queue, bool &is_end)
{
  int ret = OB_SUCCESS;
  int start_pos = -1;
  int end_pos = -1;
  Permutation snapshot_permutation(this->permutation_);
  int cmp = 0;
  is_end = false;

  if (!is_backward) {
    if (OB_FAIL(find_left_boundary_(start_key, exclude_start_key, snapshot_permutation, start_pos, is_end))) {
      TRANS_LOG(ERROR, "find left boundary failed.", K(ret), K(start_key), K(end_key), K(exclude_start_key), K(exclude_end_key), K(is_backward));
    } else if (OB_FAIL(find_right_boundary_(end_key, exclude_end_key, snapshot_permutation, end_pos, is_end))) {
      TRANS_LOG(ERROR, "find right boundary failed.", K(ret), K(start_key), K(end_key), K(exclude_start_key), K(exclude_end_key), K(is_backward));
    } else {
      // empty
    }
  } else {
    if (OB_FAIL(find_right_boundary_(start_key, exclude_start_key, snapshot_permutation, start_pos, is_end))) {
      TRANS_LOG(ERROR, "find right boundary failed.", K(ret), K(start_key), K(end_key), K(exclude_start_key), K(exclude_end_key), K(is_backward));
    } else if (OB_FAIL(find_left_boundary_(end_key, exclude_end_key, snapshot_permutation, end_pos, is_end))) {
      TRANS_LOG(ERROR, "find left boundary failed.", K(ret), K(start_key), K(end_key), K(exclude_start_key), K(exclude_end_key), K(is_backward));
    } else {
      // empty
    }
  }

  if (OB_SUCC(ret)) {
    ret = scan_(start_pos, end_pos, is_backward, kv_queue, snapshot_permutation);
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::scan(
    BtreeKey end_key, bool exclude_end_key, bool is_backward, KVQueue &kv_queue, bool &is_end)
{
  int ret = OB_SUCCESS;
  int start_pos = -1;
  int end_pos = -1;
  Permutation snapshot_permutation(this->permutation_);
  int cmp = 0;
  is_end = false;

  if (!is_backward) {
    start_pos = 0;
    ret = find_right_boundary_(end_key, exclude_end_key, snapshot_permutation, end_pos, is_end);
  } else {
    start_pos = snapshot_permutation.size() - 1;
    ret = find_left_boundary_(end_key, exclude_end_key, snapshot_permutation, end_pos, is_end);
  }

  if (OB_SUCC(ret)) {
    ret = scan_(start_pos, end_pos, is_backward, kv_queue, snapshot_permutation);
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int InternalNode<BtreeKey, BtreeVal>::search(const BtreeKey key, BtreeVal &val)
{
  int ret = OB_SUCCESS;
  int pos = -1;
  Permutation snapshot_permutation(this->permutation_);

  if (OB_FAIL(this->search_(key, snapshot_permutation, pos))) {
    TRANS_LOG(ERROR, "search_ failed.", K(ret), K(key));
  } else if (pos == -1) {
    val = leftmost_child_;
  } else {
    val = this->get_kv(pos, snapshot_permutation).val_;
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int InternalNode<BtreeKey, BtreeVal>::split_and_insert(
    BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key)
{
  int ret = OB_SUCCESS;
  int n = this->size();
  int cmp = 0;
  InternalNode *new_internal = reinterpret_cast<InternalNode *>(new_node);

  // Redistribute key-values. First copy the right half(except the middle key) of node to the new internal node.
  this->copy_to_(new_internal, n / 2 + 1, n - 1);
  // Then temporarily copy the left half(except the middle key) to the new internal node.
  this->copy_to_(new_internal, 0, n / 2 - 1);
  // Take the middle key as fence key, push it to parent.
  fence_key = this->kvs_[this->permutation_.at(n / 2)].key_;
  // Let new internal node's leftmost child to be the fence_key's corresponding val.
  new_internal->leftmost_child_ = this->kvs_[this->permutation_.at(n / 2)].val_;
  // Reset this node's permutation, let size == 0.
  this->permutation_.reset();
  // Copy back the right half (originally the left half of this node) from the new internal node.
  new_internal->copy_to_(this, (n - 1) / 2, n - 2);
  // Set the size of new leaf to n/2, retain only the left half (originally the right half of this node).
  new_internal->permutation_.set_size((n - 1) / 2);

  // insert the key-value
  if (OB_FAIL(fence_key.compare(key, cmp))) {
    TRANS_LOG(ERROR, "compare failed.", K(ret), K(fence_key), K(key));
  } else if (cmp > 0) {
    this->insert(key, val);
  } else {
    new_internal->insert(key, val);
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::find_node(
    BtreeKey &key, uint8_t level, BtreeNode *&node, Version &version, Path &path)
{
  int ret = OB_SUCCESS;
  BtreeNode *child = nullptr;
  Version child_version;
  Version node_current_version;

  do {
    path.reset();
    // First, we need to get root and its stable version. And we need to double check the
    // root to make sure that the root has not changed during this period.
    do {
      node = get_root();
      version = node->get_version().get_stable_snapshot();
    } while (node != get_root());

    while (OB_SUCC(ret) && node->get_level() > level) {
      BtreeVal val;
      if (OB_FAIL(node->search(key, val))) {
        TRANS_LOG(ERROR, "node search failed.", K(ret), K(key));
      } else {
        child = reinterpret_cast<BtreeNode *>(val);
        child_version = child->get_version().get_stable_snapshot();
        node_current_version = node->get_version().get_stable_snapshot();
        if (node_current_version.has_splitted(version)) {
          // Node has splitted while we were searching for the child. In this case, the key we're looking for
          // may not be on the node anymore, so we need to retry from the root.
          break;
        } else if (node_current_version.has_inserted(version)) {
          // Node has inserted while we were searching for the child. In this case, the key we're looking for
          // is still on the node, so we can just retry from the current node.
          version = node_current_version;
          continue;
        }
        ret = path.push(node, version);
        node = child;
        version = child_version;
      }
    }
  } while (OB_SUCC(ret) && node->get_level() > level);
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::search(BtreeKey key, BtreeVal &val)
{
  int ret = OB_SUCCESS;
  BtreeNode *leaf = nullptr;
  Version version;
  Path path;
  bool is_found = false;

  while (OB_SUCC(ret) && !is_found) {
    if (OB_FAIL(find_node(key, 0, leaf, version, path))) {
      TRANS_LOG(ERROR, "Find leaf node failed.", K(ret), K(key));
    } else if (OB_FAIL(leaf->search(key, val))) {
      if (OB_ENTRY_EXIST != ret) {
        TRANS_LOG(ERROR, "Search key on leaf failed.", K(ret), K(key));
      }
    } else if (OB_UNLIKELY(leaf->get_version().has_splitted(version))) {
      // The leaf node has splitted while we were searching the key on the leaf, in this case,
      // the key we are looking for may not be on the leaf anymore, so we need to retry from root.
    } else {
      is_found = true;
    }
  }

  return ret;
}

// TODO(shouluo): Latch guard
template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::insert(BtreeKey key, BtreeVal val)
{
  int ret = OB_SUCCESS;
  BtreeNode *leaf;
  Version version;
  Path path;
  bool is_found = false;

  while (OB_SUCC(ret) && !is_found) {
    if (OB_FAIL(find_node(key, 0, leaf, version, path))) {
      TRANS_LOG(ERROR, "Find leaf node failed.", K(ret), K(key));
    } else {
      // Before modifying the leaf, we need to get the latch.
      leaf->get_version().latch();
      if (OB_UNLIKELY(leaf->get_version().has_splitted(version))) {
        // Leaf has splitted, so we need to find the leaf again.
        leaf->get_version().unlatch();
      } else {
        is_found = true;
      }
    }
  }

  if (OB_FAIL(ret)) {
    // empty
  } else if (OB_LIKELY(!leaf->is_full())) {  // leaf is not full, insert directly
    if (OB_FAIL(leaf->insert(key, val))) {
      TRANS_LOG(ERROR, "Insert key into leaf failed.", K(ret), K(key), K(val), KP(leaf));
    }
    leaf->get_version().unlatch();
  } else if (OB_FAIL(split(leaf, key, val, path))) {  // leaf is full, do split
    TRANS_LOG(ERROR, "Split the leaf failed.", K(ret), K(key), K(val), KP(leaf));
  } else {
    // empty
  }

  if (OB_SUCC(ret)) {
    size_.inc(1);
  }

  return ret;
}

// TODO(shouluo): format, rollback
template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::split(BtreeNode *&node, BtreeKey key, BtreeVal value, Path &path)
{
  int ret = OB_SUCCESS;
  int is_finished = false;
  BtreeNode *new_node = nullptr;
  LeafNode *leaf = reinterpret_cast<LeafNode *>(node);;
  LeafNode *new_leaf = nullptr;
  BtreeNode *parent = nullptr;
  InternalNode *new_internal = nullptr;
  BtreeKey fence_key;
  Version parent_version;

  leaf->get_version().set_splitting();

  if (OB_SUCC(node_allocator_.make_leaf(new_leaf))) {
    new_node = new_leaf;
    new_node->get_version().latch();
    new_node->get_version().set_splitting();

    // Adjust the sibling pointer of the leaf node
    new_leaf->set_prev(leaf);
    if (OB_NOT_NULL(leaf->get_next())) {
      leaf->get_next()->set_prev(new_leaf);
      new_leaf->set_next(leaf->get_next());
    }
    leaf->set_next(new_leaf);
  }

  if (OB_SUCC(ret) && OB_SUCC(node->split_and_insert(new_node, key, value, fence_key))) {
    // empty
  }

  while (OB_SUCC(ret) && !is_finished) {
    key = fence_key;

    if (OB_UNLIKELY(path.empty())) {
      // need to create new root
      if (OB_SUCC(make_new_root(key, node, new_node))) {
        node->get_version().unlatch();
        new_node->get_version().unlatch();
      }
      is_finished = true;
      break;
    }

    if (OB_SUCC(ret) && OB_SUCC(path.pop(parent, parent_version))) {
      parent->get_version().latch();
    }

    while (OB_SUCC(ret) && parent->get_version().has_splitted(parent_version)) {
      // at this point parent is no longer the parent of node, so we need to
      // find the real parent
      parent->get_version().unlatch();
      if (OB_SUCC(find_node(key, parent->get_level(), parent, parent_version, path))) {
        parent->get_version().latch();
      }
    }

    if (OB_FAIL(ret)) {
      // empty
    } else if (OB_LIKELY(!parent->is_full())) {
      // insert directly
      parent->get_version().set_inserting();
      ret = parent->insert(key, reinterpret_cast<BtreeVal>(new_node));
      node->get_version().unlatch();
      new_node->get_version().unlatch();
      parent->get_version().unlatch();
      is_finished = true;
    } else {
      // chained split
      parent->get_version().set_splitting();

      // At this point the above query can't come down so node can be safely unlatched
      node->get_version().unlatch();

      // prepare new node
      if (OB_SUCC(node_allocator_.make_internal(new_internal))) {
        new_internal->get_version().latch();
        new_internal->get_version().set_splitting();
        new_internal->set_level(parent->get_level());
      }

      if (OB_SUCC(ret) &&
          OB_SUCC(parent->split_and_insert(new_internal, key, reinterpret_cast<BtreeVal>(new_node), fence_key))) {
        // empty
      }
    }
    if (OB_SUCC(ret)) {
      // going up!
      new_node->get_version().unlatch();
      node = parent;
      new_node = new_internal;
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::set_key_range(BtreeIterator<BtreeKey, BtreeVal> &iter, const BtreeKey min_key,
    const bool exclude_min_key, const BtreeKey max_key, const bool exclude_max_key)
{
  int ret = OB_SUCCESS;
  if (OB_FAIL(iter.init(*this))) {
    // do nothing
  } else if (OB_FAIL(iter.set_key_range(min_key, exclude_min_key, max_key, exclude_max_key))) {
    // do nothing
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::init(ObKeyBtree &btree)
{
  int ret = OB_SUCCESS;
  tree_ = &btree;
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::set_key_range(
    const BtreeKey min_key, const bool exclude_min_key, const BtreeKey max_key, const bool exclude_max_key)
{
  int ret = OB_SUCCESS;
  int cmp = 0;
  if(OB_FAIL(max_key.compare(min_key, cmp))) {
    TRANS_LOG(ERROR, "compare failed.", K(ret), K(max_key), K(min_key));
  } else {
    is_backward_ = (cmp < 0);
    start_key_ = min_key;
    end_key_ = max_key;
    exclude_start_key_ = exclude_min_key;
    exclude_end_key_ = exclude_max_key;
    is_end_ = false;
    ret = first_scan();
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::first_scan()
{
  int ret = OB_SUCCESS;
  Path path;
  LeafNode *next_leaf = nullptr;
  BtreeNode *node = nullptr;
  Version version;
  bool is_done = false;

  while (OB_SUCC(ret) && !is_done) {
    kv_queue_.reset();
    if (OB_FAIL(tree_->find_node(start_key_, 0, node, version, path))) {
      // empty
    } else if (FALSE_IT(leaf_ = reinterpret_cast<LeafNode *>(node))) {

    } else if (OB_UNLIKELY(leaf_->size() == 0)) {
      is_end_ = true;
      is_done = true;
    } else if (OB_FAIL(leaf_->scan(
                   start_key_, end_key_, exclude_start_key_, exclude_end_key_, is_backward_, kv_queue_, is_end_))) {
      TRANS_LOG(ERROR, "scan failed.", K(ret), K(start_key_), K(end_key_), K(exclude_start_key_), K(exclude_end_key_), K(is_backward_));
    } else {
      next_leaf = leaf_->get_next();
      if (leaf_->get_version().has_splitted(version)) {
        // retry on leaf splitted
      } else {
        is_done = true;
      }
    }
  }

  if (OB_SUCC(ret)) {
    if (is_backward_) {
      if (OB_ISNULL(leaf_->get_prev())) {
        is_end_ = true;
      }
    } else {
      leaf_ = next_leaf;  // stay at next node
      if (OB_ISNULL(leaf_)) {
        is_end_ = true;
      }
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::scan_forward()
{
  int ret = OB_SUCCESS;
  Path path;
  LeafNode *next_leaf = nullptr;
  Version version;
  bool is_done = false;

  while (OB_SUCC(ret) && !is_done) {
    version = leaf_->get_version().get_stable_snapshot();
    if (OB_FAIL(leaf_->scan(end_key_, exclude_end_key_, false, kv_queue_, is_end_))) {
      TRANS_LOG(ERROR, "scan failed.", K(ret), K(end_key_), K(exclude_end_key_));
    } else {
      next_leaf = leaf_->get_next();
      Version leaf_new_version = leaf_->get_version().get_stable_snapshot();
      if (leaf_new_version.has_splitted(version)) {
        kv_queue_.reset();
        version = leaf_new_version;
      } else {
        leaf_ = next_leaf;  // move to next node
        if (OB_ISNULL(leaf_)) {
          is_end_ = true;
        }
        is_done = true;
      }
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::scan_backward()
{
  int ret = OB_SUCCESS;
  Path path;
  LeafNode *prev_leaf = nullptr;
  Version version;
  bool is_done = false;

  while (OB_SUCC(ret) && !is_end_ && !is_done) {
    kv_queue_.reset();
    prev_leaf = leaf_->get_prev();
    version = prev_leaf->get_version().get_stable_snapshot();
    if (prev_leaf != ATOMIC_LOAD_ACQ(&leaf_->get_prev())) {
      // empty
    } else if (OB_FAIL(prev_leaf->scan(end_key_, exclude_end_key_, true, kv_queue_, is_end_))) {
      TRANS_LOG(ERROR, "scan failed.", K(ret), K(end_key_), K(exclude_end_key_));
    } else if (prev_leaf->get_version().has_splitted(version)) {
      // Prev leaf has splitted, so we may be missing some keys. Retry.
    } else {
      leaf_ = prev_leaf;  // move to prev node
      if (OB_ISNULL(leaf_->get_prev())) {
        is_end_ = true;
      }
      is_done = true;
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::get_next(BtreeKey &key, BtreeVal &val)
{
  int ret = OB_SUCCESS;
  BtreeKV kv;

  if (kv_queue_.empty()) {
    if (OB_UNLIKELY(is_end_)) {
      ret = OB_ITER_END;
    } else if (is_backward_) {
      scan_backward();
    } else {
      scan_forward();
    }
    if (OB_UNLIKELY(kv_queue_.empty())) {
      ret = OB_ITER_END;
    }
  }

  if (OB_FAIL(ret)) {
    // empty
  } else if (OB_FAIL(kv_queue_.pop(kv))) {
    TRANS_LOG(ERROR, "Pop from queue failed.", K(ret));
  } else {
    key = kv.key_;
    val = kv.val_;
  }

  return ret;
}

}  // namespace keybtreeV2
}  // namespace oceanbase
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
  int cmp;
  int mid;
  while (l < r) {
    mid = (l + r + 1) / 2;
    if (OB_FAIL(get_kv(mid, snapshot_permutation).key_.compare(key, cmp))) {
      break;
    } else if (cmp <= 0) {
      l = mid;
    } else {
      r = mid - 1;
    }
  }
  if (OB_SUCC(ret) && OB_UNLIKELY(l == 0)) {
    if (OB_FAIL(get_kv(0, snapshot_permutation).key_.compare(key, cmp))) {
      // empty
    } else if (cmp > 0) {
      l = -1;
    }
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
  int pos;
  int ret = OB_SUCCESS;
  if (OB_UNLIKELY(permutation_.size() == 0)) {
    pos = -1;
  } else {
    ret = search_(key, permutation_, pos);
  }
  if (OB_SUCC(ret)) {
    kvs_[permutation_.size()] = BtreeKV{key, val};
    permutation_.insert(pos + 1, permutation_.size());
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::search(const BtreeKey key, BtreeVal &val)
{
  // snapshot the permutation, prevent seeing inconsistent state due to
  // concurrent writers modifying permutation.
  Permutation snapshot_permutation = this->permutation_;
  int ret = OB_SUCCESS;
  int pos;
  int cmp;

  if (OB_SUCC(this->search_(key, snapshot_permutation, pos)) && pos >= 0) {
    BtreeKV kv = this->get_kv(pos, snapshot_permutation);
    if (OB_SUCC(kv.key_.compare(key, cmp))) {
      if (cmp != 0) {
        ret = OB_ENTRY_NOT_EXIST;
      } else {
        val = kv.val_;
      }
    }
  } else if (OB_SUCC(ret) && pos == -1) {
    ret = OB_ENTRY_NOT_EXIST;
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::split_and_insert(BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key)
{
  int n = this->size();
  int cmp;
  int ret = OB_SUCCESS;

  LeafNode *new_leaf = reinterpret_cast<LeafNode *>(new_node);
  this->copy_to_(new_leaf, (n + 1) / 2, n - 1);
  this->copy_to_(new_leaf, 0, (n - 1) / 2);
  this->permutation_.reset();
  new_leaf->copy_to_(this, n / 2, n - 1);
  new_leaf->permutation_.set_size(n / 2);

  fence_key = new_leaf->get_kv(0, new_leaf->permutation_).key_;
  if (OB_SUCC(fence_key.compare(key, cmp))) {
    if (cmp > 0) {
      this->insert(key, val);
    } else {
      new_leaf->insert(key, val);
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::find_left_boundary_(
    BtreeKey key, bool included, Permutation snapshot_permutation, int &boundary_pos, bool &is_end)
{
  int ret = OB_SUCCESS;
  int cmp;
  is_end = true;

  if (OB_FAIL(key.compare(this->get_kv(0, snapshot_permutation).key_, cmp))) {
    // empty
  } else if (cmp < 0) {
    // key is less than ALL keys on the node
    boundary_pos = 0;
    is_end = false;
  } else if (cmp == 0) {
    // key is equal to the minimum keys on the node
    boundary_pos = included ? 0 : 1;
  } else if (OB_FAIL(this->search_(key, snapshot_permutation, boundary_pos))) {
    // empty
  } else if (included && OB_FAIL(key.compare(this->get_kv(boundary_pos, snapshot_permutation).key_, cmp))) {
    // empty
  } else if (!included || cmp != 0) {
    boundary_pos++;
  } else {
    // empty
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::find_right_boundary_(
    BtreeKey key, bool included, Permutation snapshot_permutation, int &boundary_pos, bool &is_end)
{
  int ret = OB_SUCCESS;
  int cmp;
  is_end = true;

  if (OB_FAIL(key.compare(this->get_kv(snapshot_permutation.size() - 1, snapshot_permutation).key_, cmp))) {

  } else if (cmp > 0) {
    // key is greater than ALL keys on the node
    is_end = false;
    boundary_pos = snapshot_permutation.size() - 1;
  } else if (cmp == 0) {
    // key is equal to the maximum key on the node
    boundary_pos = included ? snapshot_permutation.size() - 1 : snapshot_permutation.size() - 2;
  } else if (OB_FAIL(this->search_(key, snapshot_permutation, boundary_pos))) {
    // empty
  } else if (boundary_pos == -1) {
    // empty
  } else if (!included && OB_FAIL(key.compare(this->get_kv(boundary_pos, snapshot_permutation).key_, cmp))) {
    // empty
  } else if (!included && cmp == 0) {
    // key is equal to the key at boundary_pos, so we don't need this key
    boundary_pos--;
  } else {
    // empty
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
    ret = kv_queue.push(this->get_kv(i, snapshot_permutation));
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int LeafNode<BtreeKey, BtreeVal>::scan(BtreeKey start_key, BtreeKey end_key, bool include_start_key,
    bool include_end_key, bool is_backward, KVQueue &kv_queue, bool &is_end)
{
  int ret = OB_SUCCESS;
  int start_pos = -1;
  int end_pos = -1;
  Permutation snapshot_permutation = this->permutation_;
  int cmp = 0;
  is_end = false;

  if (!is_backward) {
    if (OB_FAIL(find_left_boundary_(start_key, include_start_key, snapshot_permutation, start_pos, is_end))) {
      // empty
    } else if (OB_FAIL(find_right_boundary_(end_key, include_end_key, snapshot_permutation, end_pos, is_end))) {
      // empty
    } else {
      // empty
    }
  } else {
    if (OB_FAIL(find_right_boundary_(start_key, include_start_key, snapshot_permutation, start_pos, is_end))) {
      // empty
    } else if (OB_FAIL(find_left_boundary_(end_key, include_end_key, snapshot_permutation, end_pos, is_end))) {
      // empty
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
    BtreeKey end_key, bool include_end_key, bool is_backward, KVQueue &kv_queue, bool &is_end)
{
  int ret = OB_SUCCESS;
  int start_pos;
  int end_pos;
  Permutation snapshot_permutation = this->permutation_;
  int cmp = 0;
  is_end = false;

  if (!is_backward) {
    start_pos = 0;
    ret = find_right_boundary_(end_key, include_end_key, snapshot_permutation, end_pos, is_end);
  } else {
    start_pos = snapshot_permutation.size() - 1;
    ret = find_left_boundary_(end_key, include_end_key, snapshot_permutation, end_pos, is_end);
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
  int pos;
  Permutation snapshot_permutation = this->permutation_;

  if (OB_SUCC(this->search_(key, snapshot_permutation, pos))) {
    if (pos == -1) {
      val = leftmost_child_;
    } else {
      val = this->get_kv(pos, snapshot_permutation).val_;
    }
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int InternalNode<BtreeKey, BtreeVal>::split_and_insert(
    BtreeNode *new_node, BtreeKey key, BtreeVal val, BtreeKey &fence_key)
{
  int n = this->size();
  int ret = OB_SUCCESS;
  int cmp;
  InternalNode *new_internal = reinterpret_cast<InternalNode *>(new_node);

  this->copy_to_(new_internal, n / 2 + 1, n - 1);
  this->copy_to_(new_internal, 0, n / 2 - 1);
  fence_key = this->kvs_[this->permutation_.at(n / 2)].key_;
  new_internal->leftmost_child_ = this->kvs_[this->permutation_.at(n / 2)].val_;
  this->permutation_.reset();
  new_internal->copy_to_(this, (n - 1) / 2, n - 2);
  new_internal->permutation_.set_size((n - 1) / 2);

  if (OB_SUCC(fence_key.compare(key, cmp))) {
    if (cmp > 0) {
      this->insert(key, val);
    } else {
      new_internal->insert(key, val);
    }
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::find_node(
    BtreeKey &key, uint8_t level, BtreeNode *&node, Version &version, Path &path)
{
  int ret = OB_SUCCESS;
  do {
    do {
      node = get_root();
      version = node->get_version().get_stable_snapshot();
    } while (node != get_root());  // get root and its version

    while (OB_SUCC(ret) && node->get_level() > level) {
      BtreeVal val;
      if (OB_SUCC(node->search(key, val))) {
        BtreeNode *child = reinterpret_cast<BtreeNode *>(val);
        Version child_version = child->get_version().get_stable_snapshot();
        Version node_current_version = node->get_version().get_stable_snapshot();
        if (node_current_version.has_splitted(version)) {
          break;  // retry from root
        } else if (node_current_version.has_inserted(version)) {
          version = node_current_version;
          continue;  // retry from node
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
  BtreeNode *leaf;
  Version version;
  Path path;
  bool is_found = false;

  while (OB_SUCC(ret) && !is_found) {
    if (OB_FAIL(find_node(key, 0, leaf, version, path))) {
      // empty
    } else if (OB_FAIL(leaf->search(key, val))) {
      // empty
    } else if (leaf->get_version().has_splitted(version)) {
      // the leaf node has splitted, retry
    } else {
      is_found = true;
    }
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::insert(BtreeKey key, BtreeVal val)
{
  int ret = OB_SUCCESS;
  BtreeNode *leaf;
  Version version;
  Path path;
  bool is_found = false;

  while (OB_SUCC(ret) && !is_found) {
    if (OB_SUCC(find_node(key, 0, leaf, version, path))) {
      leaf->get_version().latch();
      if (leaf->get_version().has_splitted(version)) {
        leaf->get_version().unlatch();
      } else {
        is_found = true;
      }
    }
  }

  if (OB_FAIL(ret)) {
    // empty
  } else if (!leaf->is_full()) {  // leaf is not full, insert directly
    if (OB_SUCC(leaf->insert(key, val))) {
      leaf->get_version().unlatch();
    }
  } else if (OB_SUCC(split(leaf, key, val, path))) {  // leaf is full, do split
    // empty
  } else {
    // empty
  }

  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::split(BtreeNode *&node, BtreeKey key, BtreeVal value, Path &path)
{
  // TODO(shouluo): format, rollback
  int ret = OB_SUCCESS;
  int is_finished = false;
  BtreeNode *new_node;
  LeafNode *leaf;
  LeafNode *new_leaf;
  BtreeNode *parent;
  InternalNode *new_internal;
  BtreeKey fence_key;
  Version parent_version;

  leaf = reinterpret_cast<LeafNode *>(node);
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

    if (path.empty()) {
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

    if (ret != OB_SUCCESS) {
      // empty
    } else if (!parent->is_full()) {
      // insert directly
      parent->get_version().set_inserting();
      parent->insert(key, reinterpret_cast<BtreeVal>(new_node));
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

template<typename BtreeKey, typename BtreeVal>
int ObKeyBtree<BtreeKey, BtreeVal>::set_key_range(BtreeIterator<BtreeKey, BtreeVal> &iter,
    const BtreeKey min_key,
    const bool include_min_key,
    const BtreeKey max_key,
    const bool include_max_key)
{
  int ret = OB_SUCCESS;
  if (OB_FAIL(iter.init(*this))) {
    // do nothing
  } else if (OB_FAIL(iter.set_key_range(min_key, include_min_key, max_key, include_max_key))) {
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

template<typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::set_key_range(const BtreeKey min_key,
    const bool include_min_key,
    const BtreeKey max_key,
    const bool include_max_key) {
  int ret = OB_SUCCESS;
  int cmp = 0;
  ret = max_key.compare(min_key, cmp);
  is_backward_ = (cmp < 0);
  start_key_ = min_key;
  end_key_ = max_key;
  include_start_key_ = include_min_key;
  include_end_key_ = include_max_key;
  is_end_ = false;
  if(OB_SUCC(ret)) {
    ret = first_scan();
  }
  return ret;
}

template<typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::first_scan() {
  int ret = OB_SUCCESS;
  Path path;
  LeafNode *next_leaf = nullptr;
  BtreeNode *node = reinterpret_cast<BtreeNode *>(leaf_);
  Version version;

  while (OB_SUCC(ret)) {
    if (OB_FAIL(tree_->find_node(start_key_, 0, node, version, path))) {
      // empty
    } else {
      leaf_ = reinterpret_cast<LeafNode *>(node);
    }
    if (OB_FAIL(ret)) {
      // empty
    } else if (OB_FAIL(leaf_->scan(
                   start_key_, end_key_, include_start_key_, include_end_key_, is_backward_, kv_queue_, is_end_))) {
      // empty
    } else {
      next_leaf = leaf_->get_next();
      if (leaf_->get_version().has_splitted(version)) {
        // retry on leaf splitted
        kv_queue_.reset();
      } else {
        break;
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
  LeafNode *next_leaf;
  Version version;

  while (OB_SUCC(ret)) {
    version = leaf_->get_version().get_stable_snapshot();
    if (OB_FAIL(leaf_->scan(end_key_, include_end_key_, false, kv_queue_, is_end_))) {
      // empty
    } else {
      Version leaf_new_version = leaf_->get_version().get_stable_snapshot();
      next_leaf = leaf_->get_next();
      if (leaf_new_version.has_splitted(version)) {
        kv_queue_.reset();
        version = leaf_new_version;
      } else {
        leaf_ = next_leaf;  // move to next node
        if (OB_ISNULL(leaf_)) {
          is_end_ = true;
        }
        break;
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
  LeafNode *prev_leaf;
  Version version;

  while (OB_SUCC(ret) && !is_end_) {
    prev_leaf = leaf_->get_prev();
    version = prev_leaf->get_version().get_stable_snapshot();
    if (prev_leaf != ATOMIC_LOAD_ACQ(&leaf_->get_prev())) {
      // empty
    } else if (OB_FAIL(prev_leaf->scan(end_key_, include_end_key_, true, kv_queue_, is_end_))) {
      // empty
    } else if (prev_leaf->get_version().has_splitted(version)) {
      kv_queue_.reset();
    } else {
      leaf_ = prev_leaf;  // move to prev node
      if (OB_ISNULL(leaf_->get_prev())) {
        is_end_ = true;
      }
      break;
    }
  }
  return ret;
}

template <typename BtreeKey, typename BtreeVal>
int BtreeIterator<BtreeKey, BtreeVal>::get_next(BtreeKey &key, BtreeVal &val)
{
  int ret = OB_SUCCESS;
  if (kv_queue_.empty()) {
    if (is_end_) {
      ret = OB_ITER_END;
    } else if (is_backward_) {
      scan_backward();
    } else {
      scan_forward();
    }
    if (kv_queue_.empty()) {
      ret = OB_ITER_END;
    }
  }

  if (OB_SUCC(ret)) {
    BtreeKV kv;
    kv_queue_.pop(kv);
    key = kv.key_;
    val = kv.val_;
  }

  return ret;
}

}  // namespace keybtreeV2
}  // namespace oceanbase
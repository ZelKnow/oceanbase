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

#ifndef __OB_MASSTREE_H__
#define __OB_MASSTREE_H__

#include "lib/allocator/ob_retire_station.h"
#include "masstree-beta/masstree.hh"
#include "masstree-beta/masstree_insert.hh"
#include "masstree-beta/masstree_get.hh"
#include "masstree-beta/masstree_split.hh"
#include "masstree-beta/mtIndexAPI.hh"

#define BTREE_ORDER 15

namespace oceanbase
{
namespace keybtree
{
using namespace oceanbase::common;

namespace masstree
{
    


template<typename KeyType, typename ValueType>
class MassTreeIndex
{
 public:

  typedef mt_index<Masstree::default_table> MapType;

  ~MassTreeIndex() {
    delete idx;
  }

  inline void swap_endian(uint64_t &i) {
    // Note that masstree internally treat input as big-endian
    // integer values, so we need to swap here
    // This should be just one instruction
    i = __bswap_64(i);
  }

  void UpdateThreadLocal(size_t thread_num) {}
  void AssignGCID(size_t thread_id) {}
  void UnregisterThread(size_t thread_id) {}

  bool insert(KeyType key, ValueType value, threadinfo *ti) {
    swap_endian(key);
    idx->put((const char*)&key, sizeof(KeyType), (const char*)&value, 8, ti);

    return true;
  }

  int insert(KeyType key, ValueType &value) {
    if (insert(key, value, idx->get_thread_info())) {
        return OB_SUCCESS;
    } else {
        return OB_ERROR;
    }
  }


  uint64_t find(KeyType key, std::vector<uint64_t> *v, threadinfo *ti) {
    Str val;
    swap_endian(key);
    idx->get((const char*)&key, sizeof(KeyType), val, ti);

    v->clear();
    if (val.s)
      v->push_back(*(uint64_t *)val.s);

    return 0;
  }

  bool upsert(KeyType key, uint64_t value, threadinfo *ti) {
    swap_endian(key);
    idx->put((const char*)&key, sizeof(KeyType), (const char*)&value, 8, ti);
    return true;
  }

  // uint64_t scan(KeyType key, int range, threadinfo *ti) {
  //   Str val;

  //   swap_endian(key);
  //   int key_len = sizeof(KeyType);

  //   for (int i = 0; i < range; i++) {
  //     idx->dynamic_get_next(val, (char *)&key, &key_len, ti);
  //   }

  //   return 0UL;
  // }

  uint64_t scan(KeyType key, int range, threadinfo *ti) {
    Str results[range];

    swap_endian(key);
    int key_len = sizeof(KeyType);

    int resultCount = idx->get_next_n(results, (char *)&key, &key_len, range, ti);
    //printf("scan: requested: %d, actual: %d\n", range, resultCount);
    return resultCount;
  }

  int64_t getMemory() const {
    return 0;
  }

  MassTreeIndex(uint64_t kt) {
    idx = new MapType{};

    threadinfo *main_ti = threadinfo::make(threadinfo::TI_MAIN, -1);
    idx->setup(main_ti);

    return;
  }

  MapType *idx;
};

}
}
}

#endif
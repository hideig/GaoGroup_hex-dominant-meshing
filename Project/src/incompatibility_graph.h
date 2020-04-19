/*
 * HXT - Copyright (C) <2016-2018> <Université catholique de Louvain (UCL), Belgique>
 *
 * List of the contributors to the development of HXT: see AUTHORS file.
 * Description and complete License: see LICENSE file.
 * 	
 * This program (HXT) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU Lesser
 * General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING and COPYING.LESSER files).  If not, 
 * see <http://www.gnu.org/licenses/>.
 */

#ifndef HXT_COMBINE_INCOMPATIBILITY_GRAPH_H_
#define HXT_COMBINE_INCOMPATIBILITY_GRAPH_H_

#include <cstddef>
#include <utility>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>

#include "hxt_tools.h"
#include "basic_types.h"

/**
* \file incompatibility_graph.h Implementation of the API function generating the incompatibility graph
* of an array of HXTCombineCell
*
* \author Kilian Verhetsel
*/

/**
* \addtogroup Combine
* @{
*/

/**
 * 兼容性检查时的多线程安全问题，生产者消费者模型
 * The difficulty of generating the incompatibility graph in parallel comes from
 * the fact that in a naive approach, multiple threads can try to access the
 * same row of the adjacency list when adding new edges.
 *
 * To avoid this, two types of threads are used:
 *   1. producers, which find incompatible pairs and
 *   2. consumers, which add pairs to the adjacency list.
 *
 * The two types of threads share pairs using a concurrent queue. A pair (h0,
 * h1) is always handled by the same consumer determined as a function of h0 so
 * that no two threads access the same row of the adjacency list at the same
 * time.
 */

/**
 * Amount of the pairs in an individual block. The synchronization costs of the
 * data structure used to communicate between threads are amortized over all
 * elements of this value.
 *
 * Keeping this value a multiple of 512 makes each block a multiple of the page
 * size on Intel CPUs (assuming 32-bit unsigned integers).
 */
static constexpr std::size_t BLOCK_SIZE = 4096;

class IncompatibilityBlock {
  HexIndex *data_;
  std::size_t size_;

public:
  IncompatibilityBlock(): data_(nullptr), size_(0) {}
  IncompatibilityBlock(HexIndex *data): data_(data), size_(0) {}

  std::size_t size() const { return size_; }

  std::pair<HexIndex, HexIndex> operator[](std::size_t i) const {
    return std::make_pair(data_[2*i], data_[2*i+1]);
  }

  void clear() { size_ = 0; }

  /**
   * Adds a pair to the block.
   *
   * Returns true if the block is now full.
   */
  bool push_back(HexIndex h0, HexIndex h1) {
    data_[2*size_ + 0] = h0;
    data_[2*size_ + 1] = h1;
    size_++;
    return size_ == BLOCK_SIZE;
  }
};

/**
 * Data structure used to share blocks of incompatible pairs from producers to
 * consumers.
 *
 * Every time elements are added or removed from this structure, they are
 * swapped with other valid buffers that can be shared between threads.
 */
class ConcurrentQueue {
  std::vector<IncompatibilityBlock> buffer_;
  std::size_t i_, size_;

  bool closed_;

  std::mutex mutex_;
  std::condition_variable emptyVar_;
  std::condition_variable fullVar_;

public:
  ConcurrentQueue(std::size_t bufferSize, IncompatibilityBlock *blocks):
    buffer_(bufferSize),
    i_(0), size_(0),
    closed_(false)
    {
      for (std::size_t i = 0; i < bufferSize; i++) buffer_[i] = blocks[i];
    }

  void push(IncompatibilityBlock &block) {
    {
      std::unique_lock<std::mutex> lock(mutex_);

      while (size_ == buffer_.size()) fullVar_.wait(lock);

      std::swap(buffer_[(i_ + size_) % buffer_.size()], block);
      size_++;
    }

    emptyVar_.notify_one();
  }

  void close() {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      closed_ = true;
    }

    emptyVar_.notify_all();
  }

  bool pop(IncompatibilityBlock &block) {
    std::unique_lock<std::mutex> lock(mutex_);
    while (size_ == 0 && !closed_)
      emptyVar_.wait(lock);

    if (size_ > 0) {
      std::swap(block, buffer_[i_]);
      i_++;
      if (i_ == buffer_.size()) i_ = 0;
      size_--;

      lock.unlock();
      fullVar_.notify_one();

      return true;
    }
    else {
      return false;
    }
  }
};

class IncompatibilityAdder {
  std::vector<std::unique_ptr<ConcurrentQueue>> &queues_;
  std::vector<IncompatibilityBlock> blocks_;

public:
  IncompatibilityAdder(std::vector<std::unique_ptr<ConcurrentQueue>> &queues,
                       IncompatibilityBlock *blocks):
    queues_(queues), blocks_(queues.size())
    {
      for (std::size_t i = 0; i < queues.size(); i++)
        blocks_[i] = blocks[i];
    }

  void flush() {
    for (std::size_t i = 0; i < queues_.size(); i++) {
      if (blocks_[i].size() != 0) queues_[i]->push(blocks_[i]);
    }
  }

  void addPair(HexIndex h0, HexIndex h1) {
    std::size_t i = queueIndex(h0);
    IncompatibilityBlock &block = blocks_[i];
    if (block.push_back(h0, h1)) {
      queues_[i]->push(block);
      block.clear();
    }
  }

private:
  std::size_t queueIndex(CellIndex i) {
    return i % queues_.size();
  }
};

class AdjacencyListUpdater {
  IncompatibilityBlock block_;
  ConcurrentQueue &queue_;

  HexIndex **adjacencyList_;
  HexIndex *degrees_;
  HexIndex *capacities_;

public:
  AdjacencyListUpdater(IncompatibilityBlock block, ConcurrentQueue &queue,
                       HexIndex **adjacencyList,
                       HexIndex *degrees, HexIndex *capacities):
    block_(std::move(block)), queue_(queue),
    adjacencyList_(adjacencyList), degrees_(degrees), capacities_(capacities)
    {}

  void operator()() {
    IncompatibilityBlock block(std::move(block_));
    int received = 0;

    while (queue_.pop(block)) {
      for (std::size_t i = 0; i < block.size(); i++) {
        std::pair<HexIndex, HexIndex> pair = block[i];
        addEntry(pair.first, pair.second);
        received++;
      }
    }
  }

private:
  void addEntry(HXTIndex rowId, HXTIndex value) {
    increaseCapacityIfNeeded(rowId);
    adjacencyList_[rowId][degrees_[rowId]++] = value;
  }

  void increaseCapacityIfNeeded(HXTIndex rowId) {
    HXTIndex oldSize = degrees_[rowId];
    HXTIndex oldCapa = capacities_[rowId];
    if (oldSize < oldCapa) return;

    HXTIndex newCapa = 3*oldCapa / 2;
    if (hxtRealloc((void**)&adjacencyList_[rowId], newCapa*sizeof(HXTIndex)) < 0) {
      throw std::bad_alloc();
    }
    capacities_[rowId] = newCapa;
  }
};

/**
 * Runs a function in parallel in order to add edges to a graph.
 *
 * The given function is called with three parameters:
 *
 * 1. An reference to an \ref IncompatibilityAdder. \ref
 *    IncompatibilityAdder::addPair should be used to add entries to the graph.
 * 2. The integer id of the thread.
 * 3. The total number of threads for which the function was called. This is not
 *    omp_get_num_threads(), since some threads are used to move entries from
 *    intermediate buffers to the shared graph.
 *
 * \param incompatibilities Adjacency list of the graph
 * \param degrees Size of each row of the adjacency list
 * \param capacity Allocated size of each row of the adjacency list, in number
 *   of entries.
 * \param f Function executed in parallel to add entries to the graph.
 */
template<typename Function>
void addToGraph(HXTIndex **incompatibilities, HXTIndex *degrees, HXTIndex *capacities,
                Function f) {
  std::size_t nthreads = omp_get_max_threads();

  std::size_t numProducers = (nthreads * 7) / 8;
  std::size_t numConsumers = nthreads - numProducers;
  if (numProducers == 0) numProducers++;
  if (numConsumers == 0) numConsumers++;

  std::size_t bufferSize = 8;
  std::size_t numBlocks = (numProducers + bufferSize + 1) * numConsumers;
  std::vector<HexIndex> blockData(numBlocks * BLOCK_SIZE * 2);

  std::vector<IncompatibilityBlock> blocks(numBlocks);
  for (std::size_t i = 0; i < numBlocks; i++) {
    blocks[i] = IncompatibilityBlock(&blockData[2*BLOCK_SIZE*i]);
  }

  std::vector<std::unique_ptr<ConcurrentQueue>> queues(numConsumers);
  for (std::size_t i = 0; i < numConsumers; i++) {
    queues[i] = std::unique_ptr<ConcurrentQueue>(
      new ConcurrentQueue(bufferSize,
                          &blocks[(numProducers+1)*numConsumers + i*bufferSize]));
  }

  std::atomic<std::size_t> numFinished(0);

  #pragma omp parallel
  {
    std::size_t id = omp_get_thread_num();
    if (id < numProducers) {
      IncompatibilityBlock *threadBlocks = blocks.data() + id*numConsumers;
      IncompatibilityAdder adder(queues, threadBlocks);
      f(std::ref(adder), id, numProducers);
      adder.flush();

      if (numFinished.fetch_add(1, std::memory_order_acq_rel)+1 == numProducers) {
        for (auto &queue : queues) queue->close();
      }
    }
    else {
      id -= numProducers;
      AdjacencyListUpdater fn(std::move(blocks[numProducers*numConsumers + id]),
                              *queues[id], incompatibilities, degrees, capacities);
      fn();
    }
  }
}

/**
* @}
*/


#endif

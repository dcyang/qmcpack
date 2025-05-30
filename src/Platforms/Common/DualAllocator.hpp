//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// Refactored from: OMPallocator.hpp
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file
 */
#ifndef QMCPLUSPLUS_DUAL_ALLOCATOR_H
#define QMCPLUSPLUS_DUAL_ALLOCATOR_H

#include <memory>
#include <atomic>
#include <exception>
#include "allocator_traits.hpp"

namespace qmcplusplus
{
extern std::atomic<size_t> dual_device_mem_allocated;

inline size_t getDualDeviceMemAllocated() { return dual_device_mem_allocated; }

/** Generalizes the DualMemorySpace allocator
 *  This provides a limited alternative to OMPallocator for testing/benchmarking
 *  without dependence of OMPTarget/ offload.
 *  It does not provide an alternative to OMPtarget transfer semantics so many production
 *  objects will not be functional if it is used as the allocator for the data objects they depend
 *  on.
 *  If you use DualAllocator at this time you need to handle data transfer yourself.
 *
 *  \todo the OMPTarget allocation can be a "device" allocator comparable to CUDAAllocator
 *  Then OMPallocator can be replaced by a DualAllocator<T, OffloadAllocator<T>, PinnedAllocator<T>>
 */
template<typename T, class DeviceAllocator, class HostAllocator = std::allocator<T>>
struct DualAllocator : public HostAllocator
{
  using Value        = typename HostAllocator::value_type;
  using Size         = typename HostAllocator::size_type;
  using Pointer      = typename HostAllocator::pointer;
  using ConstPointer = typename HostAllocator::const_pointer;

  DualAllocator() : device_ptr_(nullptr) {};
  DualAllocator(const DualAllocator&) : device_ptr_(nullptr) {}
  DualAllocator& operator=(const DualAllocator&) { device_ptr_ = nullptr; }
  template<class U, class V>
  DualAllocator(const DualAllocator<U, V>&) : device_ptr_(nullptr)
  {}

  template<class U>
  struct rebind
  {
    using other = DualAllocator<U,
                                typename std::allocator_traits<DeviceAllocator>::template rebind_alloc<U>,
                                typename std::allocator_traits<HostAllocator>::template rebind_alloc<U>>;
  };

  Value* allocate(std::size_t n)
  {
    static_assert(std::is_same<T, Value>::value, "DualAllocator and HostAllocator data types must agree!");
    if (device_ptr_ != nullptr)
      throw std::runtime_error("DualAllocator does not support device reallocation");
    Value* host_ptr = std::allocator_traits<HostAllocator>::allocate(allocator_, n);
    device_ptr_     = std::allocator_traits<DeviceAllocator>::allocate(device_allocator_, n);
    dual_device_mem_allocated += n * sizeof(T);
    return host_ptr;
  }

  void deallocate(Value* pt, std::size_t n)
  {
    dual_device_mem_allocated -= n * sizeof(T);
    std::allocator_traits<DeviceAllocator>::deallocate(device_allocator_, device_ptr_, n);
    std::allocator_traits<HostAllocator>::deallocate(allocator_, pt, n);
    device_ptr_ = nullptr;
  }

  void attachReference(const DualAllocator& from, std::ptrdiff_t ptr_offset)
  {
    device_ptr_ = const_cast<Pointer>(from.get_device_ptr()) + ptr_offset;
  }

  T* get_device_ptr() { return device_ptr_; }
  const T* get_device_ptr() const { return device_ptr_; }

  DeviceAllocator& get_device_allocator() { return device_allocator_; }
  const DeviceAllocator& get_device_allocator() const { return device_allocator_; }

private:
  HostAllocator allocator_;
  DeviceAllocator device_allocator_;
  T* device_ptr_;
};

template<typename T, class DeviceAllocator, class HostAllocator>
struct qmc_allocator_traits<DualAllocator<T, DeviceAllocator, HostAllocator>>
{
  using DualAlloc                          = DualAllocator<T, DeviceAllocator, HostAllocator>;
  static constexpr bool is_host_accessible = true;
  static constexpr bool is_dual_space      = true;

  static void fill_n(T* ptr, size_t n, const T& value) { qmc_allocator_traits<HostAllocator>::fill_n(ptr, n, value); }

  static void attachReference(const DualAlloc& from, DualAlloc& to, std::ptrdiff_t ptr_offset)
  {
    to.attachReference(from, ptr_offset);
  }

  /** update to the device, assumes you are copying starting with the implicit host_ptr.
   *
   *  These follow the openmp target semantics where you only provide the host
   *  side of a host_ptr device_ptr pair but the verb relates to what happens on the device.
   *
   *  This is primarily for testing to reduce ifdef code and single "flavor" testing
   *
   *  This a generic API and unlikely to be the best way to handle performance critical transfers,
   *  but if you have to use it or ifdef at a level above a xxxCUDA.cu or xxxOMPTarget.hpp file
   *  thats an issue. 
   */
  static void updateTo(DualAlloc& alloc, T* host_ptr, size_t n, size_t offset = 0)
  {
    DeviceAllocator::memcpy(alloc.get_device_ptr() + offset, host_ptr + offset, n);
  }

  /** update from the device, assumes you are copying starting with the device_ptr to the implicit host_ptr.
   */
  static void updateFrom(DualAlloc& alloc, T* host_ptr, size_t n, size_t offset = 0)
  {
    DeviceAllocator::memcpy(host_ptr + offset, alloc.get_device_ptr() + offset, n);
  }

  static void deviceSideCopyN(DualAlloc& alloc, size_t to, size_t n, size_t from)
  {
    T* device_ptr = alloc.get_device_ptr();
    T* to_ptr     = device_ptr + to;
    T* from_ptr   = device_ptr + from;
    DeviceAllocator::memcpy(to_ptr, from_ptr, n);
  }
};

} // namespace qmcplusplus

#endif

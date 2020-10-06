// I need some kind of system for managing storing (caching) chunk
// groups' data on the graphics device. Take a page out of the book of
// CPU caches. This system supports:
//
// 1. Identifying and retrieving entries by a unique ivec3 name.
//
// 2. Reporting failure to retrieve entries (due to not being in the
//    cache) and supporting a FIFO (queue) of entries to be added to
//    the cache (or updated). The queue prevents starvation if there
//    is heavy pressure to add new entries.
//
// 3. Loading new requested entries in the background using worker
//    threads.  This requires "staging buffers" separate from the
//    main cache, so that partially-loaded entries are not visible.
//
// 4. Evicting least-recently-used entries when a new (fully-loaded)
//    entry gets swapped-in from a staging buffer. Entries are
//    assigned to a cache set (entry in 3D cache array) by taking the
//    entry's coordinate modulo the modulus. Each cache set contains
//    AsyncCacheArgs::associativity-many slots; the least-recently
//    used one is evicted.
//
// NOTE: This "management" is really just a high-level scheduler. This
// class never directly interacts with the graphics API; this is
// handled by the actual data types passed in as template parameters
// (which hopefully have useful destructors).

#ifndef MYRICUBE_ASYNCCACHE_HH_
#define MYRICUBE_ASYNCCACHE_HH_

#include "myricube.hh"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <stdio.h>
#include <thread>
#include <utility>
#include <vector>

namespace myricube {

struct AsyncCacheArgs
{
    size_t modulus = 4;
    size_t associativity = 12;
    size_t staging_buffers = 64;
    size_t worker_threads = 2;
    size_t condvar_timeout_ms = 10;
};

// Template class that implements the essay in the header comment.
// However, this is not usable on its own: the user must derive and
// implement the stage and swap_in functions.
//
// bool stage(StagingT*, glm::ivec3)
//
//    Function run on the worker threads. Takes an ivec3 coordinate
//    and fills in the given StagingT with the data corresponding to
//    said coordinate. Returns a success flag; the entry is silently
//    discarded if failure is reported.
//
// void swap_in(StagingT*, std::unique_ptr<EntryT>*)
//
//    Fill the EntryT with data from the given staging buffer
//    (StagingT). Runs on the "main" thread (i.e. the one calling
//    public member functions of AsyncCache).
template <typename StagingT, typename EntryT>
class AsyncCache
{
    // Staging buffer stuff; we need an atomic variable for every
    // staging buffer to coordinate swapping ownership from the "main"
    // thread (the one using the AsyncCache) and the worker threads.
    static constexpr int
        staging_unused = 0,             // Owned by main thread
        staging_unprocessed = 1,        // Owned by worker thread
        staging_ready_to_swap = 2;      // Owned by main thread

    struct StagingBuffer
    {
        StagingT data;
        glm::ivec3 coord;
        std::atomic<int> state = { staging_unused };
    };
    std::unique_ptr<StagingBuffer[]> staging_buffers;
    size_t staging_buffer_count;
    // Used to search quickly for a staging buffer.
    size_t staging_buffer_finger = 0;

    // Queue of coordinates whose entries to load. Push at the end,
    // pop from the front.
    std::vector<glm::ivec3> queue;

    // Main cache: since associativity and modulus are set at runtime,
    // I have to do pointer arithmetic myself on the cache_slots
    // vector. (Basically, it's a 4D array: 3 dimensions for xyz, 1
    // dimension for the cache associativity).
    struct CacheSlot
    {
        // Access time is never set to 0 -- allows me to use 0 as
        // a sentinel value for not-yet-used cache slots. NOTE:
        // entry_ptr need not be non-null for in-use slots!
        uint64_t nonzero_last_access = 0;
        glm::ivec3 coord;
        // Consider exception safety before replacing unique_ptr.
        std::unique_ptr<EntryT> entry_ptr = nullptr;

        bool matches_coord(glm::ivec3 coord_arg) const
        {
            return nonzero_last_access != 0 and coord == coord_arg;
        };

        void update_access_time(uint64_t* counter)
        {
            nonzero_last_access = (*counter)++;
        }
    };

    size_t modulus = 0;
    size_t associativity = 0;
    uint64_t counter = 1;
    std::vector<CacheSlot> cache_slots;

    // Worker threads and stuff to communicate with them: a flag to
    // order them to exit, and a mutex+condvar for waking them up when
    // there's work to do.
    std::unique_ptr<std::thread[]> worker_threads;
    size_t worker_thread_count;
    std::atomic<bool> thread_exit_flag = { false };
    std::mutex condvar_mutex;
    std::condition_variable condvar;
    int condvar_timeout_ms = 0;

    // Extract the bounds of a single cache set corresponding to a
    // given ivec3 coordinate. These bounds are inclusive at the start
    // and exclusive at the end (like typical STL iterators).
    void get_cache_set(
        glm::ivec3 coord,
        CacheSlot** out_begin,
        CacheSlot** out_end)
    {
        return get_cache_set(
            &cache_slots, modulus, associativity, coord, out_begin, out_end);
    }

    // Same as above, but provide the vector, modulus, and
    // associativity explicitly (needed for changing the modulus).
    static void get_cache_set(
        std::vector<CacheSlot>* p_slots,
        size_t modulus_arg,
        size_t assoc,
        glm::ivec3 coord,
        CacheSlot** out_begin,
        CacheSlot** out_end)
    {
        uint32_t modulus = uint32_t(modulus_arg);
        uint32_t x = uint32_t(coord.x ^ 0x8000'0000) % modulus;
        uint32_t y = uint32_t(coord.y ^ 0x8000'0000) % modulus;
        uint32_t z = uint32_t(coord.z ^ 0x8000'0000) % modulus;

        auto cache_set_idx = z*modulus*modulus + y*modulus + x;
        *out_end = &p_slots->at(cache_set_idx * assoc + assoc - 1) + 1;
        *out_begin = *out_end - assoc;
    }

    // Return the youngest (to evict) cache slot in a cache set.  If
    // any not-in-use cache slots are in the cache set, return one of
    // them.
    static CacheSlot* evict_slot(CacheSlot* set_begin, CacheSlot* set_end)
    {
        auto cmp = [](const CacheSlot& left, const CacheSlot& right)
        {
            return left.nonzero_last_access < right.nonzero_last_access;
        };
        return &*std::min_element(set_begin, set_end, cmp);
    }

    // Search for the cache slot corresponding to the entry with the
    // given coordinate. Return whether the entry was found in the
    // cache or not: if so, *out_slot points to that entry, and said
    // entry's access counter is updated; if not, *out_slot points to
    // the entry that may be evicted to make way.
    bool find_access_slot(glm::ivec3 coord, CacheSlot** out_slot)
    {
        CacheSlot* begin;
        CacheSlot* end;
        get_cache_set(coord, &begin, &end);

        for (CacheSlot* slot = begin; slot != end; ++slot) {
            if (slot->matches_coord(coord)) {
                slot->update_access_time(&counter);
                *out_slot = slot;
                return true;
            }
        }

        *out_slot = evict_slot(begin, end);
        return false;
    }

    // Main loop for worker threads. Search for staging buffer entries
    // that need to be loaded and load them. Note that no two worker
    // threads will try to access the staging buffer (see +=
    // worker_thread_count). Sleep and wait for the condvar if nothing
    // seems to be happening.
    void worker_loop_impl(size_t thread_idx)
    {
        while (!thread_exit_flag) {
            size_t processed_count = 0;

            for (size_t i = thread_idx;
                i < staging_buffer_count;
                i += worker_thread_count) {

                StagingBuffer* buffer = &staging_buffers[i];
                if (buffer->state == staging_unprocessed) {
                    if (thread_exit_flag) return;

                    stage(&buffer->data, buffer->coord);
                    buffer->state = staging_ready_to_swap;
                    ++processed_count;
                }
            }

            if (processed_count == 0) {
                if (thread_exit_flag) return;
                std::unique_lock<std::mutex> lock(condvar_mutex);
                auto condvar_timeout =
                    std::chrono::milliseconds(condvar_timeout_ms);
                condvar.wait_for(lock, condvar_timeout);
            }
        }
    }

    static void worker_loop(AsyncCache* self, size_t thread_idx)
    {
        self->worker_loop_impl(thread_idx);
    }

  public:
    AsyncCache(AsyncCacheArgs args)
    {
        assert(args.staging_buffers > 0);
        staging_buffer_count = args.staging_buffers;
        staging_buffers.reset(new StagingBuffer[staging_buffer_count]);

        associativity = args.associativity;
        modulus = args.modulus;
        assert(associativity > 0 and modulus > 0);
        cache_slots.resize(associativity * modulus * modulus * modulus);

        condvar_timeout_ms = args.condvar_timeout_ms;
        worker_thread_count = args.worker_threads;
        assert(args.worker_threads > 0);

        // Now that everything is set up, we can go ahead and launch
        // the worker threads with pointers back to this object.
        worker_threads.reset(new std::thread[worker_thread_count]);
        for (size_t i = 0; i < worker_thread_count; ++i) {
            worker_threads[i] = std::thread(worker_loop, this, i);
        }
    };

    virtual ~AsyncCache()
    {
        thread_exit_flag.store(true);
        for (size_t i = 0; i < worker_thread_count; ++i) {
            worker_threads[i].join();
        }
    }

    AsyncCache(AsyncCache&&) = delete;

    // Resize the cache (if new_modulus differs from the current
    // modulus), and move all the entries over (discarding some if
    // they don't fit).
    void set_modulus(size_t new_modulus)
    {
        if (new_modulus = modulus) return;
        assert(new_modulus > 0);

        std::vector<CacheSlot> new_cache_slots(
            new_modulus * new_modulus * new_modulus * associativity);

        // Just put every valid entry in the correct location in the
        // new cache.
        for (CacheSlot& slot : cache_slots) {
            CacheSlot* set_begin;
            CacheSlot* set_end;
            if (slot.nonzero_last_access == 0) continue;

            get_cache_set(
                &new_cache_slots,
                new_modulus,
                associativity,
                slot.coord, &set_begin, &set_end);
            CacheSlot* to_evict = evict_slot(set_begin, set_end);
            *to_evict = std::move(slot);
        }

        // Only modify the class state when we're guaranteed to succeed.
        cache_slots = std::move(new_cache_slots);
        modulus = new_modulus;
    };

    // Push a coordinate onto the FIFO of entries to load/update.
    // This push is not done if the coordinate is already in the
    // queue.
    void enqueue(glm::ivec3 coord)
    {
        for (glm::ivec3 queue_coord : queue) {
            if (coord == queue_coord) return;
        }
        queue.push_back(coord);
    }

    // Step 1/2 of populating the cache: Pop at most max_stage
    // coordinates from the queue and put them into staging
    // buffers. The worker threads will eventually get to loading the
    // relevant data in. Return the number actually sent to staging
    // buffers.
    size_t stage_from_queue(size_t max_stage=1) noexcept
    {
        size_t pop_count = 0;
        for (size_t i = 0; i < max_stage; ++i, ++pop_count) {
            if (pop_count >= queue.size()) break;

            // Look for an available staging buffer.
            StagingBuffer* p_buffer = nullptr;
            for (size_t j = 0; j < staging_buffer_count; ++j) {
                p_buffer = &staging_buffers[staging_buffer_finger];
                if (p_buffer->state == staging_unused) goto found;

                staging_buffer_finger =
                    staging_buffer_finger == 0 ?
                    staging_buffer_count - 1 : staging_buffer_finger - 1;
            }
            // Failure case.
            break;

            // Success -- assign coordinate to available staging buffer.
          found:
            p_buffer->coord = queue[pop_count];
            p_buffer->state = staging_unprocessed;
        }
        if (pop_count > 0) condvar.notify_all();
        queue.erase(queue.begin(), queue.begin() + pop_count);
        return pop_count;
    }

    // Step 2/2 of populating the cache: Move at most max_swap
    // fully-loaded entries from the staging buffers to the main
    // cache, evicting old stuff if needed. Return the actual number
    // of entries swapped-in.
    size_t swap_in_from_staging(size_t max_swap=1)
    {
        size_t return_count = 0;
        for (size_t i = 0; i < staging_buffer_count; ++i) {
            if (return_count >= max_swap) break;

            StagingBuffer& buffer = staging_buffers[i];
            if (buffer.state == staging_ready_to_swap) {
                CacheSlot* slot;
                find_access_slot(buffer.coord, &slot);

                // Exception safety: mark incompletely-initialized
                // slots as not used.
                slot->nonzero_last_access = 0;
                swap_in(&buffer.data, &slot->entry_ptr);
                slot->coord = buffer.coord;
                slot->update_access_time(&counter);

                buffer.state = staging_unused;
                ++return_count;
            }
        }
        return return_count;
    };

    // Search for the cache entry with the given 3D
    // coordinate. Returns nullptr if not in the cache, or the entry
    // stores a null pointer.
    EntryT* request_entry(glm::ivec3 coord)
    {
        CacheSlot* cache_slot;
        bool found = find_access_slot(coord, &cache_slot);
        return found ? cache_slot->entry_ptr.get() : nullptr;
    };

    void debug_dump()
    {
        fprintf(stderr, "================\n");
        for (size_t i = 0; i < staging_buffer_count; ++i) {
            StagingBuffer& buffer = staging_buffers[i];
            fprintf(stderr, "(%i %i %i) state=%i\n",
                buffer.coord.x, buffer.coord.y, buffer.coord.z,
                buffer.state.load());
        }
        fprintf(stderr, "================\n");
        for (glm::ivec3 coord : queue) {
            fprintf(stderr, "(%i %i %i)\n", coord.x, coord.y, coord.z);
        }
        fprintf(stderr, "================\n");

        size_t idx = 0;
        for (CacheSlot& slot : cache_slots) {
            if (idx++ % associativity == 0) {
                fprintf(stderr, "---------------\n");
            }
            fprintf(stderr, "%s (%i %i %i) %li %p\x1b[0m\n",
                slot.nonzero_last_access == 0 ? "" : "\x1b[32m",
                slot.coord.x, slot.coord.y, slot.coord.z,
                long(slot.nonzero_last_access),
                slot.entry_ptr.get());
        }
        fprintf(stderr, "---------------\n");
    }

  protected:
    virtual bool stage(StagingT*, glm::ivec3) = 0;
    virtual void swap_in(StagingT*, std::unique_ptr<EntryT>*) = 0;
};

} // end namespace myricube

#endif /* !MYRICUBE_ASYNCCACHE_HH_ */

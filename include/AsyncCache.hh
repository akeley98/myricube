// I need some kind of system for managing storing (caching) chunk
// groups' data on the graphics device. Take a page out of the book of
// CPU caches. This system supports:
//
// 1. Identifying and retrieving entries by a unique ivec3
//    name. [request_entry, which returns EntryT*]
//
// 2. Reporting failure to retrieve entries (due to not being in the
//    cache) and supporting a FIFO (queue) of entries to be added to
//    the cache (or updated). The queue prevents starvation if there
//    is heavy pressure to add new entries. [enqueue]
//
// 3. Loading new requested entries in the background using worker
//    threads.  This requires "staging buffers" separate from the
//    main cache, so that partially-loaded entries are not visible.
//    [StagingT]
//
// 4. Tracking the number of frames elapsed, and recording the last
//    frame in which an entry was used. [Call begin_frame at the start
//    of each frame. request_entry records frame number when an EntryT
//    is accessed.]
//
// 5. Evicting least-recently-used entries when a new (fully-loaded)
//    entry gets swapped-in from a staging buffer. Entries are
//    assigned to a cache set (entry in 3D cache array) by taking the
//    entry's coordinate modulo the modulus. Each cache set contains
//    AsyncCacheArgs::associativity-many slots; the least-recently
//    used one is evicted.
//
// 6. Preventing entries from being evicted until the frame number is
//    at least safety_frame_count greater than the frame number in
//    which the entry was last used; this is needed to prevent
//    still-in-use resources from being deleted, or concurrently read
//    and written.
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
    uint64_t safety_frame_count = 3;
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
// bool swap_in(StagingT*, std::unique_ptr<EntryT>*)
//
//    Fill the EntryT with data from the given staging buffer
//    (StagingT). Runs on the "main" thread (i.e. the one calling
//    public member functions of AsyncCache).
//
//    Also returns a success flag, however, instead of ignoring
//    failures, failures are retried later.
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
    size_t staging_buffer_queue_finger = 0;
    size_t staging_buffer_swap_finger = 0;

    // Queue of coordinates whose entries to load. Push at the end,
    // pop from the front.
    std::vector<glm::ivec3> queue;

    // Main cache: since associativity and modulus are set at runtime,
    // I have to do pointer arithmetic myself on the cache_slots
    // vector. (Basically, it's a 4D array: 3 dimensions for xyz, 1
    // dimension for the cache associativity).
    struct CacheSlot
    {
        // frame_counter when this entry was last accessed.
        // 0 if never used.
        uint64_t last_access_frame = 0;

        // Whether this slot corresponds to an entry for an actual
        // chunk group. NOTE: It is still not safe to evict an invalid
        // entry if it is labeled as recently used!
        bool valid = false;

        // Group coordinate of the chunk group whose data is stored here.
        glm::ivec3 coord;

        // Consider exception safety before replacing unique_ptr.
        std::unique_ptr<EntryT> entry_ptr = nullptr;

        bool matches_coord(glm::ivec3 coord_arg) const
        {
            return valid and coord == coord_arg;
        };
    };

    size_t modulus = 0;
    size_t associativity = 0;
    uint64_t frame_counter = 0;
    uint64_t safety_frame_count = 0;
    std::vector<CacheSlot> cache_slots;
    CacheSlot extra_slot;

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
    // them. If none of them can be evicted due to safety_frame_count,
    // return nullptr.
    CacheSlot* evict_slot(CacheSlot* set_begin, CacheSlot* set_end)
    {
        auto cmp = [](const CacheSlot& left, const CacheSlot& right)
        {
            return left.last_access_frame < right.last_access_frame;
        };
        CacheSlot* to_evict = &*std::min_element(set_begin, set_end, cmp);
        bool can_evict =
            to_evict->last_access_frame == 0 or
            to_evict->last_access_frame + safety_frame_count <= frame_counter;
        return can_evict ? to_evict : nullptr;
    }

    // Search for the cache slot corresponding to the entry with the
    // given coordinate. If found, mark it as recently-accessed (with
    // the current frame number). Otherwise, return nullptr.
    CacheSlot* read_access_slot(glm::ivec3 coord)
    {
        // call begin_frame() before populating the cache.
        assert(frame_counter != 0);

        CacheSlot* begin;
        CacheSlot* end;
        get_cache_set(coord, &begin, &end);

        for (CacheSlot* slot = begin; slot != end; ++slot) {
            if (slot->matches_coord(coord)) {
                slot->last_access_frame = frame_counter;
                return slot;
            }
        }
        return nullptr;
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

                    bool success = stage(&buffer->data, buffer->coord);
                    buffer->state = success
                                  ? staging_ready_to_swap : staging_unused;
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

        safety_frame_count = args.safety_frame_count;

        // Now that everything is set up, we can go ahead and launch
        // the worker threads with pointers back to this object.
        worker_threads.reset(new std::thread[worker_thread_count]);
        for (size_t i = 0; i < worker_thread_count; ++i) {
            worker_threads[i] = std::thread(worker_loop, this, i);
        }
    };

    virtual ~AsyncCache()
    {
        stop_wait_threads();
    }

    void stop_wait_threads()
    {
        thread_exit_flag.store(true);
        condvar.notify_all();
        for (size_t i = 0; i < worker_thread_count; ++i) {
            std::thread& thr = worker_threads[i];
            if (thr.joinable()) thr.join();
        }
    }

    AsyncCache(AsyncCache&&) = delete;

    // Update the frame counter; must be the first member function
    // called on the AsyncCache (aside from constructor).
    void begin_frame()
    {
        ++frame_counter;
    }

    // Resize the cache (if new_modulus differs from the current
    // modulus), and move all the entries over (discarding some if
    // they don't fit).
    //
    // We don't honor the safety_frame_count here. To avoid invalid
    // memory troubles, you need to pass a callback that stalls the
    // pipeline (glFinish or vkDeviceWaitIdle). CRUCIALLY, this only
    // affects the main cache (of EntryT), which the worker threads
    // can't access; the thread-shared staging buffer is not affected.
    template <typename GLFinish>
    void set_modulus(size_t new_modulus, GLFinish finish)
    {
        if (new_modulus == modulus) return;
        finish();
        assert(new_modulus > 0);

        std::vector<CacheSlot> new_cache_slots(
            new_modulus * new_modulus * new_modulus * associativity);

        // Just put every valid entry in the correct location in the
        // new cache.
        for (CacheSlot& slot : cache_slots) {
            CacheSlot* set_begin;
            CacheSlot* set_end;
            if (!slot.valid) continue;

            get_cache_set(
                &new_cache_slots,
                new_modulus,
                associativity,
                slot.coord, &set_begin, &set_end);
            CacheSlot* to_evict = evict_slot(set_begin, set_end);
            if (to_evict == nullptr) continue;
            *to_evict = std::move(slot);
        }

        // Only modify the class state when we're guaranteed to succeed.
        cache_slots = std::move(new_cache_slots);
        modulus = new_modulus;
    };

    // Push a coordinate onto the FIFO of entries to load/update.
    // This push is not done if the coordinate is already in the
    // queue.
    //
    // TODO: I think this queue could use some improvement: make the
    // linear search better, and find a way to cancel staging
    // something if it doesn't seem to really matter anymore.
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
    //
    // Optional Acceptor is run on every glm::ivec3 coordinate popped
    // from the queue: if Acceptor returns false, the coordinate is
    // silently discarded without being staged.
    size_t stage_from_queue(size_t max_stage=1) noexcept
    {
        auto always_accept = [] (glm::ivec3) { return true; };
        return stage_from_queue(max_stage, always_accept);
    }

    template <typename Acceptor>
    size_t stage_from_queue(size_t max_stage, Acceptor&& acceptor) noexcept
    {
        size_t pop_count = 0;
        size_t staged_count = 0;
        for (; staged_count < max_stage; ++pop_count) {
            if (pop_count >= queue.size()) break;

            // Check the next group coordinate in the queue; skip if
            // acceptor doesn't like it.
            glm::ivec3 coord = queue[pop_count];
            if (!acceptor(coord)) continue;

            // Look for an available staging buffer.
            StagingBuffer* p_buffer = nullptr;
            for (size_t j = 0; j < staging_buffer_count; ++j) {
                p_buffer = &staging_buffers[staging_buffer_queue_finger];
                if (p_buffer->state == staging_unused) goto found;

                staging_buffer_queue_finger =
                    staging_buffer_queue_finger == 0 ?
                    staging_buffer_count - 1 : staging_buffer_queue_finger - 1;
            }
            // Failure case.
            break;

            // Success -- assign coordinate to available staging buffer.
          found:
            p_buffer->coord = coord;
            p_buffer->state = staging_unprocessed;
            ++staged_count;
        }
        if (staged_count > 0) condvar.notify_all();
        queue.erase(queue.begin(), queue.begin() + pop_count);
        return staged_count;
    }

    // Step 2/2 of populating the cache: Make at most max_swap
    // attempts to move fully-loaded entries from the staging buffers
    // to the main cache, evicting old stuff if needed. Return the
    // actual number of entries swapped-in.
    size_t swap_in_from_staging(size_t max_attempts=1)
    {
        size_t attempts = 0;
        size_t swapped_count = 0;

        // This is a little weird, basically, I want to use the finger
        // to make sure all parts of the staging buffers array are
        // visited at some point, even if max_attempts is small, but I
        // also want to revisit staging buffer entries whose swap_in
        // failed ASAP.
        //
        // So, I store the index (in staging_buffers) where the first
        // failed swap_in occured, and reset it to here if it happens.
        size_t finger_to_store = size_t(-1);

        for (size_t i = 0; i < staging_buffer_count; ++i) {
            if (attempts >= max_attempts) break;
            StagingBuffer& buffer = staging_buffers[staging_buffer_swap_finger];

            if (buffer.state == staging_ready_to_swap) {
                attempts++;

                // First, find the correct cache set and an entry to evict.
                CacheSlot* set_begin;
                CacheSlot* set_end;
                get_cache_set(buffer.coord, &set_begin, &set_end);
                CacheSlot* to_evict = evict_slot(set_begin, set_end);

                // Short circuit evaluation needed here.
                // Try to swap in only if we actually have room to store
                // the newly swapped-in entry. Use the extra cache slot.
                // Avoid modifying main cache state until success is certain.
                bool success =
                    to_evict != nullptr and
                    swap_in(&buffer.data, &extra_slot.entry_ptr);


                if (!success) {
                    // Try again soon later: reset the finger to here
                    // if this was the first failure.
                    if (finger_to_store != size_t(-1)) {
                        finger_to_store = staging_buffer_swap_finger;
                        // TODO
                    }
                    continue;
                }

                buffer.state = staging_unused;

                // This shouldn't emit any exceptions.
                // Mark the slot with the coordinate and frame number,
                // then swap it into the evicted cache slot, BUT, first
                // invalidate any existing cache slot for this chunk
                // group, so that the most up-to-date entry is used.
                for (CacheSlot* s = set_begin; s != set_end; ++s) {
                    if (s->matches_coord(buffer.coord)) {
                        s->valid = false;
                    }
                }
                extra_slot.coord = buffer.coord;
                extra_slot.last_access_frame = frame_counter;
                extra_slot.valid = true;
                std::swap(extra_slot, *to_evict);

                ++swapped_count;
            }

            staging_buffer_swap_finger =
                staging_buffer_swap_finger == 0 ?
                staging_buffer_count - 1 : staging_buffer_swap_finger - 1;
        }

        if (finger_to_store != size_t(-1)) {
            staging_buffer_swap_finger = finger_to_store;
        }

        return swapped_count;
    };

    // Search for the cache entry with the given 3D
    // coordinate. Returns nullptr if not in the cache, or the entry
    // stores a null pointer.
    EntryT* request_entry(glm::ivec3 coord)
    {
        CacheSlot* cache_slot = read_access_slot(coord);
        return cache_slot ? cache_slot->entry_ptr.get() : nullptr;
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
                slot.valid ? "" : "\x1b[32m",
                slot.coord.x, slot.coord.y, slot.coord.z,
                long(slot.last_access_frame),
                slot.entry_ptr.get());
        }
        fprintf(stderr, "---------------\n");
    }

  protected:
    virtual bool stage(StagingT*, glm::ivec3) = 0;
    virtual bool swap_in(StagingT*, std::unique_ptr<EntryT>*) = 0;
};

} // end namespace myricube

#endif /* !MYRICUBE_ASYNCCACHE_HH_ */

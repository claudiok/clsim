/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id$
 *
 * @file I3CLSimStepStore.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMSTEPSTORE_H_INCLUDED
#define I3CLSIMSTEPSTORE_H_INCLUDED

/**
 * The I3CLSimStepStore class stores I3CLSimStep objects
 * indexed by a photon multiplicity. Arbitrarily sized bunches 
 * of steps can be retrieved. They will be clustered by
 * multiplicity.
 */

#include "icetray/I3TrayHeaders.h"
#include "clsim/I3CLSimStep.h"

#include <stdint.h>

#include <vector>
#include <deque>
#include <limits>

#include <boost/static_assert.hpp>
#include <boost/foreach.hpp>
#include <boost/pool/pool_alloc.hpp> 

template <typename U, class T>
class I3CLSimTemplateStore 
{
private:
    typedef std::deque<T> TDequeType;
    
    /// @cond this assert confuses doxygen
    // static_assert: U==unsigned integer (8,16,32 or 64 bit)
    BOOST_STATIC_ASSERT((std::numeric_limits<U>::digits >= 8)
                        && std::numeric_limits<U>::is_specialized
                        && std::numeric_limits<U>::is_integer
                        && !std::numeric_limits<U>::is_signed);
    /// @endcond

    // static_assert: max<U> <= max<std::size_t>
    BOOST_STATIC_ASSERT((std::numeric_limits<U>::digits <= std::numeric_limits<std::size_t>::digits));
    
public:
    I3CLSimTemplateStore(std::size_t initialSize):
    bins_(initialSize, NULL),
    currentSize_(0)
    {
        for (std::size_t i=0;i<initialSize;++i)
        {
            bins_[i] = new TDequeType();
        }
    }

    I3CLSimTemplateStore():
    currentSize_(0)
    {
    }
    
    ~I3CLSimTemplateStore()
    {
        BOOST_FOREACH(TDequeType *ptr, bins_)
        {
            delete ptr; // NULL-pointer deletion is valid (and a no-op)
        }
        
    }
    
    /**
     * inserts a copy of 
     */
    inline void insert_copy(U index, const T &value)
    {
        insert_new(index) = value;
    }
    
    /**
     * inserts a default-constructed instance of the class at a certain index and
     * returns a reference to it
     */
    inline T &insert_new(U index)
    {
        // re-size the number of bins if necessary
        if (index >= bins_.size()) 
        {
            const std::size_t oldSize = bins_.size();
            const std::size_t newSize = static_cast<std::size_t>(index+1);

            bins_.resize(newSize, NULL);
            
            for (std::size_t i = oldSize; i < newSize; ++i)
            {
                bins_[i] = new TDequeType();
            }
        }
        
        bins_[index]->push_back(T());
        ++currentSize_;
        
        return bins_[index]->back();
    }
    
    inline std::size_t size() const
    {
        return currentSize_;
    }

    inline bool empty() const
    {
        return (currentSize_==0);
    }
    
    inline TDequeType &operator[](U index)
    {
        return *(bins_[index]);
    }

    inline const TDequeType &operator[](U index) const
    {
        return *(bins_[index]);
    }
    
    inline TDequeType &at(U index)
    {
        return *(bins_.at(index));
    }

    inline const TDequeType &at(U index) const
    {
        return *(bins_.at(index));
    }
    
    /**
     * takes a number of entries, copies them into a vector
     * (sorted by index) and pops them from this container.
     * If less than the specified number of entries exist in
     * this container, it is fully emptied.
     *
     * All current entries in the vector are removed.
     */
    inline void pop_bunch_to_vector(std::size_t size, std::vector<T> &vect)
    {
        const std::size_t realSize = std::min(size, currentSize_);
        vect.clear();
        if (realSize==0) return;

        std::size_t itemsPopped=0;
        
        BOOST_FOREACH(TDequeType *ptr, bins_)
        {
            TDequeType &currentQueue = *ptr;

            for (;;)
            {
                if (currentQueue.empty()) break;

                vect.push_back(currentQueue.front());   // insert into vector
                currentQueue.pop_front();               // remove from deque
                ++itemsPopped;
                
                if (itemsPopped>=realSize) break; // are we finished yet?
            }
            if (itemsPopped>=realSize) break; // are we finished yet?
        }
        
        if (itemsPopped > currentSize_)
        {
            // something is seriously wrong
            log_fatal("Internal implementation error.");
        }
        
        currentSize_ -= itemsPopped;
    }

    /**
     * takes a number of entries, copies them into a vector
     * (sorted by index) and pops them from this container.
     * If less than the specified number of entries exist in
     * this container, the remaining entries are filled with
     * copies of a template.
     *
     * All current entries in the vector are removed.
     */
    inline void pop_bunch_to_vector(std::size_t size, std::vector<T> &vect, const T &temp)
    {
        vect.clear();
        vect.reserve(size);
        pop_bunch_to_vector(size, vect);
        const std::size_t realEntries = vect.size();

        // fill the remainder of the vector with copies
        // of the template
        for (std::size_t i=realEntries; i<size; ++i)
        {
            vect.push_back(temp);
        }
    }
    
private:
    std::vector<TDequeType *> bins_;
    std::size_t currentSize_;
};

class I3CLSimStepStore : public I3CLSimTemplateStore<uint32_t, I3CLSimStep>
{
    typedef I3CLSimStep T;
    typedef uint32_t U;
    typedef I3CLSimTemplateStore<U,T> Base;
    using I3CLSimTemplateStore<U,T>::I3CLSimTemplateStore;

public:
    /**
     * inserts a copy of 
     */
    inline void insert_copy(U index, const T &value)
    {
        insert_new(index, value.GetID()) = value;
    }
    
    /**
     * inserts a default-constructed instance of the class at a certain index and
     * returns a reference to it
     */
    inline T &insert_new(U index, uint32_t id)
    {
        pendingIds_.insert(std::make_pair(id, 0)).first->second++;
        return Base::insert_new(index);
    }

    
    inline void pop_bunch_to_vector(std::size_t size, std::vector<T> &vect)
    {
        Base::pop_bunch_to_vector(size, vect);
        for (const I3CLSimStep &step : vect) {
            auto count = pendingIds_.find(step.GetID());
            assert( count != pendingIds_.end() );
            if (--(count->second) == 0) {
                pendingIds_.erase(count);
            }
        }
    }
    
    inline void pop_bunch_to_vector(std::size_t size, std::vector<T> &vect, const T &temp)
    {
        pop_bunch_to_vector(size, vect);
        const std::size_t realEntries = vect.size();

        // fill the remainder of the vector with copies
        // of the template
        for (std::size_t i=realEntries; i<size; ++i)
        {
            vect.push_back(temp);
        }
    }
    
    uint32_t count(uint32_t identifier) const
    {
        auto it = pendingIds_.find(identifier);
        return it == pendingIds_.end() ? 0 : it->second;
    }

private:
    std::map<uint32_t, uint32_t> pendingIds_;
};

// typedef I3CLSimTemplateStore<uint32_t, I3CLSimStep> I3CLSimStepStore;

I3_POINTER_TYPEDEFS(I3CLSimStepStore);


#endif //I3CLSIMSTEPSTORE_H_INCLUDED

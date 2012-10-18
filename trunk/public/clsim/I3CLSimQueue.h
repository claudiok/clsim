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
 * @file I3CLSimQueue.h
 * @version $Revision$
 * @date $Date$
 * @author Claudio Kopper
 */

#ifndef I3CLSIMQUEUE_H_INCLUDED
#define I3CLSIMQUEUE_H_INCLUDED

#include <queue>

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

/**
 * @brief A thread-safe queue, storing objects of type T.
 * T will be copied around quite a bit, so make it light-weight.
 * (Use a shared pointer for example.)
 * 
 * Will block on Get if queue is empty and on Put if queue is full.
 * A max_size argument of 0 will get you a queue with no limit.
 */

template <typename T>
class I3CLSimQueue : private boost::noncopyable
{
public:
    I3CLSimQueue(std::size_t max_size) 
    : max_size_(max_size) {}
    
    I3CLSimQueue() {;}
    
    ~I3CLSimQueue() {;}

    void Put(const T &msg)
    {
        // lock the mutex to ensure exclusive access to the queue
        boost::unique_lock<boost::mutex> guard(mutex_);
        
        // as long as the queue is full, wait until something is taken off of it
        while ((max_size_ > 0) && (queue_.size() >= max_size_))
        {
            cond_.wait(guard);
        }

        // add the message to the queue
        queue_.push(msg);
        
        // notify the consumer thread
        cond_.notify_one();
    }
    
    
    T Get()
    {
        // lock the mutex to ensure exclusive access to the queue
        boost::unique_lock<boost::mutex> guard(mutex_);
        
        // in case the queue is empty, sleep waiting for something to be put onto it
        while (queue_.empty())
        {
            cond_.wait(guard);
        }
        
        // the queue is not empty anymore, read the value
        T msg = queue_.front();
        
        // remove the current message from the queue
        queue_.pop();
        
        // notify the producer that there is space on the queue now
        cond_.notify_one();
        
        return msg;
    }

    bool GetNonBlocking(T &value)
    {
        // lock the mutex to ensure exclusive access to the queue
        boost::unique_lock<boost::mutex> guard(mutex_);
        
        if (queue_.empty())
            return false;
        
        // the queue is not empty anymore, read the value
        value = queue_.front();
        
        // remove the current message from the queue
        queue_.pop();
        
        // notify the producer that there is space on the queue now
        cond_.notify_one();
        
        return true;
    }

    T Get(double timeout, T returnOnTimeout) // timeout in seconds
    {
        // lock the mutex to ensure exclusive access to the queue
        boost::unique_lock<boost::mutex> guard(mutex_);
        
        // in case the queue is empty, sleep waiting for something to be put onto it
        while (queue_.empty())
        {
            boost::posix_time::time_duration td = boost::posix_time::milliseconds(static_cast<long>(timeout*1000.));
            bool ret = cond_.timed_wait(guard, td);
            
            if (!ret) {
                // timeout reached, return dummy
                return returnOnTimeout;
            }
        }
        
        // the queue is not empty anymore, read the value
        T msg = queue_.front();
        
        // remove the current message from the queue
        queue_.pop();
        
        // notify the producer that there is space on the queue now
        cond_.notify_one();
        
        return msg;
    }

    bool empty() const
    {
        // lock the mutex to ensure exclusive access to the queue
        boost::unique_lock<boost::mutex> guard(mutex_);
        
        return queue_.empty();
    }

    std::size_t size() const
    {
        // lock the mutex to ensure exclusive access to the queue
        boost::unique_lock<boost::mutex> guard(mutex_);
        
        return queue_.size();
    }

    inline std::size_t max_size() const
    {
        return max_size_;
    }

private:
    mutable boost::mutex mutex_;
    boost::condition_variable_any cond_;
    std::queue<T> queue_;
    std::size_t max_size_;
};


#endif //I3CLSIMQUEUE_H_INCLUDED

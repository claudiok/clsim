
#include <boost/date_time/posix_time/posix_time.hpp>

#include "clsim/I3CLSimServer.h"

#include <list>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/math/common_factor_rt.hpp>


namespace {

std::tuple<zmq::message_t, std::list<zmq::message_t> >
read_message(zmq::socket_t &socket)
{
    zmq::message_t address;
    std::list<zmq::message_t> body;
    
    socket.recv(&address);
    if (!address.more()) {
        body.push_back(std::move(address));
        address.rebuild();
        return std::make_tuple(std::move(address), std::move(body));
    }
    
    do {
        zmq::message_t part;
        socket.recv(&part);
        body.push_back(std::move(part));
    } while ((--body.end())->more());
    
    // remove the delimiter added by REQ or REP sockets
    if (body.front().size() == 0)
        body.pop_front();
        
    return std::make_tuple(std::move(address), std::move(body));
}

// NB: for simple struct-like types such as I3CLSimStep and I3CLSimPhoton, serializing
// to a properly-sized vector<char> is only a few percent slower than a straight memcpy().
// The extra copy to the output message_t 
template<typename T>
zmq::message_t serialize(const T &value)
{
    std::vector<char> output;
    {
        boost::iostreams::filtering_ostream ostream;
        ostream.push(boost::iostreams::back_inserter(output));
        icecube::archive::portable_binary_oarchive poa(ostream);

        poa << value;
    }

    return zmq::message_t(output.data(), output.size());
}

template <typename T>
T deserialize(const zmq::message_t &message)
{
    T value;

    boost::iostreams::filtering_istream istream;
    istream.push(boost::iostreams::array_source(message.data<char>(), message.size()));
    icecube::archive::portable_binary_iarchive pia(istream);

    pia >> value;

    return std::move(value);
}

}

I3CLSimServer::I3CLSimServer(const std::string &address, const std::vector<I3CLSimStepToPhotonConverterPtr> &converters) :
    converters_(converters),
    workgroupSize_(0), maxBunchSize_(0),
    context_(1),
    frontend_(context_, ZMQ_ROUTER),
    backend_(context_, ZMQ_ROUTER),
    control_(context_, ZMQ_PUB),
    heartbeat_(context_, ZMQ_REP)
{
    frontend_.bind(address);
    backend_.bind("inproc://worker");
    control_.bind("inproc://control");
    heartbeat_.bind("inproc://heartbeat");
    
    if (converters_.empty())
        log_fatal("Need at least 1 I3CLSimStepToPhotonConverter");

    // Harmonize bunch sizes
    for (auto &converter : converters_) {
        if (!converter || !converter->IsInitialized())
            log_fatal("All I3CLSimStepToPhotonConverters must be initialized");
        
        if (workgroupSize_ == 0) {
            workgroupSize_ = converter->GetWorkgroupSize();
        } else {
            workgroupSize_ = boost::math::lcm(workgroupSize_, converter->GetWorkgroupSize());
        }
        
        if (maxBunchSize_ == 0) {
            maxBunchSize_ = converter->GetMaxNumWorkitems();
        } else {
            std::size_t newMaxBunchSize = std::min(maxBunchSize_, converter->GetMaxNumWorkitems());
            std::size_t newMaxBunchSizeWithGranularity = newMaxBunchSize - newMaxBunchSize%workgroupSize_;
            
            if (newMaxBunchSizeWithGranularity == 0)
                log_fatal("maximum bunch sizes are incompatible with kernel work group sizes.");
            maxBunchSize_ = newMaxBunchSizeWithGranularity;
        }
    }
    
    serverThread_ = std::thread(&I3CLSimServer::ServerThread, this);
    {
        // Wait for thread to come online so that shutdown messages
        // in the destructor will have an effect
        zmq::message_t ping;
        heartbeat_.recv(&ping);
        heartbeat_.send(ping);
    }
    
    unsigned queueDepth = 5;

    for (unsigned i=0; i < converters_.size(); i++)
        for (unsigned j=0; j < queueDepth; j++) {
            workerThreads_.emplace_back(&I3CLSimServer::WorkerThread, this, i, j);
            
            zmq::message_t ping;
            heartbeat_.recv(&ping);
            heartbeat_.send(ping);
        }
    
}

I3CLSimServer::~I3CLSimServer()
{
    // send shutdown message to threads
    try {
        control_.send(zmq::message_t());
    } catch (...) {}
    log_debug("Sent shutdown message");
    serverThread_.join();
    for(auto &thread : workerThreads_)
        thread.join();
}

template<typename T>
inline bool send_msg(
  zmq::socket_t &socket,
  const T& payload,
  bool will_continue = false
)
{
  zmq::message_t msg_payload( sizeof(T) );
  std::memcpy(msg_payload.data(), &payload, sizeof(T) );
  return socket.send(msg_payload, will_continue ? ZMQ_SNDMORE : 0); // last message part
}

void I3CLSimServer::ServerThread()
{
    // Listen for control messages
    zmq::socket_t control(context_, ZMQ_SUB);
    control.connect("inproc://control");
    control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    // Signal main thread that we're alive
    zmq::socket_t heartbeat(context_, ZMQ_REQ);
    heartbeat.connect("inproc://heartbeat");
    {
        zmq::message_t ping;
        heartbeat.send(ping);
        heartbeat.recv(&ping);
    }
    
    zmq::pollitem_t items[] = {
      { static_cast<void *>(frontend_), 0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(backend_),  0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(control),   0, ZMQ_POLLIN, 0 },
    };
    
    log_trace("Server thread started");
    
    while (true) {
        
        if (int rc = zmq::poll(&items[0], 3, 0) < -1) {
            log_error_stream("ZMQ polling error " << rc);
            break;
        }
        
        
        // Message from client
        // Process only if there are available workers
        if ((items[0].revents & ZMQ_POLLIN) && workers_.size() > 0) {
            
            zmq::message_t address;
            std::list<zmq::message_t> body;
            std::tie(address, body) = read_message(frontend_);
            
            if ( body.size() == 2 ) {
                // assign an internal ID for later reply to client
                uint32_t client_id(0);
                if (clients_.size() > 0)
                    client_id = (--clients_.end())->first+1;
                if (clients_.find(client_id) != clients_.end())
                    log_fatal("Repeated client ID");
                clients_.emplace(client_id,
                    std::array<zmq::message_t,2>({{std::move(address), std::move(body.back())}}));
                // drop external ID
                body.pop_back();
            
                // forward to a worker
                log_trace_stream("Forwarding "<<(body.size())<<" packets to worker");
                backend_.send(workers_.front(), ZMQ_SNDMORE);
                backend_.send(zmq::message_t(), ZMQ_SNDMORE);
                for (auto &packet : body)
                    backend_.send(packet, ZMQ_SNDMORE);
                backend_.send(serialize(client_id), 0);
            
                workers_.pop();
            } else if ( body.size() == 1 ) {
                // a client has just connected. send the desired workgroup batching
                frontend_.send(address, ZMQ_SNDMORE);
                frontend_.send(serialize(std::make_pair(workgroupSize_, maxBunchSize_)), 0);
            } else {
                log_fatal("Unknown message format!");
            }
        }
        
        // Message from worker
        if (items[1].revents & ZMQ_POLLIN) {
            
            zmq::message_t address;
            std::list<zmq::message_t> body;
            std::tie(address, body) = read_message(backend_);
            
            workers_.push(std::move(address));
            
            // A "ready" message; no payload
            if (body.size() == 1 && body.front().size() == 0)
                continue;
            
            // Otherwise, look up address of originating client
            // and return the payload
            uint32_t client_id = deserialize<uint32_t>(body.back());
            body.pop_back();
            auto destination = clients_.find(client_id);
            if (destination == clients_.end()) {
                log_error_stream("Unknown client ID " << client_id);
                continue;
            }
            
            frontend_.send(destination->second[0], ZMQ_SNDMORE);
            frontend_.send(zmq::message_t(), ZMQ_SNDMORE);
            for (auto &packet : body)
                frontend_.send(packet, ZMQ_SNDMORE);
            frontend_.send(destination->second[1], 0);
            
            clients_.erase(destination);
        }
        
        // Shutdown signal from main thread
        if (items[2].revents & ZMQ_POLLIN) {
            zmq::message_t dummy;
            control.recv(&dummy);
            log_debug("Server thread shutting down");
            break;
        }
        

    }
}

void I3CLSimServer::WorkerThread(unsigned index, unsigned buffer_id)
{
    // Listen for work
    zmq::socket_t backend(context_, ZMQ_REQ);
    backend.connect("inproc://worker");
    // Listen for control messages
    zmq::socket_t control(context_, ZMQ_SUB);
    control.connect("inproc://control");
    control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    // Signal main thread that we're alive
    zmq::socket_t heartbeat(context_, ZMQ_REQ);
    heartbeat.connect("inproc://heartbeat");
    {
        zmq::message_t ping;
        heartbeat.send(ping);
        heartbeat.recv(&ping);
    }
    
    // Register with server thread
    backend.send(zmq::message_t());
    
    zmq::pollitem_t items[] = {
      { static_cast<void *>(backend), 0, ZMQ_POLLIN, 0 },
      { static_cast<void *>(control), 0, ZMQ_POLLIN, 0 },
    };
    
    while (true) {
        
        if (int rc = zmq::poll(&items[0], 2, -1) < 0) {
            log_error_stream("ZMQ polling error " << rc);
            break;
        }
        
        // Message from client
        if (items[0].revents & ZMQ_POLLIN) {
            {
                zmq::message_t steps_msg, client_msg;
                backend.recv(&steps_msg);
                i3_assert( steps_msg.more() );
                backend.recv(&client_msg);
                i3_assert( !client_msg.more() );
                uint32_t client_id = deserialize<uint32_t>(client_msg);
                
                // Submit job
                auto steps = deserialize<I3CLSimStepSeriesPtr>(steps_msg);
                i3_assert( steps );
                // boost::posix_time::ptime now(boost::posix_time::microsec_clock::local_time());
                // log_info("[%s %u,%u] enqueuing %zu steps with id [%u]", to_simple_string(now).c_str(), index, buffer_id, steps->size(), client_id);
                converters_[index]->EnqueueSteps(steps, client_id);
                // log_info("[%s %u,%u] enqueued id [%u]", to_simple_string(now).c_str(), index, buffer_id, client_id);
                
            }
            
            // Get next result (not necessarily from the batch we just enqueued)
            I3CLSimStepToPhotonConverter::ConversionResult_t result =
                converters_[index]->GetConversionResult();
            
            // boost::posix_time::ptime now(boost::posix_time::microsec_clock::local_time());
            // log_info("[%s %u,%u] got photons with id [%u], queue size %zu", to_simple_string(now).c_str(), index, buffer_id,
            //     result.identifier, converters_[index]->QueueSize());
            
            // Send result back to server thread
            if (result.photons)
                backend.send(serialize(result.photons), ZMQ_SNDMORE);
            if (result.photonHistories)
                backend.send(serialize(result.photonHistories), ZMQ_SNDMORE);
            backend.send(serialize(result.identifier), 0);
        }
        
        // Shutdown signal from main thread
        if (items[1].revents & ZMQ_POLLIN) {
            zmq::message_t dummy;
            control.recv(&dummy);
            break;
        }
        
    }
}

std::map<std::string, double> I3CLSimServer::GetStatistics() const
{
    std::map<std::string, double> summary;
    
    for (std::size_t i=0; i<converters_.size(); ++i)
    {
        const std::string postfix = (converters_.size()==1)?"":"_"+boost::lexical_cast<std::string>(i);
        for (auto &v : converters_[i]->GetStatistics()) {
            summary[v.first + postfix] = v.second;
        }
    }
    
    return summary;
}

I3CLSimClient::I3CLSimClient(const std::string &server_address) : context_(1),
    socket_(context_, ZMQ_DEALER), workgroupSize_(0), maxBunchSize_(0), pending_(0)
{
    socket_.connect(server_address);
    
    zmq::message_t reply;
    std::string greeting = "servus";
    socket_.send(zmq::message_t(greeting.data(), greeting.size()), 0);
    socket_.recv(&reply);
    std::pair<std::size_t,std::size_t> bunchSize = deserialize<std::pair<std::size_t,std::size_t> >(reply);
    workgroupSize_ = bunchSize.first;
    maxBunchSize_ = bunchSize.second;
}

void I3CLSimClient::EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier)
{
    socket_.send(serialize(steps), ZMQ_SNDMORE);
    socket_.send(serialize(identifier), 0);
    pending_++;
}

I3CLSimStepToPhotonConverter::ConversionResult_t
I3CLSimClient::GetConversionResult()
{
    I3CLSimStepToPhotonConverter::ConversionResult_t result;
    if (pending_ != 0) {
        
        zmq::message_t address;
        std::list<zmq::message_t> body;
        std::tie(address, body) = read_message(socket_);
        
        // End of message packet is always the identifier
        i3_assert( body.size() > 0 );
        result.identifier = deserialize<decltype(result.identifier)>(body.back());
        body.pop_back();
        
        // Optional: photons
        if (body.size() > 0) {
            result.photons = deserialize<decltype(result.photons)>(body.front());
            body.pop_front();
        }
        if (body.size() > 0) {
            result.photonHistories = deserialize<decltype(result.photonHistories)>(body.front());
            body.pop_front();
        }
        i3_assert( body.size() == 0 );
            
        pending_--;
    }
    
    return result;
}
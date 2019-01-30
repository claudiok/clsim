
#ifndef CLSIM_I3CLSIMSERVER_H_INCLUDED
#define CLSIM_I3CLSIMSERVER_H_INCLUDED

#include "icetray/I3Logging.h"
#include "clsim/I3CLSimStep.h"
#include "clsim/I3CLSimStepToPhotonConverter.h"

#include <thread>
#include <memory>
#include <map>
#include <list>
#include <queue>
#include <zmq.hpp>

class I3CLSimServer {
public:
    I3CLSimServer(const std::string &address, const std::vector<I3CLSimStepToPhotonConverterPtr> &converters);
    ~I3CLSimServer();
    
    std::map<std::string, double> GetStatistics() const;
private:
    std::vector<I3CLSimStepToPhotonConverterPtr> converters_;
    std::size_t workgroupSize_, maxBunchSize_;
    
    zmq::context_t context_;
    zmq::socket_t frontend_, backend_, control_, heartbeat_;
    
    SET_LOGGER("I3CLSimServer");
    
    void ServerThread();
    void WorkerThread(unsigned index, unsigned buffer_id);
    
    std::thread serverThread_;
    std::vector<std::thread> workerThreads_;
    
    /// Addresses of idle workers
    std::queue<zmq::message_t > workers_;
    
    /// Map from internal task ID to address of originating client and external ID
    std::map<uint32_t, std::array<zmq::message_t, 2> > clients_;
};

class I3CLSimClient {
public:
    I3CLSimClient(const std::string &server_address);
    
    void EnqueueSteps(I3CLSimStepSeriesConstPtr steps, uint32_t identifier);
    I3CLSimStepToPhotonConverter::ConversionResult_t GetConversionResult();
    
    std::size_t GetWorkgroupSize() const { return workgroupSize_; }
    std::size_t GetMaxNumWorkitems() const { return maxBunchSize_; }
    
private:
    zmq::context_t context_;
    zmq::socket_t socket_;
    std::size_t workgroupSize_, maxBunchSize_;

    SET_LOGGER("I3CLSimClient");

    uint32_t pending_;
};

#endif // CLSIM_I3CLSIMSERVER_H_INCLUDED



#include "clsim/I3CLSimLightSourcePropagatorFromI3PropagatorService.h"
#include "clsim/I3CLSimLightSource.h"

I3CLSimLightSourcePropagatorFromI3PropagatorService::I3CLSimLightSourcePropagatorFromI3PropagatorService(I3ParticleTypePropagatorServiceMapPtr particleToPropagatorServiceMap)
    : particleToPropagatorServiceMap_(particleToPropagatorServiceMap), initialized_(false)
{
    if (!particleToPropagatorServiceMap_)
        log_fatal("(null) propagator map is invalid");
    for (auto &pair : *particleToPropagatorServiceMap_) {
        if (!pair.second)
            log_fatal("(null) propagator in map is invalid");
        I3Particle p;
        p.SetType(pair.first);
        log_info("Propagating Particle %s with %s",p.GetTypeString().c_str(),
                 typeid(*(pair.second)).name());
    }
}

I3CLSimLightSourcePropagatorFromI3PropagatorService::~I3CLSimLightSourcePropagatorFromI3PropagatorService()
{}

void I3CLSimLightSourcePropagatorFromI3PropagatorService::SetRandomService(I3RandomServicePtr rng)
{
  for (auto &pair : *particleToPropagatorServiceMap_){
    pair.second->SetRandomNumberGenerator(rng);
  }
}

bool I3CLSimLightSourcePropagatorFromI3PropagatorService::IsValidForLightSource(const I3CLSimLightSource &source)
{
    return (source.GetType() == I3CLSimLightSource::Particle)
        && (particleToPropagatorServiceMap_->find(source.GetParticle().GetType()) != particleToPropagatorServiceMap_->end());
}

void I3CLSimLightSourcePropagatorFromI3PropagatorService::Convert(I3CLSimLightSourceConstPtr &lightSource, uint32_t identifier,
    secondary_callback emitSecondary, step_callback)
{
    std::deque<std::pair<I3Particle, I3PropagatorServicePtr>> queue;
    
    auto iter = particleToPropagatorServiceMap_->find(lightSource->GetParticle().GetType());
    assert( iter != particleToPropagatorServiceMap_->end());
    queue.emplace_back(lightSource->GetParticle(), iter->second);

    auto emitNonDark = [&](const I3Particle &particle) {
        if (particle.GetShape() != I3Particle::Dark) {
            I3CLSimLightSourceConstPtr secondary = boost::make_shared<I3CLSimLightSource>(particle);
            emitSecondary(secondary, identifier);
        }
    };
    
    while (!queue.empty()) {
        // retrieve the first entry
        auto &item = queue.front();

        // propagate it!
        for (auto &particle : item.second->Propagate(item.first, NULL, NULL)) {
            
            auto iter = particleToPropagatorServiceMap_->find(particle.GetType());
            if (iter->second != item.second && iter != particleToPropagatorServiceMap_->end()) {
                queue.emplace_back(particle, iter->second);
            } else {
                emitNonDark(particle);
            }
        }
        // emit parent if the propagator didn't set its shape to dark
        emitNonDark(item.first);

        queue.pop_front();
    }
}


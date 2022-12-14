
#include <mixers/Mixer.hpp>

using namespace sisi4s;

template <typename F>
Mixer<F>::Mixer(Algorithm *algorithm_)
    : algorithm(algorithm_) {}

template <typename F>
Mixer<F>::~Mixer() {}

// instantiate
template class sisi4s::Mixer<sisi4s::Float64>;
template class sisi4s::Mixer<sisi4s::Complex64>;

template <typename F>
std::map<std::string, std::function<PTR(Mixer<F>)(Algorithm *algorithm)>>
    *MixerFactory<F>::mixerMap;

// instantiate
template class sisi4s::MixerFactory<sisi4s::Float64>;
template class sisi4s::MixerFactory<sisi4s::Complex64>;

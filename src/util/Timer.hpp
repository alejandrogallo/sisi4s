#ifndef TIMER_DEFINED
#define TIMER_DEFINED

#include <util/Time.hpp>

namespace sisi4s {
/**
 * Timer class providing timing functionality.
 * The object can be created in a scope whose lifetime
 * is then measured.
 */
class Timer {
public:
  Timer(Time *time);
  ~Timer();

protected:
  Time *time;
  Time start;
};
} // namespace sisi4s

#endif

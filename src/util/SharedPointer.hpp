#ifndef SHARED_POINTER_DEFINED
#define SHARED_POINTER_DEFINED

#include <memory>

// provide convenience macros for shared pointers

#define PTR(TYPE) std::shared_ptr<TYPE>
#define WEAK_PTR(TYPE) std::weak_ptr<TYPE>
#define NEW(TYPE, ...) std::make_shared<TYPE>(__VA_ARGS__)

#endif

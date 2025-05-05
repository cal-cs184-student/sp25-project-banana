#ifndef CGL_LINUX_TYPES_H
#define CGL_LINUX_TYPES_H

// This header ensures Linux kernel type definitions are available
// when they're needed by system headers

#ifdef __linux__
// Include Linux kernel types
#include <linux/types.h>

// If for some reason the types aren't defined even with the include,
// define them here as fallbacks
#ifndef __u32
typedef unsigned int __u32;
#endif

#ifndef __s32
typedef signed int __s32;
#endif

#ifndef __u64
typedef unsigned long long __u64;
#endif
#endif // __linux__

#endif // CGL_LINUX_TYPES_H

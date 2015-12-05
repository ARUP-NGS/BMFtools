#ifndef COMPILER_UTIL_H
#define COMPILER_UTIL_H

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#ifdef __GNUC__
#	ifndef UNUSED
#		define UNUSED(x) x __attribute__((unused))
#	endif
#	ifndef UNUSED_FUNC
#		define UNUSED_FUNC(x) __attribute__((__unused__)) x
#	endif
#else
#	define UNUSED(x) x
#	define UNUSED_FUNC(x) x
#endif

#endif

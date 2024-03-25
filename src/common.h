#ifndef PHASE2_COMMON_H
#define PHASE2_COMMON_H

#include <stdint.h>

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"

typedef int8_t	i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t	 u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float  f32;
typedef double f64;

typedef f64 fl;

enum eno {
	OK = 0,
	EMPI,
	EMEM,
	ESIZE,
};

#endif // PHASE2_COMMON_H

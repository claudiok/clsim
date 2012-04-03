#ifndef IEEE_HALF_PRECISION_H_
#define IEEE_HALF_PRECISION_H_

// Conversion betwen doubles/floats and 16-bit floating point.

#ifdef __cplusplus
extern "C" {
#endif

int singles2halfp(void *target, const void *source, int numel);
int doubles2halfp(void *target, const void *source, int numel);
int halfp2singles(void *target, const void *source, int numel);
int halfp2doubles(void *target, const void *source, int numel);

#ifdef __cplusplus
}
#endif

#endif

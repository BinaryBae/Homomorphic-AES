#ifndef _KEYSWITCH_H
#define _KEYSWITCH_H

#include <tfhe.h>

struct BaseBKeySwitchKey
{
	int32_t n;
	int32_t t;
	int32_t basebit;
	int32_t b;
	int32_t base;
	const TLweParams *out_params;
	TLweSample *ks0_raw;
	TLweSample ***ks;

	BaseBKeySwitchKey(int32_t n, int32_t t, int32_t basebit, int32_t b, const TLweParams *out_params, TLweSample *ks0_raw);
	~BaseBKeySwitchKey();
	BaseBKeySwitchKey(const BaseBKeySwitchKey &) = delete;
	void operator=(const BaseBKeySwitchKey &) = delete;
};

BaseBKeySwitchKey *alloc_BaseBKeySwitchKey();

void free_BaseBKeySwitchKey(BaseBKeySwitchKey *ptr);

void init_BaseBKeySwitchKey(BaseBKeySwitchKey *obj, int n, int t, int basebit, int b, const TLweParams *out_params);

void destroy_BaseBKeySwitchKey(BaseBKeySwitchKey *obj);

BaseBKeySwitchKey *new_BaseBKeySwitchKey(int n, int t, int basebit, int b, const TLweParams *out_params);

void delete_BaseBKeySwitchKey(BaseBKeySwitchKey *obj);

struct BaseBExtra
{
public:
	static void CreateKeySwitchKey_fromArray(
		TLweSample ***result,
		const TLweKey *out_key,
		const double out_alpha,
		const int *in_key,
		const int n,
		const int t,
		const int basebit,
		const int mod);

	static void KeySwitchTranslate_fromArray(
		TLweSample *result,
		const TLweSample ***ks,
		const TLweParams *params,
		const Torus32 *ai,
		const int n,
		const int t,
		const int basebit);

	static void KeySwitchTranslate_fromArray_Generic(
		TLweSample *result,
		const TLweSample ***ks,
		const TLweParams *params,
		const Torus32 **as,
		const int n,
		const int t,
		const int basebit,
		const int M);

	static void CreateKeySwitchKey(
		BaseBKeySwitchKey *result,
		const LweKey *in_key,
		const TLweKey *out_key);

	static void KeySwitch(
		TLweSample *result,
		const BaseBKeySwitchKey *ks,
		const LweSample *sample);

	static void KeySwitch_Id(
		TLweSample *result,
		const BaseBKeySwitchKey *ks,
		std::vector<LweSample *> samples);
};

#endif

#include "homo_aes.h"
#include "tree_mvb.h"
#include <time.h>

#define NUM_THREADS 8

void AddRoundKey_fhe(vector<LweSample *> a[4][MAXBC],
					 vector<LweSample *> rk[4][MAXBC],
					 TFheGateBootstrappingSecretKeySet *gk,
					 BaseBKeySwitchKey *ks_key)
{
	int i, j, k;
#pragma omp parallel for default(shared) private(k, i, j) num_threads(NUM_THREADS)
	for (k = 0; k < 16; ++k)
	{
		i = k >> 2;
		j = k % 4;
		XOR_fhe(a[i][j], rk[i][j], gk, ks_key);
	}
}

void AddRoundKey_nextSubBytes_fhe(vector<LweSample *> a[4][MAXBC],
								  vector<LweSample *> rk[4][MAXBC],
								  TFheGateBootstrappingSecretKeySet *gk,
								  BaseBKeySwitchKey *ks_key)
{
	int i, j, k;

#pragma omp parallel for default(shared) private(k, i, j) num_threads(NUM_THREADS)
	for (k = 0; k < 16; ++k)
	{
		vector<LweSample *> tmp;
		tmp.push_back(new_LweSample(gk->lwe_key->params));
		tmp.push_back(new_LweSample(gk->lwe_key->params));
		tmp.push_back(new_LweSample(gk->lwe_key->params));
		i = k >> 2;
		j = k % 4;
		XOR_fhe(a[i][j], rk[i][j], gk, ks_key);
		lweCopy(tmp[0], a[i][j][0], gk->lwe_key->params);
		lweCopy(tmp[1], a[i][j][1], gk->lwe_key->params);
		lweCopy(tmp[2], a[i][j][1], gk->lwe_key->params);
		deref_boot_opti(a[i][j], gk, tmp, ks_key, Sbox);
	}
}

void ShiftRows_fhe(vector<LweSample *> a[4][MAXBC],
				   word8 d,
				   int BC,
				   TFheGateBootstrappingSecretKeySet *gk)
{
	vector<LweSample *> tmp_shift[MAXBC];
	for (int w = 0; w < MAXBC; ++w)
	{
		tmp_shift[w].push_back(new_LweSample(gk->lwe_key->params));
		tmp_shift[w].push_back(new_LweSample(gk->lwe_key->params));
	}
	int i, j;
	if (d == 0)
	{
		for (i = 1; i < 4; i++)
		{
			for (j = 0; j < BC; j++)
			{
				lweCopy(tmp_shift[j][0], a[i][(j + shifts[0][i]) % BC][0], gk->lwe_key->params);
				lweCopy(tmp_shift[j][1], a[i][(j + shifts[0][i]) % BC][1], gk->lwe_key->params);
			}
			for (j = 0; j < BC; j++)
			{
				a[i][j].swap(tmp_shift[j]);
			}
		}
	}
	else
	{
		for (i = 1; i < 4; i++)
		{
			for (j = 0; j < BC; j++)
			{
				lweCopy(tmp_shift[j][0], a[i][(BC + j - shifts[0][i]) % BC][0], gk->lwe_key->params);
				lweCopy(tmp_shift[j][1], a[i][(BC + j - shifts[0][i]) % BC][1], gk->lwe_key->params);
			}
			for (j = 0; j < BC; j++)
			{
				a[i][j].swap(tmp_shift[j]);
			}
		}
	}
}

void MixColumns_fhe(vector<LweSample *> a[4][MAXBC],
					TFheGateBootstrappingSecretKeySet *gk,
					BaseBKeySwitchKey *ks_key)
{
	vector<LweSample *> b[4][MAXBC];
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			b[i][j].push_back(new_LweSample(gk->lwe_key->params));
			b[i][j].push_back(new_LweSample(gk->lwe_key->params));
		}
	}

	int i, j, k;
#pragma omp parallel for default(shared) private(k, i, j) num_threads(NUM_THREADS)
	for (k = 0; k < 16; ++k)
	{
		j = k >> 2;
		i = k % 4;
		vector<LweSample *> tmp(2);
		tmp[0] = new_LweSample(gk->lwe_key->params);
		tmp[1] = new_LweSample(gk->lwe_key->params);

		mul2_fhe(b[i][j], a[i][j], gk, ks_key);
		mul3_fhe(tmp, a[(i + 1) % 4][j], gk, ks_key);
		XOR_fhe(b[i][j], tmp, gk, ks_key);
		XOR_fhe(b[i][j], a[(i + 2) % 4][j], gk, ks_key);
		XOR_fhe(b[i][j], a[(i + 3) % 4][j], gk, ks_key);
		tmp.clear();
	}
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			a[i][j].swap(b[i][j]);
		}
}

int KeyExpansion(word8 k[4][MAXKC],
				 word8 W[MAXROUNDS + 1][4][MAXBC],
				 int BC,
				 int KC,
				 int ROUNDS)
{
	int i, j, t, RCpointer = 1;
	word8 tk[4][MAXKC];
	for (j = 0; j < KC; j++)
		for (i = 0; i < 4; i++)
			tk[i][j] = k[i][j];
	t = 0;
	for (j = 0; (j < KC) && (t < (ROUNDS + 1) * BC); j++, t++)
		for (i = 0; i < 4; i++)
			W[t / BC][i][t % BC] = tk[i][j];

	while (t < (ROUNDS + 1) * BC)
	{
		tk[0][0] ^= S[tk[1][3]];
		tk[1][0] ^= S[tk[2][3]];
		tk[2][0] ^= S[tk[3][3]];
		tk[3][0] ^= S[tk[0][3]];

		tk[0][0] ^= RC[RCpointer++];
		if (KC <= 6)
			for (j = 1; j < KC; j++)
				for (i = 0; i < 4; i++)
					tk[i][j] ^= tk[i][j - 1];
		else
		{
			for (j = 1; j < 4; j++)
				for (i = 0; i < 4; i++)
					tk[i][j] ^= tk[i][j - 1];
			for (i = 0; i < 4; i++)
				tk[i][4] ^= S[tk[i][3]];
			for (j = 5; j < KC; j++)
				for (i = 0; i < 4; i++)
					tk[i][j] ^= tk[i][j - 1];
		}
		for (j = 0; (j < KC) && (t < (ROUNDS + 1) * BC); j++, t++)
			for (i = 0; i < 4; i++)
				W[t / BC][i][t % BC] = tk[i][j];
	}
	return 0;
}

int Encrypt_fhe(vector<LweSample *> a[4][MAXBC],
				vector<LweSample *> rk[MAXROUNDS + 1][4][MAXBC],
				int ROUNDS,
				TFheGateBootstrappingSecretKeySet *gk,
				BaseBKeySwitchKey *ks_key)
{
	printf("NUM_THREADS = %d\n", NUM_THREADS);
	int r;
	AddRoundKey_nextSubBytes_fhe(a, rk[0], gk, ks_key);
	for (r = 1; r < ROUNDS; r++)
	{
		ShiftRows_fhe(a, 0, 4, gk);
		MixColumns_fhe(a, gk, ks_key);
		AddRoundKey_nextSubBytes_fhe(a, rk[r], gk, ks_key);
	}
	ShiftRows_fhe(a, 0, 4, gk);
	AddRoundKey_fhe(a, rk[ROUNDS], gk, ks_key);
	return 0;
}

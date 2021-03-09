#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

#pragma comment(linker, "/stack:200000000")
#pragma GCC optimize("Ofast,unroll-loops,no-stack-protector,fast-math")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

#include "q2_trace.h"

using ll = long long;
#define mp make_pair
#define eb emplace_back
#define all(v) (v).begin(), (v).end()

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
#define shuf(v) shuffle((v).begin(), (v).end(), rng)

template<typename T>
struct Fraction {
	T x;
	int e;
	Fraction(T a, int b) {
		x = a, e = b;
		simplify();
	}
	Fraction(const T& i) {
		x = i, e = 0;
	}
	void simplify() {
		while (x && !(x & 1))
			x >>= 1, e--;
	}
	void update(const Fraction& oth) {
		if (x == 0 || oth.x == 0) 
			x = 0, e = 0;
		else
			x *= oth.x, e += oth.e - 1;
	}
	bool operator<(const Fraction& oth) {
		T tmp = T(1) << abs(e - oth.e);
		if (e < oth.e)
			return abs(tmp * x) < abs(oth.x);
		else
			return abs(x) < abs(tmp * oth.x);
	}
	bool operator==(Fraction& oth) {
		return (abs(x) == abs(oth.x) && e == oth.e);
	}
	friend ostream& operator<<(ostream& os, const Fraction& a) {
		if (a.e <= 0) return os << (a.x * (T(1) << (-a.e)));
		return os << a.x << '/' << (T(1) << a.e);
	}
};
using frac = Fraction<ll>;

int NUM_BITS, NUM_STAGES;
int UNITS, SBOX_BITS, SBOX_SZ;
vector<int> SBOX, PBOX;
vector<int> PINV_MEM;

vector<vector<frac>> bias, dp;
vector<vector<int>> parent;

string binary(int p) {
	string s;
	for (int i = SBOX_BITS - 1; i >= 0; i--)
		if ((p >> i) & 1)
			s += '1';
		else
			s += '0';
	return s;
}
int permute(int plain, vector<int>& p) {
	int ret = 0, n = p.size();
	for (int j = 0; j < n; j++) {
		if (!((plain >> j) & 1)) continue;
		int c = p[NUM_BITS - j - 1];
		c = NUM_BITS - c - 1;
		ret ^= (1 << c);
	}
	return ret;
}
void print(int plain) {
	vector<int> sboxs(UNITS);
	for (int i = UNITS - 1; i >= 0; i--) {
		sboxs[i] = plain & (SBOX_SZ - 1);
		plain >>= SBOX_BITS;
	}
	for (int j = 0; j < UNITS; j++) {
		cout << binary(sboxs[j]);
		if (j < UNITS - 1) cout << '|';
	}
	cout << endl;
}
void draw(int plain) {
	string sep(NUM_BITS + UNITS - 1, '*');
	cout << sep << endl;
	string dash(SBOX_BITS, '-');
	for (int i = 0; i < NUM_STAGES - 1; i++) {
		print(plain);
		for (int j = 0; j < UNITS; j++) {
			cout << dash;
			if (j < UNITS - 1) cout << '+';
		}
		cout << endl;
		plain = parent[i][plain];
		print(PINV_MEM[plain]);
		cout << endl;
	}
	print(plain);
	cout << sep << endl;
}
void random_init() {
	NUM_BITS = 16, NUM_STAGES = 5, SBOX_BITS = 4;
	assert(NUM_BITS % SBOX_BITS == 0);
	UNITS = NUM_BITS / SBOX_BITS;
	SBOX_SZ = 1 << SBOX_BITS;
	SBOX.resize(SBOX_SZ);
	PBOX.resize(NUM_BITS);
	iota(all(SBOX), 0), iota(all(PBOX), 0);
	shuf(SBOX), shuf(PBOX);
}
void custom_init() {
	SBOX = {14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7};
	PBOX = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};
	// SBOX = {8, 15, 2, 10, 12, 11, 6, 13, 14, 5, 4, 0, 1, 3, 9, 7};
	// PBOX = {9, 8, 2, 6, 1, 7, 11, 3, 10, 4, 0, 5};
	NUM_STAGES = 5;
	NUM_BITS = PBOX.size();
	SBOX_SZ = SBOX.size();
	SBOX_BITS = int(log2(SBOX_SZ));
	UNITS = NUM_BITS / SBOX_BITS;
	assert((1 << SBOX_BITS) == SBOX_SZ);
	assert(NUM_BITS % SBOX_BITS == 0);
}
int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	cout.tie(NULL);

	// random_init();
	custom_init();

	trace(NUM_BITS, NUM_STAGES, SBOX_BITS);
	trace(SBOX);
	trace(PBOX);

	// initialize inverse boxes
	PINV_MEM.resize(1 << NUM_BITS);
	for (int i = 0; i < (1 << NUM_BITS); i++)
		PINV_MEM[permute(i, PBOX)] = i;

	// calculate biases
	bias.resize(SBOX_SZ);
	bias[0].assign(SBOX_SZ, 0);
	for (int i = 1; i < SBOX_SZ; i++) {
		bias[i].assign(SBOX_SZ, 0);
		for (int j = 1; j < SBOX_SZ; j++) {
			int cnt = 0, cur;
			for (int k = 0; k < SBOX_SZ; k++) {
				cur = __builtin_popcount(SBOX[k] & j);
				cur += __builtin_popcount(k & i);
				if (!(cur & 1)) cnt++;
			}
			assert(!(cnt & 1));
			cnt = (cnt << 1) - SBOX_SZ;
			if (cnt != 0)
				bias[i][j] = frac(cnt, SBOX_BITS + 1);
		}
	}
	
	// initialise dp and parent
	parent.resize(NUM_STAGES), dp.resize(NUM_STAGES);
	for (int i = 0; i < NUM_STAGES - 1; i++) {
		dp[i].assign(1 << NUM_BITS, 0);
		parent[i].assign(1 << NUM_BITS, -1);
	}
	dp[NUM_STAGES - 1].assign(1 << NUM_BITS, frac(1, 1));
	iota(all(parent[NUM_STAGES - 1]), 0);

	// declare auxiliary variables outside BIG for loop
	frac max_bias(0), current_total_bias(0);
	vector<int> nxt_sboxes_mask(UNITS);

	int rev_mask, bits, stage;
	int sub_mask, cur_mask, bad;
	int inp_mask, out_mask;
	int nxt_mask, mask, ptr;

	// calculate dp and parents
	for (int s = NUM_STAGES - 2; s >= 0; s--) {
		max_bias = frac(0);
		for (nxt_mask = 1; nxt_mask < (1 << NUM_BITS); nxt_mask++) {
			if (dp[s + 1][nxt_mask].x == 0) continue;
			rev_mask = PINV_MEM[nxt_mask], bits = 0;
			for (ptr = UNITS - 1; ptr >= 0; ptr--) {
				nxt_sboxes_mask[ptr] = rev_mask & ((1 << SBOX_BITS) - 1);
				rev_mask >>= SBOX_BITS;
				if (nxt_sboxes_mask[ptr]) bits += SBOX_BITS;
			}
			for (mask = 1; mask < (1 << bits); mask++) {
				sub_mask = mask, cur_mask = 0, bad = 0;
				current_total_bias = dp[s + 1][nxt_mask];
				for (ptr = UNITS - 1; ptr >= 0; ptr--) {
					if (nxt_sboxes_mask[ptr] == 0) continue;
					inp_mask = sub_mask & ((1 << SBOX_BITS) - 1);
					out_mask = nxt_sboxes_mask[ptr];
					sub_mask >>= SBOX_BITS;
					cur_mask += inp_mask << ((UNITS - ptr - 1) * SBOX_BITS);
					if (bias[inp_mask][out_mask].x == 0) bad = 1;
					current_total_bias.update(bias[inp_mask][out_mask]);
				}
				if (!bad && dp[s][cur_mask] < current_total_bias) {
					dp[s][cur_mask] = current_total_bias;
					parent[s][cur_mask] = nxt_mask;
					if (max_bias < current_total_bias)
						max_bias = current_total_bias;
				}
			}
		}
		stage = NUM_STAGES - s;
		if (max_bias.x == 0) {
			trace(stage, max_bias);
			return 0;
		}
		trace(stage, max_bias);
	}

	vector<int> final_parent(1 << NUM_BITS, -1);
	int min_cnt = UNITS + 1, plaintext;
	for (int i = 0; i < (1 << NUM_BITS); i++) {
		int j = i, cnt = 0;
		for (int k = 0; k < NUM_STAGES - 1; k++)
			j = parent[k][j];
		final_parent[i] = j;
		if (!(max_bias == dp[0][i])) continue;
		for (int k = UNITS - 1; k >= 0; k--) {
			if (j & (SBOX_SZ - 1)) cnt++;
			j >>= SBOX_BITS;
		}
		if (cnt < min_cnt)
			plaintext = i, min_cnt = cnt;
	}

	draw(plaintext);

	cout << "Active SBoxes : ";
	cout << min_cnt << endl;

	cout << "Plaintext bits :";
	for (int i = 0; i < NUM_BITS; i++)
		if ((plaintext >> i) & 1) cout << " " << i;
	cout << endl;

	cout << "Ciphertext bits :";
	for (int i = 0; i < NUM_BITS; i++)
		if ((final_parent[plaintext] >> i) & 1)
			cout << " " << i;
	cout << endl;
}

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#pragma comment(linker, "/stack:200000000")
#pragma GCC optimize("Ofast,unroll-loops,no-stack-protector,fast-math")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

using namespace std;
#ifdef LOCAL
#include "q2_trace.h"
#else
#define trace(args...)
#endif

using ll = long long;
using ld = long double;
using pii = pair<int, int>;
using vi = vector<int>;
#define mp make_pair
#define ub upper_bound
#define lb lower_bound
#define fi first
#define se second
#define pb push_back
#define eb emplace_back
#define all(v) (v).begin(), (v).end()
#define rall(v) (v).rbegin(), (v).rend()

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
#define shuf(v) shuffle((v).begin(), (v).end(), rng)

struct frac {
	ll x;
	int e;
	frac(int a, int b) {
		assert(b >= 0);
		x = a, e = b;
		simplify();
	}
	frac(const int& i) {
		x = i, e = 0;
	}
	void simplify() {
		while (x && !(x & 1))
			x >>= 1, e--;
	}
	void update(const frac& oth) {
		x *= oth.x; e += oth.e - 1;
	}
	bool operator<=(frac& oth) {
		ll tmp = 1LL << abs(e - oth.e);
		if (e < oth.e)
			return abs(tmp * x) <= abs(oth.x);
		else
			return abs(x) <= abs(tmp * oth.x);
	}
	bool operator==(frac& oth) {
		return (x == oth.x && e == oth.e);
	}
	friend ostream& operator<<(ostream& os, frac const& a) {
		if (a.e == 0) return os << a.x;
		if (a.e < 0) return os << (a.x * (1LL << (-a.e)));
		return os << a.x << '/' << (1LL << a.e);
	}
};

int NUM_BITS, NUM_STAGES;
int UNITS, SBOX_SZ;
vector<int> SBOX, PBOX;
vector<int> SINV, PINV;

vector<vector<frac>> bias, dp;
vector<vector<int>> parent;

string binary(int p) {
	string s;
	for (int i = SBOX_SZ - 1; i >= 0; i--)
		if ((p >> i) & 1)
			s += '1';
		else
			s += '0';
	return s;
}
int permute(int plain, vector<int>& p) {
	int ret = 0, n = p.size();
	for (int j = 0; j < n; j++)
		if ((plain >> j) & 1)
			ret ^= (1 << p[j]);
	return ret;
}
void print(int plain) {
	vector<int> sboxs(UNITS);
	for (int i = UNITS - 1; i >= 0; i--) {
		sboxs[i] = plain % (1 << SBOX_SZ);
		plain >>= SBOX_SZ;
	}
	for (int j = 0; j < UNITS; j++) {
		cerr << binary(sboxs[j]);
		if (j < UNITS - 1) cerr << '|';
	}
	cerr << endl;
}
void draw(int plain) {
	string sep(NUM_BITS + UNITS - 1, '*');
	cerr << sep << endl;
	string dash(SBOX_SZ, '-');
	for (int i = 0; i < NUM_STAGES - 1; i++) {
		print(plain);
		for (int j = 0; j < UNITS; j++) {
			cerr << dash;
			if (j < UNITS - 1) cerr << '+';
		}
		cerr << endl;
		plain = permute(parent[i][plain], PINV);
		print(plain);
		plain = permute(plain, PBOX);
		cerr << endl;
	}
	print(plain);
	cerr << sep << endl;
}
void random_init() {
	SBOX_SZ = 4, NUM_BITS = 16;
	assert(NUM_BITS % SBOX_SZ == 0);
	UNITS = NUM_BITS / SBOX_SZ;
	NUM_STAGES = 5;

	SBOX.resize(1 << SBOX_SZ);
	PBOX.resize(NUM_BITS);

	iota(all(SBOX), 0), iota(all(PBOX), 0);
	shuf(SBOX), shuf(PBOX);
}
void custom_init() {
	SBOX = {8, 15, 2, 10, 12, 11, 6, 13, 14, 5, 4, 0, 1, 3, 9, 7};
	PBOX = {9, 8, 2, 6, 1, 7, 11, 3, 10, 4, 0, 5};
	NUM_STAGES = 7;
	SBOX_SZ = int(log2(SBOX.size()));
	NUM_BITS = PBOX.size();
	UNITS = NUM_BITS / SBOX_SZ;
	assert((1 << SBOX_SZ) == SBOX.size());
	assert(NUM_BITS % SBOX_SZ == 0);
}
int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	cout.tie(NULL);

	// random_init();
	custom_init();

	// initialize inverse boxes
	SINV.resize(1 << SBOX_SZ);
	PINV.resize(NUM_BITS);
	for (int i = 0; i < (1 << SBOX_SZ); i++)
		SINV[SBOX[i]] = i;
	for (int i = 0; i < NUM_BITS; i++)
		PINV[PBOX[i]] = i;

	// calculate biases
	bias.resize(1 << SBOX_SZ);
	bias[0].assign(1 << SBOX_SZ, frac(0));
	for (int i = 1; i < (1 << SBOX_SZ); i++) {
		bias[i].assign(1 << SBOX_SZ, frac(0, SBOX_SZ));
		for (int j = 1; j < (1 << SBOX_SZ); j++) {
			frac& b = bias[i][j];
			for (int k = 0; k < (1 << SBOX_SZ); k++) {
				int cur = __builtin_popcount((k & i) ^ (SBOX[k] & j));
				if (cur % 2 == 0) b.x++;
			}
			assert(b.x % 2 == 0);
			b.x = 2 * b.x - (1LL << b.e); b.e++;
			if (b.x == 0) b.e = 0;
			b.simplify();
		}
	}

	// trace(bias[0b1111][0b1100]);
	// trace(bias[0b1100][0b0100]);
	// trace(bias[0b1000][0b1000]);
	// trace(bias[0b0100][0b0001]);
	// trace(bias[0b0010][0b1001]);
	// trace(bias[0b0101][0b1111]);

	trace(NUM_BITS, NUM_STAGES, SBOX_SZ);
	trace(SBOX);
	trace(PBOX);

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
		for (nxt_mask = 1; nxt_mask < (1 << NUM_BITS); nxt_mask++) {
			if (dp[s + 1][nxt_mask].x == 0) continue;
			rev_mask = permute(nxt_mask, PINV), bits = 0;
			for (ptr = UNITS - 1; ptr >= 0; ptr--) {
				nxt_sboxes_mask[ptr] = rev_mask & ((1 << SBOX_SZ) - 1);
				rev_mask >>= SBOX_SZ;
				if (nxt_sboxes_mask[ptr]) bits += SBOX_SZ;
			}
			for (mask = 1; mask < (1 << bits); mask++) {
				sub_mask = mask, cur_mask = 0, bad = 0;
				current_total_bias = dp[s + 1][nxt_mask];
				for (ptr = UNITS - 1; ptr >= 0; ptr--) {
					if (nxt_sboxes_mask[ptr] == 0) continue;
					inp_mask = sub_mask & ((1 << SBOX_SZ) - 1);
					out_mask = nxt_sboxes_mask[ptr];
					sub_mask >>= SBOX_SZ;
					cur_mask += inp_mask << ((UNITS - ptr - 1) * SBOX_SZ);
					if (bias[inp_mask][out_mask].x == 0) bad = 1;
					current_total_bias.update(bias[inp_mask][out_mask]);
				}
				if (!bad && dp[s][cur_mask] <= current_total_bias) {
					if (s == 0 && mask == 15) {
						rev_mask = permute(nxt_mask, PINV);
						print(cur_mask); print(rev_mask);
						trace(current_total_bias);
					}
					dp[s][cur_mask] = current_total_bias;
					parent[s][cur_mask] = nxt_mask;
				}
			}
		}
		stage = NUM_STAGES - s;
		max_bias = frac(0);
		for (int i = 0; i < (1 << NUM_BITS); i++)
			if (max_bias <= dp[s][i])
				max_bias = dp[s][i];
		trace(stage, max_bias);
	}

	if (max_bias.x == 0) {
		cerr << "No path found" << endl;
		return 0;
	}

	vector<int> final_parent(1 << NUM_BITS, -1);
	int min_cnt = UNITS, plaintext;
	for (int i = 0; i < (1 << NUM_BITS); i++) {
		int j = i, cnt = 0;
		for (int k = 0; k < NUM_STAGES - 1; k++)
			j = parent[k][j];
		final_parent[i] = j;
		if (!(max_bias == dp[0][i])) continue;
		for (int k = UNITS - 1; k >= 0; k--) {
			if (j % (1 << SBOX_SZ)) cnt++;
			j >>= SBOX_SZ;
		}
		if (cnt < min_cnt)
			plaintext = i, min_cnt = cnt;
	}

	draw(plaintext);

	cerr << "Active SBoxes : ";
	cerr << min_cnt << endl;
	
	cerr << "Plaintext bits :";
	for (int i = 0; i < NUM_BITS; i++)
		if ((plaintext >> i) & 1) cerr << " " << i;
	cerr << endl;

	cerr << "Ciphertext bits :";
	for (int i = 0; i < NUM_BITS; i++)
		if ((final_parent[plaintext] >> i) & 1) 
			cerr << " " << i;
	cerr << endl;
}

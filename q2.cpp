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
#define popcnt __builtin_popcount
#define pii pair<int, int>
#define mp make_pair
#define eb emplace_back
#define all(v) (v).begin(), (v).end()
#define rall(v) (v).rbegin(), (v).rend()

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
int KEY_BITS;
int plain_bits, cipher_bits;

vector<int> SBOX, PBOX;
vector<int> SINV, PINV;

vector<vector<frac>> bias, dp;
vector<vector<int>> parent;
vector<int> table;
vector<pii> data;

string binary(int p) {
	string s(SBOX_BITS, '0');
	for (int i = 0; i < SBOX_BITS; i++)
		if ((p >> i) & 1)
			s[SBOX_BITS - i - 1] = '1';
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
int pass_box(int plain, int sz, vector<int>& box) {
	vector<int> vec(sz);
	for (int i = 0; i < sz; i++) {
		vec[i] = box[plain & (SBOX_SZ - 1)];
		plain >>= SBOX_BITS;
	}
	int ret = 0;
	for (int i = sz - 1; i >= 0; i--)
		ret = (ret << SBOX_BITS) + vec[i];
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
		print(permute(plain, PINV));
		cout << endl;
	}
	print(plain);
	cout << sep << endl;
}
void random_init() {
	NUM_BITS = 12, NUM_STAGES = 5, SBOX_BITS = 4;
	assert(NUM_BITS % SBOX_BITS == 0);
	UNITS = NUM_BITS / SBOX_BITS;
	SBOX_SZ = 1 << SBOX_BITS;
	SBOX.resize(SBOX_SZ);
	PBOX.resize(NUM_BITS);
	iota(all(SBOX), 0), iota(all(PBOX), 0);
	shuf(SBOX), shuf(PBOX);
}
void custom_init() {
	SBOX = {4, 3, 6, 11, 12, 14, 13, 10, 5, 15, 8, 7, 1, 2, 9, 0};
	PBOX = {1, 7, 9, 11, 6, 0, 8, 3, 2, 10, 4, 5};
	NUM_STAGES = 5;
	NUM_BITS = PBOX.size();
	SBOX_SZ = SBOX.size();
	SBOX_BITS = int(log2(SBOX_SZ));
	UNITS = NUM_BITS / SBOX_BITS;
	assert((1 << SBOX_BITS) == SBOX_SZ);
	assert(NUM_BITS % SBOX_BITS == 0);
}
void generate(int key) {
	data.clear();
	for (int plain = 0; plain < (1 << NUM_BITS); plain++) {
		int cipher = plain;
		for (int s = 0; s < NUM_STAGES - 1; s++) {
			cipher ^= key;
			cipher = pass_box(cipher, UNITS, SBOX);
			cipher = permute(cipher, PBOX);
		}
		cipher ^= key;	
		cipher = pass_box(cipher, UNITS, SBOX);
		cipher ^= key;
		data.eb(plain, cipher);
	}
	shuf(data);
}
int get_sbox_mask(int p, int mask) {
	vector<int> sbox_mask;
	for (int i = UNITS - 1; i >= 0; i--) {
		int c = mask & (SBOX_SZ - 1);
		int d = p & (SBOX_SZ - 1);
		if (c > 0) sbox_mask.eb(d);
		p >>= SBOX_BITS; 
		mask >>= SBOX_BITS; 
	}
	int active = sbox_mask.size(), ret = 0;
	for (int i = active - 1; i >= 0; i--)
		ret = (ret << SBOX_BITS) + sbox_mask[i];
	return ret;
}
void update_scores(int plain, int cipher) {
	int new_cipher = get_sbox_mask(cipher, cipher_bits);
	int new_cipher_bits = get_sbox_mask(cipher_bits, cipher_bits);
	int active = 0, tmp_cipher_bits = cipher_bits;
	for (int i = 0; i < UNITS; i++) {
		if (tmp_cipher_bits & (SBOX_SZ - 1)) active++;
		tmp_cipher_bits >>= SBOX_BITS;
	}
	for (int k = 0; k < (1 << KEY_BITS); k++) {
		int bits = pass_box(new_cipher ^ k, active, SINV);
		bits = popcnt(bits & new_cipher_bits);
		bits += popcnt(plain & plain_bits);
		if (bits % 2 == 0) table[k]++;
	}
}
int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	cout.tie(NULL);

	random_init();
	// custom_init();

	trace(NUM_BITS, NUM_STAGES, SBOX_BITS);
	trace(SBOX);
	trace(PBOX);

	// initialize inverse boxes
	PINV.resize(NUM_BITS); SINV.resize(SBOX_SZ);
	for (int i = 0; i < NUM_BITS; i++)
		PINV[PBOX[i]] = i;
	for (int i = 0; i < SBOX_SZ; i++)
		SINV[SBOX[i]] = i;

	// calculate biases
	bias.resize(SBOX_SZ);
	bias[0].assign(SBOX_SZ, 0);
	for (int i = 1; i < SBOX_SZ; i++) {
		bias[i].assign(SBOX_SZ, 0);
		for (int j = 1; j < SBOX_SZ; j++) {
			int cnt = 0, cur;
			for (int k = 0; k < SBOX_SZ; k++) {
				cur = popcnt(SBOX[k] & j);
				cur += popcnt(k & i);
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
			rev_mask = permute(nxt_mask, PINV), bits = 0;
			for (ptr = UNITS - 1; ptr >= 0; ptr--) {
				nxt_sboxes_mask[ptr] = rev_mask & (SBOX_SZ - 1);
				rev_mask >>= SBOX_BITS;
				if (nxt_sboxes_mask[ptr]) bits += SBOX_BITS;
			}
			for (mask = 1; mask < (1 << bits); mask++) {
				sub_mask = mask, cur_mask = 0, bad = 0;
				current_total_bias = dp[s + 1][nxt_mask];
				for (ptr = UNITS - 1; ptr >= 0; ptr--) {
					if (nxt_sboxes_mask[ptr] == 0) continue;
					inp_mask = sub_mask & (SBOX_SZ - 1);
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
	int min_cnt = UNITS + 1;
	for (int i = 0; i < (1 << NUM_BITS); i++) {
		int j = i, cnt = 0;
		for (int k = 0; k < NUM_STAGES - 1; k++)
			j = parent[k][j];
		final_parent[i] = j;
		if (!(max_bias == dp[0][i])) continue;
		for (int k = 0; k < UNITS; k++) {
			if (j & (SBOX_SZ - 1)) cnt++;
			j >>= SBOX_BITS;
		}
		if (cnt < min_cnt)
			plain_bits = i, min_cnt = cnt;
	}

	draw(plain_bits);
	cipher_bits = final_parent[plain_bits];

	cout << "Active SBoxes : ";
	cout << min_cnt << endl;

	cout << "Plaintext bits :";
	for (int i = NUM_BITS - 1; i >= 0; i--)
		if ((plain_bits >> i) & 1) cout << " " << (NUM_BITS - i - 1);
	cout << endl;

	cout << "Ciphertext bits :";
	for (int i = NUM_BITS - 1; i >= 0; i--)
		if ((cipher_bits >> i) & 1) cout << " " << (NUM_BITS - i - 1);
	cout << endl;

	KEY_BITS = min_cnt * SBOX_BITS;
	table.assign(1 << KEY_BITS, 0);

	int key = abs((ll)rng()) % (1 << NUM_BITS);
	key = 2843;
	generate(key);	

	int samples = data.size();

	for (int i = 0; i < samples; i++) {
		int plaintext = data[i].first;
		int ciphertext = data[i].second;
		update_scores(plaintext, ciphertext);
	}

	vector<pii> candidates(1 << KEY_BITS);
	for (int i = 0; i < (1 << KEY_BITS); i++) {
		candidates[i].first = abs(2 * table[i] - samples);
		candidates[i].second = i;
	}
	sort(rall(candidates));

	int subkey = get_sbox_mask(key, cipher_bits);
	trace(key, subkey);
	trace(candidates);
}

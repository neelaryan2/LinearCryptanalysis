#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

#pragma comment(linker, "/stack:200000000")
#pragma GCC optimize("Ofast,unroll-loops,no-stack-protector,fast-math")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

using ll = long long;
#define popcnt __builtin_popcount
#define pii pair<int, int>
#define mp make_pair
#define eb emplace_back
#define all(v) (v).begin(), (v).end()
#define rall(v) (v).rbegin(), (v).rend()

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
	for (int i = 0; i < NUM_STAGES; i++) {
		print(plain);
		for (int j = 0; j < UNITS; j++) {
			cout << dash;
			if (j < UNITS - 1) cout << '+';
		}
		cout << endl;
		plain = parent[i][plain];
		if (i == NUM_STAGES - 1)
			print(plain);
		else
			print(permute(plain, PINV));
		cout << endl;
	}
	print(plain);
	cout << sep << endl;
}
vector<int> get_bits(int p) {
	vector<int> ret;
	for (int i = NUM_BITS - 1; i >= 0; i--) {
		if ((p >> i) & 1) {
			ret.eb(NUM_BITS - i - 1);
		}
	}
	return ret;
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
void init() {
	cin >> NUM_STAGES;
	cin >> NUM_BITS;
	PBOX.resize(NUM_BITS);
	for (int i = 0; i < NUM_BITS; i++)
		cin >> PBOX[i];
	cin >> SBOX_BITS;
	SBOX_SZ = (1 << SBOX_BITS);
	SBOX.resize(SBOX_SZ);
	for (int i = 0; i < SBOX_SZ; i++)
		cin >> SBOX[i];
	assert(NUM_BITS % SBOX_BITS == 0);
	UNITS = NUM_BITS / SBOX_BITS;

	// initialize inverse boxes
	PINV.resize(NUM_BITS); SINV.resize(SBOX_SZ);
	for (int i = 0; i < NUM_BITS; i++)
		PINV[PBOX[i]] = i;
	for (int i = 0; i < SBOX_SZ; i++)
		SINV[SBOX[i]] = i;

	// initialise dp and parent
	parent.resize(NUM_STAGES), dp.resize(NUM_STAGES);
	for (int i = 0; i < NUM_STAGES; i++) {
		dp[i].assign(1 << NUM_BITS, 0);
		parent[i].assign(1 << NUM_BITS, -1);
	}
	dp[NUM_STAGES - 1].assign(1 << NUM_BITS, frac(1, 1));
	iota(all(parent[NUM_STAGES - 1]), 0);
}
int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	cout.tie(NULL);

	init();

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

	// declare auxiliary variables outside BIG for loop
	frac max_bias(0), current_total_bias(0);
	vector<int> nxt_sboxes_mask(UNITS);

	int rev_mask, bits, stage;
	int sub_mask, cur_mask, bad;
	int inp_mask, out_mask;
	int nxt_mask, mask, ptr;

	// calculate dp and parents
	for (int s = NUM_STAGES - 2; s >= 0; s--) {
		max_bias = frac(0), nxt_mask = 1;
		for (mask = 1; mask < (1 << NUM_BITS); mask++)
			if (dp[s + 1][nxt_mask] < dp[s + 1][mask])
				nxt_mask = mask;
		vector<int> nxts;
		for (mask = 1; mask < (1 << NUM_BITS); mask++)
			if (dp[s + 1][nxt_mask] == dp[s + 1][mask])
				nxts.eb(mask);
		for (int nxt_mask : nxts) {
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
	}

	if (max_bias.x == 0) {
		cout << "Every path has 0 bias. Exiting...." << endl;
		return 0;
	}

	vector<int> final_parent(1 << NUM_BITS, -1);
	int min_cnt = UNITS + 1;
	for (int i = 1; i < (1 << NUM_BITS); i++) {
		if (!(max_bias == dp[0][i])) continue;
		int j = i, cnt = 0;
		for (int k = 0; k < NUM_STAGES - 1; k++)
			j = parent[k][j];
		final_parent[i] = j;
		for (int k = 0; k < UNITS; k++) {
			if (j & (SBOX_SZ - 1)) cnt++;
			j >>= SBOX_BITS;
		}
		if (cnt < min_cnt)
			plain_bits = i, min_cnt = cnt;
	}

	cipher_bits = final_parent[plain_bits];

	vector<int> vec = get_bits(plain_bits);
	for (int i : vec) cout << 'P' << i << ", ";

	int keybits = plain_bits;
	for (int k = 0; k < NUM_STAGES; k++) {
		vec = get_bits(keybits);
		for (int i : vec) cout << 'K' << k << i << ", ";
		keybits = parent[k][keybits];
	}
	vec = get_bits(cipher_bits);
	for (int i : vec) {
		cout << 'C' << i;
		if (i != vec.back()) 
			cout << ", ";
		else
			cout << endl;
	}
	cout << "Bias = " << max_bias << endl;

	cout << endl;
	draw(plain_bits);
}

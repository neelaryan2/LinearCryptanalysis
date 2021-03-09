template <class T>
ostream& operator<<(ostream& os, vector<T> V) {
    os << "[ ";
    for (auto v : V) os << v << ", ";
    return os << "]";
}
template <class L, class R>
ostream& operator<<(ostream& os, pair<L, R> P) {
    return os << "(" << P.first << "," << P.second << ")";
}
#define trace(...) __f(#__VA_ARGS__, __VA_ARGS__)
template <typename Arg1>
void __f(const char* name, Arg1&& arg1) {
    cout << name << " : " << arg1 << endl;
}
template <typename Arg1, typename... Args>
void __f(const char* names, Arg1&& arg1, Args&&... args) {
    const char* comma = strchr(names + 1, ',');
    cout.write(names, comma - names) << " : " << arg1 << " | ";
    __f(comma + 1, args...);
}
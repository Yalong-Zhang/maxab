#include <bits/stdc++.h>
typedef long long ll;
using namespace std;
const int INF = 2000000000;
const double epsilon = 0.0000000001;
inline bool double_equal(double x, double y) { return x >= y - epsilon && x <= y + epsilon; }

struct Timer {
	double start_time, end_time;
	void start() { start_time = clock(); }
	void end() { end_time = clock(); }
	double time() { return (end_time - start_time) / CLOCKS_PER_SEC; }
};
template <class T>
struct Set {
	T* nodes; bool* in; int size = -1;
	Set() {}
	Set(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(T)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void alloc(int sz) { size = 0; nodes = (T*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); }
	void insert(T x) { nodes[size++] = x; in[x] = true; }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false; size = 0; }
	~Set() { free(nodes), free(in); }
};
template <class T>
struct Map {
	int* nodes; bool* in; int size = -1; T* value;
	Map() {}
	Map(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void alloc(int sz) { size = 0; nodes = (int*)malloc(sz * sizeof(int)); in = (bool*)malloc(sz * sizeof(int)); memset(in, 0, sz * sizeof(int)); value = (T*)malloc(sz * sizeof(T)); memset(value, 0, sz * sizeof(T)); }
	void freememory() { free(nodes), free(in), free(value); }
	void clear() { for (int i = 0; i < size; i++) in[nodes[i]] = false, value[nodes[i]] = 0; size = 0; }
	T& operator[](int x) { if (!in[x]) nodes[size++] = x, in[x] = true; return value[x]; }
	~Map() { free(nodes), free(in), free(value); }
};
template <class T>
struct Queue {
	T* nodes; int head, tail;
	Queue() {}
	Queue(int size) { nodes = (T*)malloc(size * sizeof(T)); head = tail = 0; }
	void alloc(int sz) { head = tail = 0; nodes = (int*)malloc(sz * sizeof(int)); }
	~Queue() { free(nodes); }
	bool empty() { return head == tail; }
	int pop() { return nodes[head++]; }
	void push(T x) { nodes[tail++] = x; }
	void clear() { head = tail = 0; }
};

inline int read_number(FILE* in) {
	int x = 0; char ch = 0; while (ch < '0' || ch > '9') ch = fgetc(in); while (ch >= '0' && ch <= '9') { x = x * 10 + (ch - '0'); ch = fgetc(in); } return x;
}
inline void check(bool flag, const char* message) {
	if (!flag) {
		printf("!!!!! CHECK ERROR !!!!!\n");
		printf("Error message: %s\n", message);
		assert(0);
	}
}

enum Algorithm {
	SINGLE, MULTIPLE, MULTIPLE1, MULTIPLE2, ALL
};
Algorithm algorithm_used;

struct Edge { int u, v, to; };
struct Graph {
	int M, U, V, N, dmax_U, dmax_V;
	Edge* e;
	int* undeg, * indeg;
	int** adj;
	void read_graph_from_dataset(char* dataset_address);
	int* undeg_order_U, * undeg_order_V;

	inline bool in_U(int x) { return x < U; }
	inline bool in_V(int x) { return x >= U; }

	void get_core();

	Set<int> C;
	void get_a_core(int core_alpha, int core_beta, Set<int>& C);

	void output_information(Set<int>& S);

	int* node2p_primary, * node2p_secondary, * p2node_primary, * p2node_secondary, * d_start_primary, * d_start_secondary, * tem_undeg, * core;
	int* core_alpha, * core_beta, * core_start_alpha, * core_start_beta;
	int get_a_column(bool now_alpha, int now_p);

	int degeneracy;
	void get_single_core(Set<int>& C);

	double eps, gamma; int core_count;
	double now_max_density;
	void get_a_line(double k);

	void check_correctness();
	const int OUTPUT_NODES_LIMIT = 10;
};
Graph G;

bool cmp(int x, int y) { return G.undeg[x] < G.undeg[y]; }
void Graph::read_graph_from_dataset(char* dataset_address) {
	FILE* in = fopen(dataset_address, "r");
	check(in != NULL, "Can not open file dataset_address\n");

	int num1, num2, num3, count;
	char line[256];
	fgets(line, sizeof(line), in);
	count = sscanf(line, "%d %d %d", &num1, &num2, &num3);
	if (count == 2) {
		if (num2 > num1) swap(num1, num2);
		M = num1, U = V = num2;
	}
	else if (count == 3)
		M = num1, U = num2, V = num3;
	else
		check(0, "dataset format error");

	U++, V++; N = U + V;

	e = (Edge*)malloc(M * sizeof(Edge));
	undeg = (int*)malloc(N * sizeof(int)); indeg = (int*)malloc(N * sizeof(int));
	memset(undeg, 0, N * sizeof(int)); memset(indeg, 0, N * sizeof(int));
	adj = (int**)malloc(N * sizeof(int*));
	p2node_primary = (int*)malloc(N * sizeof(int)), p2node_secondary = (int*)malloc(N * sizeof(int)), node2p_primary = (int*)malloc(N * sizeof(int)), node2p_secondary = (int*)malloc(N * sizeof(int)), core = (int*)malloc(N * sizeof(int)), tem_undeg = (int*)malloc(N * sizeof(int));
	core_alpha = (int*)malloc(N * sizeof(int)), core_beta = (int*)malloc(N * sizeof(int));
	undeg_order_U = (int*)malloc(N * sizeof(int)); undeg_order_V = (int*)malloc(N * sizeof(int));
	C.alloc(N);

	for (int i = 0; i < M; i++) {
		e[i].u = read_number(in), e[i].v = read_number(in) + U;
		undeg[e[i].u]++, undeg[e[i].v]++;
	}
	for (int i = 0; i < N; i++) {
		if (in_U(i)) dmax_U = max(dmax_U, undeg[i]);
		else dmax_V = max(dmax_V, undeg[i]);
		adj[i] = (int*)malloc(undeg[i] * sizeof(Edge));
	}
	memset(undeg, 0, N * sizeof(int));
	for (int i = 0; i < M; i++) {
		int u = e[i].u, v = e[i].v;
		adj[u][undeg[u]++] = i;
		adj[v][undeg[v]++] = i;
	}
	d_start_primary = (int*)malloc((max(dmax_U, dmax_V) + 5) * sizeof(int)), d_start_secondary = (int*)malloc((max(dmax_U, dmax_V) + 5) * sizeof(int));
	core_start_alpha = (int*)malloc((max(dmax_U, dmax_V) + 5) * sizeof(int)), core_start_beta = (int*)malloc((max(dmax_U, dmax_V) + 5) * sizeof(int));

	for (int i = 0; i < U; i++) undeg_order_U[i] = i;
	for (int i = U; i < N; i++) undeg_order_V[i] = i;
	sort(undeg_order_U, undeg_order_U + U, &cmp);
	sort(undeg_order_V + U, undeg_order_V + N, &cmp);

	return;
}
int Graph::get_a_column(bool now_alpha, int now_p) {
	memcpy(tem_undeg, undeg, N * sizeof(int));
	memset(core, -1, N * sizeof(int));
	/*
	for (int x = 0; x < N; x++) {
		if (now_alpha ? in_U(x) : in_V(x)) p2node_primary[x] = x;
		else p2node_secondary[x] = x;
	}
	sort(p2node_primary + (now_alpha ? 0 : U), p2node_primary + (now_alpha ? U : N), &cmp);
	sort(p2node_secondary + (now_alpha ? U : 0), p2node_secondary + (now_alpha ? N : U), &cmp);*/
	if (now_alpha) {
		memcpy(p2node_primary, undeg_order_U, U * sizeof(int));
		memcpy(p2node_secondary + U, undeg_order_V + U, V * sizeof(int));
	}
	else {
		memcpy(p2node_primary + U, undeg_order_V + U, V * sizeof(int));
		memcpy(p2node_secondary, undeg_order_U, U * sizeof(int));
	}
	int nowr = -1, pointer, len_primary, len_secondary;
	pointer = now_alpha ? 0 : U, len_primary = now_alpha ? U : N;
	while (pointer < len_primary) {
		int x = p2node_primary[pointer];
		node2p_primary[x] = pointer;
		if (tem_undeg[x] > nowr) d_start_primary[++nowr] = pointer;
		else pointer++;
	}
	nowr = -1, pointer = now_alpha ? U : 0, len_secondary = now_alpha ? N : U;
	while (pointer < len_secondary) {
		int x = p2node_secondary[pointer];
		node2p_secondary[x] = pointer;
		if (tem_undeg[x] > nowr) d_start_secondary[++nowr] = pointer;
		else pointer++;
	}
	nowr = 0; int cnt = 0, pointer_primary = now_alpha ? 0 : U, pointer_secondary = now_alpha ? U : 0;
	int now_U = U, now_V = V, now_E = M;
	while (cnt < N) {
		bool flag = false;
		int now = pointer_primary < len_primary ? p2node_primary[pointer_primary] : 0;
		if (pointer_primary < len_primary && tem_undeg[now] <= nowr) {
			if (now_alpha) now_U--; else now_V--;
			core[now] = nowr, flag = true, cnt++, pointer_primary++;
			for (int j = 0; j < undeg[now]; j++) {
				Edge& ne = e[adj[now][j]];
				int tar = now_alpha ? ne.v : ne.u;
				if (core[tar] != -1) continue;
				now_E--;
				int lp = d_start_primary[tem_undeg[now]], rp = node2p_primary[now], ln = p2node_primary[lp], rn = now;
				node2p_primary[ln] = rp;
				node2p_primary[rn] = lp;
				p2node_primary[lp] = rn;
				p2node_primary[rp] = ln;
				d_start_primary[tem_undeg[now]]++;
				tem_undeg[now]--;

				lp = d_start_secondary[tem_undeg[tar]], rp = node2p_secondary[tar], ln = p2node_secondary[lp], rn = tar;
				node2p_secondary[ln] = rp;
				node2p_secondary[rn] = lp;
				p2node_secondary[lp] = rn;
				p2node_secondary[rp] = ln;
				d_start_secondary[tem_undeg[tar]]++;
				tem_undeg[tar]--;
			}
		}
		now_max_density = max(now_max_density, now_E / (sqrt(now_U) * sqrt(now_V)));
		now = pointer_secondary < len_secondary ? p2node_secondary[pointer_secondary] : 0;
		if (pointer_secondary < len_secondary && tem_undeg[now] < now_p) {
			if (now_alpha) now_V--; else now_U--;
			core[now] = nowr, flag = true, cnt++, pointer_secondary++;
			for (int j = 0; j < undeg[now]; j++) {
				Edge& ne = e[adj[now][j]];
				int tar = now_alpha ? ne.u : ne.v;
				if (core[tar] != -1) continue;
				now_E--;
				int lp = d_start_secondary[tem_undeg[now]], rp = node2p_secondary[now], ln = p2node_secondary[lp], rn = now;
				node2p_secondary[ln] = rp;
				node2p_secondary[rn] = lp;
				p2node_secondary[lp] = rn;
				p2node_secondary[rp] = ln;
				d_start_secondary[tem_undeg[now]]++;
				tem_undeg[now]--;

				lp = d_start_primary[tem_undeg[tar]], rp = node2p_primary[tar], ln = p2node_primary[lp], rn = tar;
				node2p_primary[ln] = rp;
				node2p_primary[rn] = lp;
				p2node_primary[lp] = rn;
				p2node_primary[rp] = ln;
				d_start_primary[tem_undeg[tar]]++;
				tem_undeg[tar]--;
			}
		}
		now_max_density = max(now_max_density, now_E / (sqrt(now_U) * sqrt(now_V)));
		if (!flag)
			nowr++;
	}
	int* core_x = now_alpha ? core_alpha : core_beta, * core_start = now_alpha ? core_start_alpha : core_start_beta;
	cnt = 0, nowr = 0, pointer_primary = now_alpha ? 0 : U, pointer_secondary = now_alpha ? U : 0;
	while (pointer_primary < len_primary && pointer_secondary < len_secondary) {
		core_start[nowr] = cnt;
		while (pointer_primary < len_primary && core[p2node_primary[pointer_primary]] <= nowr)
			core_x[cnt++] = p2node_primary[pointer_primary++];
		while (pointer_secondary < len_secondary && core[p2node_secondary[pointer_secondary]] <= nowr)
			core_x[cnt++] = p2node_secondary[pointer_secondary++];
		nowr++;
	}
	return nowr - 1;
}
void Graph::get_a_core(int core_alpha, int core_beta, Set<int>& C) {
	C.clear();

	// a node --> will_be_deleted_nodes --> deleted_nodes
	Set<int> will_be_deleted_nodes; will_be_deleted_nodes.alloc(N); int pointer = 0;
	Set<int> deleted_nodes; deleted_nodes.alloc(N);
	int* temporary_undeg = (int*)malloc(N * sizeof(int)); memcpy(temporary_undeg, undeg, N * sizeof(int));

	for (int x = 0; x < N; x++)
		if (undeg[x] < (in_U(x) ? core_alpha : core_beta))
			will_be_deleted_nodes.insert(x);

	while (pointer < will_be_deleted_nodes.size) {
		int x = will_be_deleted_nodes.nodes[pointer++];
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			int y = in_U(x) ? ne.v : ne.u;
			if (deleted_nodes.in[y]) continue;
			if (temporary_undeg[y] == (in_U(y) ? core_alpha : core_beta)) {
				will_be_deleted_nodes.insert(y);
			}
			temporary_undeg[x]--;
			temporary_undeg[y]--;
		}
		deleted_nodes.insert(x);
	}

	for (int x = 0; x < N; x++)
		if (!deleted_nodes.in[x])
			C.insert(x);

	free(temporary_undeg);
	return;
}
void Graph::get_core() {
	if (algorithm_used == MULTIPLE)
		gamma = max(U, V);
	Timer timer;
	if (algorithm_used == SINGLE) {
		timer.start();
		get_single_core(C);
		timer.end();
		printf("- %-20s: %d\n", "Degeneracy", degeneracy);
	}
	else if (algorithm_used == MULTIPLE || algorithm_used == MULTIPLE1) {
		timer.start();
		now_max_density = 0.0;
		double t = eps * eps + 4 * eps + 2 + sqrt(eps * eps * eps * eps + 8 * eps * eps * eps + 20 * eps * eps + 16 * eps); check(t > 2, "t error");
		int max_exp = 1;
		while (pow(t / 2, max_exp) < gamma) max_exp += 2;
		printf("LOG: max_exp = %d\n", max_exp);
		for (int exp = max_exp; exp >= 3; exp -= 2)
			get_a_line(pow(t / 2, exp - 1));
		get_a_line(1.0);
		for (int exp = 3; exp <= max_exp; exp += 2)
			get_a_line(pow(2 / t, exp - 1));
		timer.end();
	}
	else if (algorithm_used == MULTIPLE2) {
		timer.start();
		now_max_density = 0.0;
		double t = eps * eps + 4 * eps + 2 + sqrt(eps * eps * eps * eps + 8 * eps * eps * eps + 20 * eps * eps + 16 * eps); check(t > 2, "t error");
		int number_of_round = core_count;
		for (int exp = number_of_round; exp >= 3; exp -= 2)
			get_a_line(pow(t / 2, exp - 1));
		get_a_line(1.0);
		for (int exp = 3; exp <= number_of_round; exp += 2)
			get_a_line(pow(2 / t, exp - 1));
		timer.end();
	}
	else if (algorithm_used == ALL) {
		timer.start();
		now_max_density = 0.0;
		int now_p = 1, maxa, maxb, degeneracy;
		ll maxab = 0;
		while (true) {
			printf("LOG: now_p = %d\n", now_p);
			int alpha_max, beta_max;
			alpha_max = get_a_column(true, now_p);
			if ((ll)alpha_max * now_p > maxab) maxab = (ll)alpha_max * now_p, maxa = alpha_max, maxb = now_p;
			beta_max = get_a_column(false, now_p);
			if ((ll)now_p * beta_max > maxab) maxab = (ll)now_p * beta_max, maxa = now_p, maxb = beta_max;
			if (alpha_max < now_p || beta_max < now_p) {
				degeneracy = now_p - 1;
				break;
			}
			now_p++;
		}
		printf("- %-20s: %d, %d\n", "maxa, maxb", maxa, maxb);
		printf("- %-20s: %lf\n", "Max density", now_max_density);
		get_a_core(maxa, maxb, C);
		timer.end();
	}
	output_information(C);
	printf("- %-20s: %lf\n", "Get core time", timer.time());
	return;
}
void Graph::get_single_core(Set<int>& C) {
	int* p2node = (int*)malloc(N * sizeof(int)), * node2p = (int*)malloc(N * sizeof(int)), * d_start = (int*)malloc(N * sizeof(int));
	for (int x = 0; x < N; x++)
		tem_undeg[x] = undeg[x], p2node[x] = x, core[x] = -1;
	sort(p2node, p2node + N, &cmp);
	int nowr = -1, pointer = 0;
	while (pointer < N)
	{
		node2p[p2node[pointer]] = pointer;
		if (tem_undeg[p2node[pointer]] > nowr)
			d_start[++nowr] = pointer;
		else
			pointer++;
	}
	nowr = 0;
	int now_U = U, now_V = V, now_E = M, now_max_index = 0;
	double now_density = 0.0;
	for (pointer = 0; pointer < N; pointer++)
	{
		if (now_E / (sqrt(now_U) * sqrt(now_V)) > now_density)
			now_density = now_E / (sqrt(now_U) * sqrt(now_V)), now_max_index = pointer;
		int now = p2node[pointer];
		if (core[now] != -1) continue;
		nowr = max(nowr, tem_undeg[now]);
		core[now] = nowr;
		if (in_U(now)) now_U--; else now_V--;
		for (int j = 0; j < undeg[now]; j++)
		{
			Edge& ne = e[adj[now][j]];
			int tar = ne.u == now ? ne.v : ne.u;
			if (core[tar] != -1) continue;
			now_E--;
			int lp = d_start[tem_undeg[now]], rp = node2p[now], ln = p2node[lp], rn = now;
			node2p[ln] = rp;
			node2p[rn] = lp;
			p2node[lp] = rn;
			p2node[rp] = ln;
			d_start[tem_undeg[now]]++;
			tem_undeg[now]--;

			lp = d_start[tem_undeg[tar]], rp = node2p[tar], ln = p2node[lp], rn = tar;
			node2p[ln] = rp;
			node2p[rn] = lp;
			p2node[lp] = rn;
			p2node[rp] = ln;
			d_start[tem_undeg[tar]]++;
			tem_undeg[tar]--;
		}
	}
	degeneracy = nowr;
	C.clear();
	for (int i = now_max_index; i < N; i++)
		C.insert(p2node[i]);
	return;
}
void Graph::get_a_line(double k) {
	memcpy(tem_undeg, undeg, N * sizeof(int));
	memset(core, -1, N * sizeof(int));
	for (int x = 0; x < U; x++) p2node_primary[x] = x;
	for (int x = U; x < N; x++) p2node_secondary[x] = x;

	memcpy(p2node_primary, undeg_order_U, U * sizeof(int));
	memcpy(p2node_secondary + U, undeg_order_V + U, V * sizeof(int));
	int nowd = -1, pointer = 0;
	while (pointer < U) {
		int x = p2node_primary[pointer];
		node2p_primary[x] = pointer;
		if (tem_undeg[x] > nowd) d_start_primary[++nowd] = pointer;
		else pointer++;
	}
	nowd = -1, pointer = U;
	while (pointer < N) {
		int x = p2node_secondary[pointer];
		node2p_secondary[x] = pointer;
		if (tem_undeg[x] > nowd) d_start_secondary[++nowd] = pointer;
		else pointer++;
	}
	int cnt = 0, pointer_primary = 0, pointer_secondary = U, now_alpha = 0, now_beta = 0;
	int now_U = U, now_V = V, now_E = M, now_max_index_primary = 0, now_max_index_secondary = U; bool density_updated = false;
	double max_density_of_this_line = 0.0;
	while (cnt < N) {
		bool flag = false;
		int now = pointer_primary < U ? p2node_primary[pointer_primary] : 0;
		if (pointer_primary < U && tem_undeg[now] < now_alpha) {
			now_U--;
			core[now] = now_alpha, flag = true, cnt++, pointer_primary++;
			for (int j = 0; j < undeg[now]; j++) {
				Edge& ne = e[adj[now][j]];
				int tar = now_alpha ? ne.v : ne.u;
				if (core[tar] != -1) continue;
				now_E--;
				int lp = d_start_primary[tem_undeg[now]], rp = node2p_primary[now], ln = p2node_primary[lp], rn = now;
				node2p_primary[ln] = rp;
				node2p_primary[rn] = lp;
				p2node_primary[lp] = rn;
				p2node_primary[rp] = ln;
				d_start_primary[tem_undeg[now]]++;
				tem_undeg[now]--;

				lp = d_start_secondary[tem_undeg[tar]], rp = node2p_secondary[tar], ln = p2node_secondary[lp], rn = tar;
				node2p_secondary[ln] = rp;
				node2p_secondary[rn] = lp;
				p2node_secondary[lp] = rn;
				p2node_secondary[rp] = ln;
				d_start_secondary[tem_undeg[tar]]++;
				tem_undeg[tar]--;
			}
		}
		if (now_max_density < now_E / (sqrt(now_U) * sqrt(now_V)))
			density_updated = true, now_max_density = now_E / (sqrt(now_U) * sqrt(now_V)), now_max_index_primary = pointer_primary, now_max_index_secondary = pointer_secondary;
		max_density_of_this_line = max(max_density_of_this_line, now_E / (sqrt(now_U) * sqrt(now_V)));
		now = pointer_secondary < N ? p2node_secondary[pointer_secondary] : 0;
		if (pointer_secondary < N && tem_undeg[now] < now_beta) {
			now_V--;
			core[now] = now_beta, flag = true, cnt++, pointer_secondary++;
			for (int j = 0; j < undeg[now]; j++) {
				Edge& ne = e[adj[now][j]];
				int tar = now_alpha ? ne.u : ne.v;
				if (core[tar] != -1) continue;
				now_E--;
				int lp = d_start_secondary[tem_undeg[now]], rp = node2p_secondary[now], ln = p2node_secondary[lp], rn = now;
				node2p_secondary[ln] = rp;
				node2p_secondary[rn] = lp;
				p2node_secondary[lp] = rn;
				p2node_secondary[rp] = ln;
				d_start_secondary[tem_undeg[now]]++;
				tem_undeg[now]--;

				lp = d_start_primary[tem_undeg[tar]], rp = node2p_primary[tar], ln = p2node_primary[lp], rn = tar;
				node2p_primary[ln] = rp;
				node2p_primary[rn] = lp;
				p2node_primary[lp] = rn;
				p2node_primary[rp] = ln;
				d_start_primary[tem_undeg[tar]]++;
				tem_undeg[tar]--;
			}
		}
		if (now_max_density < now_E / (sqrt(now_U) * sqrt(now_V)))
			density_updated = true, now_max_density = now_E / (sqrt(now_U) * sqrt(now_V)), now_max_index_primary = pointer_primary, now_max_index_secondary = pointer_secondary;
		max_density_of_this_line = max(max_density_of_this_line, now_E / (sqrt(now_U) * sqrt(now_V)));
		if (!flag) {
			if (double_equal(now_beta, k * now_alpha))
				now_alpha++, now_beta++;
			else if (now_beta > k * now_alpha)
				now_alpha++;
			else if (now_beta < k * now_alpha)
				now_beta++;
			else
				check(0, "now_alpha now_beta error");
		}
	}
	if (density_updated) {
		C.clear();
		for (int i = now_max_index_primary; i < U; i++)
			C.insert(p2node_primary[i]);
		for (int i = now_max_index_secondary; i < N; i++)
			C.insert(p2node_secondary[i]);
	}
	printf("LOG: k = %e, max_density = %lf\n", k, max_density_of_this_line);
	return;
}
void Graph::output_information(Set<int>& S) {
	int Ucnt = 0, Vcnt = 0, Ecnt = 0;
	for (int i = 0; i < S.size; i++) {
		int x = S.nodes[i];
		if (in_U(x)) Ucnt++;
		else { Vcnt++; continue; }
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (!S.in[ne.v]) continue;
			Ecnt++;
		}
	}
	printf("- %-20s: %d, %d, %d\n", "Ucnt, Vcnt, Ecnt", Ucnt, Vcnt, Ecnt);
	printf("- %-20s: %lf\n", "Approximate density", Ecnt / (sqrt(Ucnt) * sqrt(Vcnt)));
	return;
}
void Graph::check_correctness() {
	return;
}

int main(int argc, char** argv) {
	if (argc < 3) {
	argument_error:
		printf("Usage: ./main <dataset_address> <algorithm> <core_reduction>\n");
		printf("algorithm:\n");
		printf("-s: compute single core\n");
		printf("-m: set gamma as n");
		printf("-m1 <eps> <C>: compute multiple core\n");
		printf("-m2 <eps> <core_count>: compute multiple core\n");
		printf("-a: compute all core\n");
		return 0;
	}
	Timer timer;
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);
	if (strcmp(argv[2], "-s") == 0) algorithm_used = SINGLE;
	else if (strcmp(argv[2], "-m") == 0) {
		if (argc != 4) goto argument_error;
		algorithm_used = MULTIPLE;
		G.eps = stod(argv[3]);
		check(G.eps > 0, "eps error");
	}
	else if (strcmp(argv[2], "-m1") == 0) {
		if (argc != 5) goto argument_error;
		algorithm_used = MULTIPLE1;
		G.eps = stod(argv[3]);
		G.gamma = stod(argv[4]);
		check(G.eps > 0 && G.gamma > 0, "eps or gamma error");
	}
	else if (strcmp(argv[2], "-m2") == 0) {
		if (argc != 5) goto argument_error;
		algorithm_used = MULTIPLE2;
		G.eps = stod(argv[3]);
		G.core_count = stoi(argv[4]);
		check(G.eps > 0 && G.core_count >= 1 && G.core_count % 2 == 1, "eps or core_count error");
	}
	else if (strcmp(argv[2], "-a") == 0) algorithm_used = ALL;
	else goto argument_error;

	double runtime;
	printf("----------Now processing %s----------\n", dataset_address);

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d, %d\n", "|E|, |U|, |V|", G.M, G.U, G.V);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	G.get_core();

	return 0;
}
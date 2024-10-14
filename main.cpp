#include <bits/stdc++.h>
using namespace std;
const int INF = 2000000000;

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
	ALL, DIVIDE
};
Algorithm algorithm_used;

struct Edge { int u, v, to; };
struct Graph {
	int M, U, V, N, dmax_U, dmax_V;
	Edge* e;
	int* undeg, * indeg;
	int** adj;
	void read_graph_from_dataset(char* dataset_address);

	inline bool in_U(int x) { return x < U; }
	inline bool in_V(int x) { return x >= U; }
	inline bool in_S(int x) { return indeg[x] < (in_U(x) ? alpha : beta); }
	inline bool in_T(int x) { return indeg[x] > (in_U(x) ? alpha : beta); }

	void initialize_orientation();
	int alpha, beta;

	Set<int> D; int U_num, V_num;
	void test();

	int maxab = 0, maxa, maxb;
	void get_maxab_all();

	int pseudo, line_cnt, test_cnt, binary_cnt, flow_cnt;
	void get_pseudo(); double get_pseudo_time;
	int* upper_alpha, * upper_beta;
	void get_maxab_divide();
	void get_maxab_dividea();
	void get_maxab_divideb();
	void get_subgraph(); Set<int> reached;
	void divide(int alpha_l, int alpha_u, bool now_alpha);

	double line_a(int alpha_1, int beta_1, int alpha_2, int beta_2, int alphap) {
		check(alpha_1 != alpha_2, "line_a error");
		return double((alphap - alpha_2) * beta_1 + (alpha_1 - alphap) * beta_2) / (alpha_1 - alpha_2);
	}
	double line_b(int alpha_1, int beta_1, int alpha_2, int beta_2, int betap) {
		check(beta_1 != beta_2, "line_b error");
		return double((betap - beta_2) * alpha_1 + (beta_1 - betap) * alpha_2) / (beta_1 - beta_2);
	}
	void update_maxab(int a, int b) {
		if (a * b > maxab) {
			maxab = a * b, maxa = a, maxb = b;
			printf("LOG: new maxab = %d * %d = %d\n", a, b, maxab);
		}
	}

	Set<int> T_node; Queue<int> Q;
	Map<int> parent; Set<int> vis; Map<int> dist, cur;

	int now_alpha_max_indegree, now_beta_max_indegree;
	bool DinicBFS();
	bool DinicDFS(int x);

	void check_correctness();
	const int OUTPUT_NODES_LIMIT = 10;
	void output_core(Set<int>& C);
	void output_dense_subgraph(char*);
	void display();
};
Graph G;

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
	upper_alpha = (int*)malloc(N * sizeof(int)), upper_beta = (int*)malloc(N * sizeof(int));

	T_node.alloc(N); Q.alloc(N); parent.alloc(N); vis.alloc(N); D.alloc(N); dist.alloc(N); cur.alloc(N); reached.alloc(N);

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
}
void Graph::initialize_orientation() {
	for (int i = 0; i < M; i++) {
		if (indeg[e[i].u] < indeg[e[i].v]) e[i].to = e[i].u, indeg[e[i].u]++;
		else e[i].to = e[i].v, indeg[e[i].v]++;
	}
}
void Graph::test() {
	flow_cnt++;
	int pre_alpha = alpha, pre_beta = beta;
	if (alpha > beta) alpha--; else beta--;
	while (DinicBFS()) {
		for (int i = 0; i < T_node.size; i++) {
			int x = T_node.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x);
		}
	}
	alpha = pre_alpha, beta = pre_beta;
	if (T_node.size == 0) {
		D.clear();
		return;
	}
	else {
		while (DinicBFS()) {
			for (int i = 0; i < T_node.size; i++) {
				int x = T_node.nodes[i];
				if (in_T(x))
					parent[x] = -2, cur[x] = 0, DinicDFS(x);
			}
		}
		if (T_node.size != 0) {
			D.clear(), D.insert(0);
			return;
		}
		D.clear(), Q.clear(), reached.clear();
		for (int x = 0; x < N; x++) {
			if (in_S(x)) {
				Q.push(x), reached.insert(x);
			}
		}
		while (!Q.empty()) {
			int x = Q.pop();
			for (int j = 0; j < undeg[x]; j++) {
				Edge& ne = e[adj[x][j]];
				if (ne.to == x) continue;
				if (reached.in[ne.to]) continue;
				Q.push(ne.to), reached.insert(ne.to);
			}
		}
		for (int x = 0; x < N; x++) {
			if (!reached.in[x]) {
				D.insert(x);
				return;
			}
		}
		return;
	}
}
bool Graph::DinicBFS() {
	int dist_t = INF;

	now_alpha_max_indegree = now_beta_max_indegree = 0;
	Q.clear(), dist.clear(), parent.clear(), cur.clear(), T_node.clear();
	for (int x = 0; x < N; x++) {
		if (in_U(x)) now_alpha_max_indegree = max(now_alpha_max_indegree, indeg[x]);
		else now_beta_max_indegree = max(now_beta_max_indegree, indeg[x]);
		if (in_T(x))
			dist[x] = 1, Q.push(x), T_node.insert(x);
	}

	bool break_loop = false;
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to != x) continue;
			int from = in_U(ne.to) ? ne.v : ne.u;
			if (in_S(from)) {
				dist_t = dist[x] + 2; break_loop = true; break;
			}
			if (dist.in[from]) continue;
			dist[from] = dist[x] + 1;
			Q.push(from);
		}
		if (break_loop) break;
	}
	return dist_t != INF;
}
bool Graph::DinicDFS(int x) {
	if (in_S(x)) {
		indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
		return true;
	}
	for (int& j = cur[x]; j < undeg[x]; j++) {
		Edge& ne = e[adj[x][j]];
		if (ne.to != x) continue;
		int from = in_U(ne.to) ? ne.v : ne.u;
		if ((dist[from] != dist[x] + 1) && !in_S(from)) continue;
		parent[from] = adj[x][j];
		if (DinicDFS(from)) {
			if (parent[x] == -2) {
				if (indeg[x] == (in_U(x) ? alpha : beta)) return true;
				continue;
			}
			indeg[x]++, indeg[e[parent[x]].to]--, e[parent[x]].to = x;
			return true;
		}
	}
	return false;
}
void Graph::get_subgraph() {
	while (DinicBFS()) {
		for (int i = 0; i < T_node.size; i++) {
			int x = T_node.nodes[i];
			if (in_T(x))
				parent[x] = -2, cur[x] = 0, DinicDFS(x);
		}
	}
	D.clear(), Q.clear(), reached.clear(); U_num = V_num = 0;
	for (int x = 0; x < N; x++) {
		if (in_S(x)) {
			Q.push(x), reached.insert(x);
		}
	}
	while (!Q.empty()) {
		int x = Q.pop();
		for (int j = 0; j < undeg[x]; j++) {
			Edge& ne = e[adj[x][j]];
			if (ne.to == x) continue;
			if (reached.in[ne.to]) continue;
			Q.push(ne.to), reached.insert(ne.to);
		}
	}
	for (int x = 0; x < N; x++) {
		if (!reached.in[x]) {
			D.insert(x);
			if (in_U(x)) U_num++; else V_num++;
		}
	}
}
void Graph::get_maxab_all() {
	flow_cnt = 0;
	Timer timer;
	int p = 1;
	while (true) {
		int beta_l = pseudo, beta_u = dmax_V;
		while (beta_u > beta_l) {
			int beta_mid = (beta_l + beta_u) / 2;
			alpha = p, beta = beta_mid;
			test();
			if (D.size == 0) beta_u = beta_mid;
			else beta_l = beta_mid + 1;
		}
		beta_l--;
		update_maxab(p, beta_l);
		int alpha_l = pseudo, alpha_u = dmax_V;
		while (alpha_u > alpha_l) {
			int alpha_mid = (alpha_l + alpha_u) / 2;
			alpha = alpha_mid, beta = p;
			test();
			if (D.size == 0) alpha_u = alpha_mid;
			else alpha_l = alpha_mid + 1;
		}
		alpha_l--;
		update_maxab(alpha_l, p);
		printf("- %-20s: %d, %d, %d\n", "p, beta_max, alpha_max", p, beta_l, alpha_l);
		if (beta_l <= p || alpha_l <= p)
			break;
		p++;
	}
	printf("- %-20s: %d\n", "flow count", flow_cnt);
}
void Graph::get_pseudo() {
	int now_indegree, pre_indegree = INF;
	while (true) {
		for (int i = 0; i < M; i++) {
			int from = e[i].to == e[i].u ? e[i].v : e[i].u;
			if (indeg[e[i].to] - indeg[from] >= 2)
				indeg[e[i].to]--, indeg[from]++, e[i].to = from;
		}
		now_indegree = 0;
		for (int x = 0; x < N; x++)
			now_indegree = max(now_indegree, indeg[x]);
		if (now_indegree == pre_indegree) break;
		pre_indegree = now_indegree;
	}
	int appro_pseudo = now_indegree;
	while (true) {
		alpha = appro_pseudo, beta = appro_pseudo;
		test();
		if (D.size != 0) break;
		appro_pseudo--;
	}
	pseudo = appro_pseudo;
	printf("LOG: appro_pseudo = pseudo + %d\n", now_indegree - pseudo);
	return;
}
void Graph::divide(int c_l, int c_u, bool now_alpha) {
	if (c_l > c_u - 2) return;
	int* upper = now_alpha ? upper_alpha : upper_beta;
	int c_m = (c_l + c_u) / 2;
	int d_m;
	d_m = maxab / c_m + 1;
	if (now_alpha) upper[c_m] = line_a(c_l, upper[c_l] + 1, c_u, upper[c_u] + 1, c_m);
	else upper[c_m] = line_b(upper[c_l] + 1, c_l, upper[c_u] + 1, c_u, c_m);
	int condition = 0;
	if (upper[c_m] < d_m) {
		printf("LOG: %c = %d, line %d < %d\n", now_alpha ? 'a' : 'b', c_m, upper[c_m], d_m), line_cnt++;
	}
	else {
		if (now_alpha) alpha = c_m, beta = d_m;
		else alpha = d_m, beta = c_m;
		test();
		if (D.size == 0) {
			upper[c_m] = d_m - 1;
			printf("LOG: %c = %d, test %d\n", now_alpha ? 'a' : 'b', c_m, d_m), test_cnt++;
		}
		else {
			int d_l = d_m + 1, d_u = upper[c_m] + 1;
			printf("LOG: %c = %d, binary [%d, %d]\n", now_alpha ? 'a' : 'b', c_m, d_l, d_u), binary_cnt++;
			while (d_l < d_u) {
				int d_mid = (d_l + d_u) / 2;
				if (now_alpha) beta = d_mid;
				else alpha = d_mid;
				test();
				if (D.size == 0) d_u = d_mid;
				else d_l = d_mid + 1;
			}
			d_l--, upper[c_m] = d_l;
			
			if (now_alpha) beta = d_l; else alpha = d_l;
			get_subgraph();
			double A = 0.5 * (alpha + 1.0 * beta * V_num / U_num), B = 0.5 * (1.0 * alpha * U_num / V_num + beta);
			update_maxab(A, B), update_maxab(alpha, beta);
			if (now_alpha) {
				if (A < c_m) condition = 1;
				else condition = -1;
			}
			else {
				if (B < c_m) condition = 1;
				else condition = -1;
			}
			int a1, a2, a3, a4, b1, b2, b3, b4;
			a2 = A, b2 = B;
			a1 = a2 - 1, b1 = 1.0 * alpha * U_num / V_num + beta - 1.0 * U_num * a1 / V_num;
			a3 = a2 + 1, b3 = 1.0 * alpha * U_num / V_num + beta - 1.0 * U_num * a3 / V_num;
			a4 = a2 + 2, b4 = 1.0 * alpha * U_num / V_num + beta - 1.0 * U_num * a4 / V_num;
			printf("LOG: a1b1 = %d, a2b2 = %d, a3b3 = %d, a4b4 = %d\n", a1 * b1, a2 * b2, a3 * b3, a4 * b4);
			update_maxab(a1, b1), update_maxab(a3, b3), update_maxab(a4, b4);
		}
	}
	if (condition == 1) {
		divide(c_l, c_m, now_alpha);
		divide(c_m, c_u, now_alpha);
	}
	else {
		divide(c_m, c_u, now_alpha);
		divide(c_l, c_m, now_alpha);
	}
	return;
}
void Graph::get_maxab_divide() {
	Timer timer;
	line_cnt = test_cnt = binary_cnt = flow_cnt = 0;
	timer.start();
	get_pseudo();
	alpha = beta = pseudo;
	get_subgraph();
	timer.end();
	printf("- %-20s: %lf\n", "Get pseudoarboricity time", get_pseudo_time = timer.time());
	update_maxab(pseudo, pseudo);
	if (0.5 * (alpha + 1.0 * beta * V_num / U_num) <= pseudo) {
		get_maxab_dividea();
		get_maxab_divideb();
	}
	else {
		get_maxab_divideb();
		get_maxab_dividea();
	}
	return;
}
void Graph::get_maxab_dividea() {
	upper_alpha[pseudo + 1] = pseudo, upper_alpha[0] = dmax_V;
	divide(0, pseudo + 1, true);
}
void Graph::get_maxab_divideb() {
	upper_beta[pseudo + 1] = pseudo, upper_beta[0] = dmax_U;
	divide(0, pseudo + 1, false);
}
void Graph::output_core(Set<int>& C) {
	printf("- %-20s: %d\n", "Core size", C.size);

	if (true) {
		Set<int> output_nodes; output_nodes.alloc(N);
		sort(C.nodes, C.nodes + C.size);
		int i;
		for (i = 0; i < C.size && C.nodes[i] < U; i++)
			output_nodes.insert(C.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Core node in U");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Core node in U are too many, with size of", output_nodes.size);
		}

		output_nodes.clear();
		for (; i < C.size; i++)
			output_nodes.insert(C.nodes[i]);
		if (output_nodes.size <= OUTPUT_NODES_LIMIT) {
			printf("- %-20s: ", "Core node in V");
			for (int j = 0; j < output_nodes.size; j++)
				printf("%d ", output_nodes.nodes[j]);
			printf("\n");
		}
		else {
			printf("- %-20s: %d\n", "Core node in V are too many, with size of", output_nodes.size);
		}
	}
}
void Graph::output_dense_subgraph(char* dataset_address) {
	alpha = maxa, beta = maxb;
	get_subgraph();

	sort(D.nodes, D.nodes + D.size);
	int cnt_U = 0, cnt_V = 0, cnt_E = 0;
	for (int i = 0; i < D.size; i++) {
		int u = D.nodes[i];
		if (!in_U(u)) {
			cnt_V++;
			continue;
		}
		cnt_U++;
		for (int j = 0; j < undeg[u]; j++) {
			Edge& ne = e[adj[u][j]];
			int v = ne.v;
			if (!D.in[v]) continue;
			cnt_E++;
		}
	}
	printf("- %-20s: %d %d %d\n", "E, U, V of maxab subgraph", cnt_E, cnt_U, cnt_V);
	printf("- %-20s: %lf\n", "Density of maxab subgraph", cnt_E / (sqrt(cnt_U) * sqrt(cnt_V)));

	char dataset[256];

	const char* lastSlash = strrchr(dataset_address, '/');

	if (lastSlash != nullptr) {
		strncpy(dataset, lastSlash + 1, sizeof(dataset) - 1);
		dataset[sizeof(dataset) - 1] = '\0';
	}
	else {
		strncpy(dataset, dataset_address, sizeof(dataset) - 1);
		dataset[sizeof(dataset) - 1] = '\0';
	}

	const char* prefix = "./output/maxab_";
	memmove(dataset + strlen(prefix), dataset, strlen(dataset) + 1);
	memcpy(dataset, prefix, strlen(prefix));
	FILE* out = fopen(dataset, "w");
	check(out != NULL, "cannot open output file, maybe create a folder named \'output\'");

	fprintf(out, "%d %d\n", cnt_U, cnt_V);
	fprintf(out, "in U:\n");
	int i;
	for (i = 0; ; i++) {
		int x = D.nodes[i];
		if (in_V(x)) break;
		fprintf(out, "%d\n", x);
	}
	fprintf(out, "in V:\n");
	for (; i < D.size; i++) {
		int x = D.nodes[i];
		fprintf(out, "%d\n", x - U);
	}
	fclose(out);
}
void Graph::display() {
	for (int i = 0; i < M; i++) {
		Edge& ne = e[i];
		int from = in_U(ne.to) ? ne.v : ne.u;
		printf("%d %d\n", from, ne.to);
	}
}
void Graph::check_correctness() {
	return;
}

int main(int argc, char** argv) {
	if (argc != 3) {
	argument_error:
		printf("Usage: ./main <dataset_address> <algorithm> <core_reduction>\n");
		printf("algorithm:\n");
		printf("-a: get all max subgraphs\n");
		printf("-d: divide and conquer\n");
		return 0;
	}
	Timer timer;
	char dataset_address[1000]; strcpy(dataset_address, argv[1]);
	if (strcmp(argv[2], "-a") == 0) algorithm_used = ALL;
	else if (strcmp(argv[2], "-d") == 0) algorithm_used = DIVIDE;
	else goto argument_error;

	double runtime;
	printf("----------Now processing %s----------\n", dataset_address);

	timer.start();
	G.read_graph_from_dataset(dataset_address);
	timer.end(); runtime = timer.time();
	printf("- %-20s: %d, %d, %d\n", "|E|, |U|, |V|", G.M, G.U, G.V);
	printf("- %-20s: %lf\n", "Read graph time", runtime);

	timer.start();
	G.initialize_orientation();
	timer.end(); runtime = timer.time();
	printf("- %-20s: %lf\n", "Initialization time", runtime);

	if (algorithm_used == ALL) {
		timer.start();
		G.get_maxab_all();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf\n", "Get maxab time", runtime);
		printf("- %-20s: %d, %d, %d\n", "maxab = maxa * maxb", G.maxab, G.maxa, G.maxb);
	}
	else if (algorithm_used == DIVIDE) {
		timer.start();
		G.get_maxab_divide();
		timer.end(); runtime = timer.time();
		printf("- %-20s: %lf %d\n", "Get pseudoarboricity time", G.get_pseudo_time, G.pseudo);
		printf("- %-20s: %lf\n", "Get maxab time", runtime);
		printf("- %-20s: %d, %d, %d, %d\n", "line, test, binary, flow count", G.line_cnt, G.test_cnt, G.binary_cnt, G.flow_cnt);
		printf("- %-20s: %d, %d, %d\n", "maxab = maxa * maxb", G.maxab, G.maxa, G.maxb);
		G.output_dense_subgraph(dataset_address);
	}
	return 0;
}

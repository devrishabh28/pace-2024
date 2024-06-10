#include <iostream>
#include <signal.h>
#include <unistd.h>
#include <cstring>
#include <vector>
#include <cmath>

#include <fstream>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

class SegmentTree
{
private:
    std::vector<int> tree;
    int size;

    void build(const std::vector<int> &arr, int v, int tl, int tr)
    {
        if (tl == tr)
        {
            tree[v] = arr[tl];
        }
        else
        {
            int tm = (tl + tr) / 2;
            build(arr, v * 2, tl, tm);
            build(arr, v * 2 + 1, tm + 1, tr);
            tree[v] = tree[v * 2] + tree[v * 2 + 1];
        }
    }

    int query(int v, int tl, int tr, int l, int r)
    {
        if (l > r)
            return 0;
        if (l == tl && r == tr)
            return tree[v];
        int tm = (tl + tr) / 2;
        return query(v * 2, tl, tm, l, std::min(r, tm)) +
               query(v * 2 + 1, tm + 1, tr, std::max(l, tm + 1), r);
    }

    void update(int v, int tl, int tr, int pos, int new_val)
    {
        if (tl == tr)
        {
            tree[v] = new_val;
        }
        else
        {
            int tm = (tl + tr) / 2;
            if (pos <= tm)
                update(v * 2, tl, tm, pos, new_val);
            else
                update(v * 2 + 1, tm + 1, tr, pos, new_val);
            tree[v] = tree[v * 2] + tree[v * 2 + 1];
        }
    }

public:
    SegmentTree(const std::vector<int> &arr)
    {
        size = arr.size();
        tree.resize(4 * size);
        build(arr, 1, 0, size - 1);
    }

    int query(int l, int r)
    {
        return query(1, 0, size - 1, l, r);
    }

    void update(int pos, int new_val)
    {
        update(1, 0, size - 1, pos, new_val);
    }
};

class PaceGraph
{
private:
public:
    vector<int> right;
    map<int, int> right_order;
    int left;

    vector<pair<int, int>> edgeset;
    map<int, vector<int>> edges;

    map<pair<int, int>, unsigned long long> crossing_table;

    int preprocess_vertex;

    PaceGraph(int a, int b, vector<pair<int, int>> edges = {}, vector<int> order = {});
    vector<int> neighbours(int u);
    void set_order(vector<int> order = {}, bool check = false);
    bool cross(int a, int b, int c, int d);
    void swap(int a, int b);
    unsigned long long countcrossings_segtree(volatile sig_atomic_t *tle = NULL);
    unsigned long long crossing_number(int u, int v, volatile sig_atomic_t *tle = NULL);
    unsigned long long count_crossings_partition(int low, int high, volatile sig_atomic_t *tle = NULL);
    void set_preprocess_vertex(int zero_vertex) { preprocess_vertex = zero_vertex; }

    // Getter for right side
    const vector<int> &get_right() const
    {
        return right;
    }

    // Setter for right side
    void set_right(const vector<int> &new_right)
    {
        right = new_right;
        right_order.clear();
        for (size_t p = 0; p < right.size(); ++p)
            right_order[right[p]] = p;
    }

    // Getter for left side
    int get_left() const
    {
        return left;
    }

    // Setter for left side
    void set_left(int new_left)
    {
        left = new_left;
    }

    // Getter for edges
    const map<int, vector<int>> &get_edges() const
    {
        return edges;
    }

    // Setter for edges
    void set_edges(const map<int, vector<int>> &new_edges)
    {
        edges = new_edges;
        edgeset.clear();
        for (const auto &pair : edges)
        {
            int u = pair.first;
            for (int v : pair.second)
            {
                edgeset.push_back({u, v});
            }
        }
    }

    static PaceGraph from_gr(ifstream &gr, vector<int> order = {})
    {
        int a = 0, b = 0;
        bool pfound = false;
        vector<pair<int, int>> edges;
        string line;
        while (getline(gr, line))
        {
            if (line[0] == 'p')
            {
                sscanf(line.c_str(), "p ocr %d %d", &a, &b);
                pfound = true;
            }
            else if (line[0] == 'c')
            {
                continue;
            }
            else if (pfound)
            {
                int u, v;
                sscanf(line.c_str(), "%d %d", &u, &v);
                edges.push_back({u, v});
            }
            else
            {
                throw invalid_argument("ERROR: Encountered edge before p-line.");
            }
        }
        return PaceGraph(a, b, edges, order);
    }

    static PaceGraph from_gr(vector<int> order = {})
    {
        int a = 0, b = 0, c;
        bool pfound = false;
        vector<pair<int, int>> edges;
        string line;

        while (getline(cin, line))
        {
            if (line[0] == 'p')
            {
                sscanf(line.c_str(), "p ocr %d %d %d", &a, &b, &c);
                pfound = true;
            }
            else if (line[0] == 'c')
            {
                continue;
            }
            else if (pfound)
            {
                int u, v;
                sscanf(line.c_str(), "%d %d", &u, &v);
                edges.push_back({u, v});
            }
            else
            {
                throw invalid_argument("ERROR: Encountered edge before p-line.");
            }

            if (edges.size() >= c)
                break;
        }
        return PaceGraph(a, b, edges, order);
    }

    string to_gr() const
    {
        string newline = "\n";
        string thestring = "p ocr " + to_string(left) + " " + to_string(right.size()) + " " + to_string(edgeset.size()) + "\n";
        for (const auto &e : edgeset)
            thestring += to_string(e.first) + " " + to_string(e.second) + newline;
        return thestring;
    }
};

PaceGraph::PaceGraph(int a, int b, vector<pair<int, int>> edges, vector<int> order)
{

    preprocess_vertex = 0;

    if (!order.empty())
    {
        set<int> oelements(order.begin(), order.end());
        set<int> relements;
        for (int i = a + 1; i <= a + b; ++i)
            relements.insert(i);

        set<int> diff;
        set_difference(oelements.begin(), oelements.end(), relements.begin(), relements.end(), inserter(diff, diff.begin()));

        if (!diff.empty())
            throw invalid_argument("Error! 'order' contains elements not present in the graph.");

        right = order;
    }
    else
    {
        for (int i = a + 1; i <= a + b; ++i)
            right.push_back(i);
    }

    for (size_t p = 0; p < right.size(); ++p)
        right_order[right[p]] = p;

    left = a;

    if (!edges.empty())
    {
        sort(edges.begin(), edges.end());
        for (const auto &e : edges)
        {
            int u = (e.first < e.second) ? e.first : e.second;
            int v = (e.second < e.first) ? e.first : e.second;

            this->edges[u].push_back(v);
            this->edges[v].push_back(u);
            edgeset.push_back({u, v});
        }
    }
}

vector<int> PaceGraph::neighbours(int u)
{
    return edges.at(u);
}

void PaceGraph::set_order(vector<int> order, bool check)
{
    if (check)
    {
        set<int> oelements(order.begin(), order.end());
        set<int> relements;
        for (size_t i = 0; i < right.size(); ++i)
            relements.insert(i + left + 1);

        set<int> diff;
        set_difference(oelements.begin(), oelements.end(), relements.begin(), relements.end(), inserter(diff, diff.begin()));

        if (!diff.empty())
            throw invalid_argument("Error! 'order' contains elements not present in the graph.");

        right = order;
    }
    else
    {
        right = order;
    }

    for (size_t p = 0; p < right.size(); ++p)
        right_order[right[p]] = p;
}

bool PaceGraph::cross(int a, int b, int c, int d)
{
    if (a > b)
        swap(a, b);

    if (c > d)
        swap(c, d);

    b = right_order.at(b);
    d = right_order.at(d);

    return ((a < c && b > d) || (c < a && d > b));
}

void PaceGraph::swap(int a, int b)
{
    int apos = right_order.at(a);
    int bpos = right_order.at(b);

    std::swap(right[apos], right[bpos]);
    right_order[a] = bpos;
    right_order[b] = apos;
}

unsigned long long PaceGraph::countcrossings_segtree(volatile sig_atomic_t *tle)
{
    int size = right.size();
    vector<pair<int, int>> edges(edgeset);
    sort(edges.begin(), edges.end(), [&](const auto &a, const auto &b)
         {
            if (a.first == b.first)
                return right_order.at(a.second) < right_order.at(b.second);
            return a.first < b.first; });

    if (tle && *tle)
        return 0;

    vector<int> arr(size + 1, 0);
    SegmentTree t(arr);

    unsigned long long crossings = 0;
    for (const auto &edge : edges)
    {
        int old = t.query(right_order.at(edge.second), right_order.at(edge.second));
        t.update(right_order.at(edge.second), old + 1);

        if (tle && *tle)
            return 0;

        if (size > right_order.at(edge.second) + 1)
        {
            int crossings_found = t.query(right_order.at(edge.second) + 1, size);
            crossings += crossings_found;
        }
    }
    return crossings;
}

unsigned long long PaceGraph::crossing_number(int u, int v, volatile sig_atomic_t *tle)
{
    if (edges.find(u) == edges.end() || edges.find(v) == edges.end())
        return 0;
    if (u == v)
        return 0;

    const vector<int> &adj_u = edges.at(u);
    const vector<int> &adj_v = edges.at(v);

    unsigned long long cr = 0;

    for (int i = adj_u.size() - 1; i >= 0; i--)
    {
        int w = adj_u[i];
        auto it = lower_bound(adj_v.begin(), adj_v.end(), w);

        if (tle && *tle)
            return 0;

        cr += distance(adj_v.begin(), it);
        if (it == adj_v.begin())
            break;
    }

    return cr;
}

unsigned long long PaceGraph::count_crossings_partition(int low, int high, volatile sig_atomic_t *tle)
{
    unsigned long long total_crossings = 0;
    unsigned long long crossings = 0;

    for (int i = low; i < high; ++i)
    {
        for (int j = i + 1; j <= high; ++j)
        {
            pair<int, int> key(right[i] - left, right[j] - left);
            if (crossing_table.find(key) == crossing_table.end())
            {
                crossings = crossing_number(right[i], right[j], tle);
                if (tle && *tle)
                    break;
                crossing_table[key] = crossings;
            }
            if (tle && *tle)
                break;
            total_crossings += crossing_table[key];
        }
    }
    return total_crossings;
}

void write_graph(const PaceGraph &graph, const string &filename)
{
    ofstream f(filename);
    f << graph.to_gr();
}

void write_solution(const vector<int> &order, const string &filename)
{
    ofstream f(filename);
    for (int n : order)
        f << n << endl;
}

PaceGraph read_graph_from_file(const string &filename, const vector<int> &order = {})
{
    ifstream f(filename);
    return PaceGraph::from_gr(f, order);
}

vector<int> read_solution(const string &filename)
{
    vector<int> order;
    ifstream f(filename);
    int n;
    while (f >> n)
        order.push_back(n);
    return order;
}

/*--------------------------------------- Methods --------------------------------------- */
vector<int> median(const PaceGraph &graph, volatile sig_atomic_t *tle = NULL);
vector<int> barycenter(const PaceGraph &graph, volatile sig_atomic_t *tle = NULL);
vector<int> local_search(PaceGraph &graph, bool multi_mode = false, volatile sig_atomic_t *tle = NULL, unsigned long epochs = 0, vector<int> order = {});
vector<int> modified_local_search(PaceGraph &graph, bool multi_mode = false, volatile sig_atomic_t *tle = NULL, unsigned long epochs = 0, vector<int> order = {});
vector<int> backtracking_local_search(PaceGraph &graph, bool multi_mode = false, volatile sig_atomic_t *tle = NULL, unsigned long epochs = 0, vector<int> order = {});
// vector<int> local_search_multi_swap(PaceGraph &graph, int epochs = 0, vector<int> order = {});
/*--------------------------------------- Methods --------------------------------------- */

/*--------------------------------------- Helpers --------------------------------------- */
double median_helper(vector<int> nums, volatile sig_atomic_t *tle = NULL);
double barycenter_helper(vector<int> nums, volatile sig_atomic_t *tle = NULL);
void get_best_initial_order(PaceGraph &graph, volatile sig_atomic_t *tle = NULL);
vector<int> backtracking_local_search_helper(PaceGraph &graph, bool &changes, int separation = 1, bool multi_mode = false, volatile sig_atomic_t *tle = NULL, unsigned long epochs = 0, vector<int> order = {});

// Custom comparator function to sort based on all elements of the inner vectors
auto comparator = [](const std::vector<double> &a, const std::vector<double> &b)
{
    return lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
};
/*--------------------------------------- Helpers --------------------------------------- */

/*------------------------------- Methods Implementation -------------------------------- */
vector<int> median(PaceGraph &graph, volatile sig_atomic_t *tle)
{
    vector<vector<double>> positions;

    for (auto &&u : graph.right)
    {
        try
        {
            positions.push_back({median_helper(graph.edges.at(u), tle), (double)(graph.edges.at(u).size() % 2 == 0), (double)u});
        }
        catch (exception &e)
        {
            positions.push_back({0, 1, (double)u});
        }
    }

    sort(positions.begin(), positions.end(), comparator);

    int zero_vertex = 0;
    vector<int> order;
    for (auto &&i : positions)
    {
        order.push_back(i[2]);
        if (i[0] == 0)
            zero_vertex++;
    }

    /* Preprocessing */
    graph.set_preprocess_vertex(zero_vertex);

    return order;
}

vector<int> barycenter(PaceGraph &graph, volatile sig_atomic_t *tle)
{
    vector<vector<double>> positions;

    for (auto &&u : graph.right)
    {
        try
        {
            positions.push_back({barycenter_helper(graph.edges.at(u), tle), (double)u});
        }
        catch (exception &e)
        {
            positions.push_back({0, (double)u});
        }
    }

    sort(positions.begin(), positions.end(), comparator);

    int zero_vertex = 0;
    vector<int> order;
    for (auto &&i : positions)
    {
        order.push_back(i[1]);
        if (i[0] == 0)
            zero_vertex++;
    }

    /* Preprocessing */
    graph.set_preprocess_vertex(zero_vertex);

    return order;
}

vector<int> local_search(PaceGraph &graph, bool multi_mode, volatile sig_atomic_t *tle, unsigned long epochs, vector<int> order)
{

    if (order.empty())
    {
        get_best_initial_order(graph, tle);
        order = graph.right;
    }
    else
    {
        graph.set_order(order);
        order = graph.right;
    }

    if (epochs <= 0)
    {
        epochs = order.size() * order.size();
    }

    set<pair<int, int>> visited_pairs;
    pair<int, int> suited_pair;
    bool swapped = false;

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        if (tle && *tle)
            break;

        swapped = false;
        for (int i = graph.preprocess_vertex; i < order.size() - 1; i++)
        {
            if (tle && *tle)
                break;

            if (order[i] <= order[i + 1])
            {
                suited_pair = {order[i], order[i + 1]};
            }
            else
            {
                suited_pair = {order[i + 1], order[i]};
            }

            if (visited_pairs.find(suited_pair) == visited_pairs.end())
            {
                if (graph.crossing_number(order[i], order[i + 1], tle) > graph.crossing_number(order[i + 1], order[i], tle))
                {
                    visited_pairs.insert(suited_pair);
                    swap(order[i], order[i + 1]);
                    swapped = true;
                    if (!multi_mode)
                        break;
                }
            }
        }

        if (!swapped)
            break;
    }

    if (tle && *tle)
        return order;

    graph.set_order(order);
    return order;
}

vector<int> modified_local_search(PaceGraph &graph, bool multi_mode, volatile sig_atomic_t *tle, unsigned long epochs, vector<int> order)
{

    if (order.empty())
    {
        get_best_initial_order(graph, tle);
        order = graph.right;
    }
    else
    {
        graph.set_order(order);
        order = graph.right;
    }

    if (epochs <= 0)
    {
        epochs = order.size() * order.size();
    }

    set<pair<int, int>> visited_pairs;
    pair<int, int> suited_pair;
    bool swapped = false;
    unsigned long long crossings_uv;
    unsigned long long crossings_vu;

    for (int separation = 1; separation < order.size(); separation++)
    {
        for (int epoch = 0; epoch < epochs; epoch++)
        {
            if (tle && *tle)
                break;

            swapped = false;
            for (int i = graph.preprocess_vertex; i < order.size() - separation; i++)
            {
                if (tle && *tle)
                    break;

                if (order[i] <= order[i + separation])
                {
                    suited_pair = {order[i], order[i + separation]};
                }
                else
                {
                    suited_pair = {order[i + separation], order[i]};
                }

                if (visited_pairs.find(suited_pair) == visited_pairs.end())
                {

                    crossings_uv = graph.count_crossings_partition(i, i + separation, tle);

                    if (tle && *tle)
                        break;

                    graph.swap(order[i], order[i + separation]);
                    crossings_vu = graph.count_crossings_partition(i, i + separation, tle);

                    if (crossings_uv > crossings_vu)
                    {
                        visited_pairs.insert(suited_pair);
                        swapped = true;
                        if (!multi_mode)
                            break;
                    }
                    else
                    {
                        graph.swap(order[i], order[i + separation]);
                    }

                    order = graph.right;
                }
            }

            // cout << "Epoch: " << epoch << ", Sep: " << separation <<  ", cr:" << graph.countcrossings_segtree() << endl;

            if (!swapped)
                break;
        }
        if (tle && *tle)
            break;
    }

    if (tle && *tle)
        return order;

    graph.set_order(order);
    return order;
}

vector<int> backtracking_local_search(PaceGraph &graph, bool multi_mode, volatile sig_atomic_t *tle, unsigned long epochs, vector<int> order)
{

    if (order.empty())
    {
        get_best_initial_order(graph, tle);
        order = graph.right;
    }
    else
    {
        graph.set_order(order);
        order = graph.right;
    }

    if (epochs <= 0)
    {
        epochs = order.size() * order.size();
    }

    bool changes = false;

    for (int separation = 1; separation < order.size(); separation++)
    {
        changes = false;
        order = backtracking_local_search_helper(graph, changes, separation, multi_mode, tle, epochs, order);

        if (changes)
            separation -= 2;

        if (tle && *tle)
            break;
    }

    if (tle && *tle)
        return order;

    graph.set_order(order);
    return order;
}

/* vector<int> local_search_multi_swap(PaceGraph &graph, int epochs, vector<int> order)
{

    if (order.empty())
    {
        get_best_initial_order(graph);
        order = graph.right;
    }

    if (epochs <= 0)
    {
        epochs = order.size() * order.size();
    }

    set<pair<int, int>> visited_pairs;
    pair<int, int> suited_pair;

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        for (int i = 0; i < order.size() - 1; i++)
        {
            if (order[i] <= order[i + 1])
            {
                suited_pair = {order[i], order[i + 1]};
            }
            else
            {
                suited_pair = {order[i + 1], order[i]};
            }

            if (visited_pairs.find(suited_pair) == visited_pairs.end())
            {
                if (graph.crossing_number(order[i], order[i + 1]) > graph.crossing_number(order[i + 1], order[i]))
                {
                    visited_pairs.insert(suited_pair);
                    swap(order[i], order[i + 1]);
                    swapped = true;
                }
            }
        }
    }

    return order;
} */

/*------------------------------- Methods Implementation -------------------------------- */

/*------------------------------- Helpers Implementation -------------------------------- */

double median_helper(vector<int> nums, volatile sig_atomic_t *tle)
{
    size_t size = nums.size();
    nth_element(nums.begin(), nums.begin() + size / 2, nums.end());
    double median;
    if (size % 2 == 0)
    {
        int val1 = nums[size / 2 - 1];
        int val2 = nums[size / 2];
        median = (val1 + val2) / 2.0;
    }
    else
    {
        median = nums[size / 2];
    }

    return median;
}

double barycenter_helper(vector<int> nums, volatile sig_atomic_t *tle)
{
    if (nums.empty())
    {
        return 0;
    }

    unsigned int sum = 0;
    for (int num : nums)
    {
        sum += num;
    }

    return (double)(sum) / nums.size();
}

void get_best_initial_order(PaceGraph &graph, volatile sig_atomic_t *tle)
{

    unsigned long long best_crossing_score = graph.countcrossings_segtree(tle);
    vector<int> best_order = graph.right;

    // Median heuristic
    vector<int> current_order = median(graph, tle);
    graph.set_order(current_order);
    unsigned long long current_crossing_score = graph.countcrossings_segtree(tle);

    if (current_crossing_score <= best_crossing_score)
    {
        best_crossing_score = current_crossing_score;
        best_order = current_order;
    }
    else
    {
        graph.set_order(best_order);
    }

    // Barycenter heuristic
    current_order = barycenter(graph, tle);
    graph.set_order(current_order);
    current_crossing_score = graph.countcrossings_segtree(tle);

    if (current_crossing_score <= best_crossing_score)
    {
        best_crossing_score = current_crossing_score;
        best_order = current_order;
    }
    else
    {
        graph.set_order(best_order);
    }
}

vector<int> backtracking_local_search_helper(PaceGraph &graph, bool &changes, int separation, bool multi_mode, volatile sig_atomic_t *tle, unsigned long epochs, vector<int> order)
{

    if (order.empty())
    {
        get_best_initial_order(graph, tle);
        order = graph.right;
    }
    else
    {
        graph.set_order(order);
        order = graph.right;
    }

    if (epochs <= 0)
    {
        epochs = order.size() * order.size();
    }

    set<pair<int, int>> visited_pairs;
    pair<int, int> suited_pair;
    bool swapped = false;
    unsigned long long crossings_uv;
    unsigned long long crossings_vu;

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        if (tle && *tle)
            break;

        swapped = false;
        for (int i = graph.preprocess_vertex; i < order.size() - separation; i++)
        {
            if (tle && *tle)
                break;

            if (order[i] <= order[i + separation])
            {
                suited_pair = {order[i], order[i + separation]};
            }
            else
            {
                suited_pair = {order[i + separation], order[i]};
            }

            if (visited_pairs.find(suited_pair) == visited_pairs.end())
            {
                crossings_uv = graph.count_crossings_partition(i, i + separation, tle);
                if (tle && *tle)
                    break;
                graph.swap(order[i], order[i + separation]);
                crossings_vu = graph.count_crossings_partition(i, i + separation, tle);
                if (crossings_uv > crossings_vu)
                {
                    visited_pairs.insert(suited_pair);
                    swapped = true;
                    changes = true;
                    if (!multi_mode)
                        break;
                }
                else
                {
                    graph.swap(order[i], order[i + separation]);
                }
                order = graph.right;
            }
        }

        // cout << "Epoch: " << epoch << ", Sep: " << separation <<  ", cr:" << graph.countcrossings_segtree() << endl;
        if (!swapped)
            break;
    }

    if (tle && *tle)
        return order;

    graph.set_order(order);
    return order;
}

/*------------------------------- Helpers Implementation -------------------------------- */

volatile sig_atomic_t tle = 0;
int n;
float worker;
vector<int> order;

void term(int signum)
{
    tle = 1;
}

void print_solution()
{
    for (auto &&i : order)
    {
        printf("%d\n", i);
    }
    exit(0);
}

void test()
{
    PaceGraph gr = read_graph_from_file("./dataset/10.gr");
    cout << "Crossings: " << gr.countcrossings_segtree(&tle) << endl;

    // get_best_initial_order(gr, &tle);
    // cout << "Crossings: " << gr.countcrossings_segtree(&tle) << endl;

    order = backtracking_local_search(gr, false, &tle);
    cout << "Crossings: " << gr.countcrossings_segtree() << endl;
}

int main()
{
    ios::sync_with_stdio(0);
    cin.tie(0);

    // to handle SIGTERM
    struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);
    sigaction(SIGINT, &action, NULL);

    // test();

    PaceGraph gr = PaceGraph::from_gr();
    order = backtracking_local_search(gr, false, &tle);
    // cout << "Crossings: " << gr.countcrossings_segtree() << endl;

    print_solution();

    return 0;
}
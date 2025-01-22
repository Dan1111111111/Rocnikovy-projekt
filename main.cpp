#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <algorithm>
#include <iomanip>
#include <cstdlib>

using namespace std;

static const long double EPS = 1e-12;
static const long double INF = 1e18L;

struct Circle {
    long double x, y, r;
};

struct Point {
    long double x, y;
    int circle_index;
};

struct Edge {
    int index;
    long double length;
};

long double distance(const Point &a, const Point &b)
{
    return hypotl(a.x - b.x, a.y - b.y);
}

bool inside_circle(long double px, long double py, const Circle &C)
{
    long double dx = px - C.x;
    long double dy = py - C.y;

    return (dx * dx + dy * dy < C.r * C.r - EPS);
}

bool inside_any_circle(const Point &p, const vector<Circle> &circles)
{
    for (size_t i = 0; i < circles.size(); i++)
    {
        if (inside_circle(p.x, p.y, circles[i]))
            return true;
    }
    return false;
}

bool segment_intersects_circle(const Point &p1, const Point &p2, const Circle &C)
{
    long double dx = p2.x - p1.x;
    long double dy = p2.y - p1.y;
    long double seg_length = dx * dx + dy * dy;

    if (seg_length < EPS)
        return inside_circle(p1.x, p1.y, C);

    long double t = ((C.x - p1.x) * dx + (C.y - p1.y) * dy) / seg_length;
    t = max(0.0L, min(1.0L, t));
    long double cx = p1.x + t * dx;
    long double cy = p1.y + t * dy;

    return inside_circle(cx, cy, C);
}

bool can_use_segment(const Point &p1, const Point &p2, const vector<Circle> &circles)
{
    for (size_t i = 0; i < circles.size(); i++)
    {
        if (segment_intersects_circle(p1, p2, circles[i]))
            return false;
    }
    return true;
}

vector<Point> find_tangent_points(const Point &P, int circle_id, const Circle &C)
{
    vector<Point> answer;
    long double dx = C.x - P.x;
    long double dy = C.y - P.y;
    long double d2 = dx * dx + dy * dy;
    long double r = C.r;

    if (d2 < r * r + EPS)
        return answer;

    long double d = sqrtl(d2);
    long double beta = acosl(r/d);
    long double gamma = atan2l(dy, dx);

    int signs[2] = { -1, 1 };
    for (int i = 0; i < 2; i++)
    {
        int s = signs[i];
        long double final_angle = gamma + s * beta;
        Point T;
        T.x = C.x - r * cosl(final_angle);
        T.y = C.y - r * sinl(final_angle);
        T.circle_index = circle_id;
        answer.push_back(T);
    }
    return answer;
}

long double adjust_angle(long double a)
{
    while (a < 0)
        a += 2.0L * M_PI;
    while (a >= 2.0L * M_PI)
        a -= 2.0L * M_PI;
    return a;
}

bool can_use_arc(const Point &circle_center, long double R, long double start_angle, long double end_angle, bool ccw, const vector<Circle> &circles)
{
    long double arc = 0.0L;
    if (ccw)
    {
        if (end_angle < start_angle)
            end_angle += 2.0L * M_PI;

        arc = end_angle - start_angle;
    }
    else
    {
        if (start_angle < end_angle)
            start_angle += 2.0L * M_PI;

        arc = start_angle - end_angle;
    }

    const int arc_sample_count = 15;
    for (int i = 0; i <= arc_sample_count; i++)
    {
        long double fraction = (long double)i / arc_sample_count;
        long double direction;
        if (ccw)
            direction = 1.0L;
        else
            direction = -1.0L;

        long double alpha = start_angle + fraction * (direction * arc);
        long double px = circle_center.x + R * cosl(alpha);
        long double py = circle_center.y + R * sinl(alpha);

        for (size_t j = 0; j < circles.size(); j++)
        {
            long double dx = px - circles[j].x;
            long double dy = py - circles[j].y;
            long double d2 = dx * dx + dy * dy;

            if (d2 < circles[j].r * circles[j].r - EPS)
                return false;
        }
    }
    return true;
}

void add_edge(int u, int v, long double length, vector<vector<Edge> > &graph)
{
    Edge e1;
    e1.index = v;
    e1.length = length;
    graph[u].push_back(e1);

    Edge e2;
    e2.index = u;
    e2.length = length;
    graph[v].push_back(e2);
}

void run_unit_tests();

int main(int argc, char** argv)
{
    for (int i = 1; i < argc; i++)
    {
        string arg = argv[i];
        if (arg == "-test")
        {
            run_unit_tests();
            return 0;
        }
    }

    cin.tie(0)->sync_with_stdio(0);

    int N;
    Point A, B;
    cin >> A.x >> A.y >> B.x >> B.y;
    cin >> N;

    A.circle_index = -1;
    B.circle_index = -1;

    vector<Circle> circles(N);
    for (int i = 0; i < N; i++)
        cin >> circles[i].x >> circles[i].y >> circles[i].r;

    if (inside_any_circle(A, circles) || inside_any_circle(B, circles))
    {
        cout << "Path between A and B does not exist.\n";
        return 0;
    }

    vector<Point> vertices;
    vertices.push_back(A);
    vertices.push_back(B);

    for (int i = 0; i < N; i++)
    {
        long double distance_A = hypotl(A.x - circles[i].x, A.y - circles[i].y);
        long double distance_B = hypotl(B.x - circles[i].x, B.y - circles[i].y);

        if (fabsl(distance_A - circles[i].r) < EPS)
        {
            Point A_copy = A;
            A_copy.circle_index = i;
            vertices.push_back(A_copy);
        }
        if (fabsl(distance_B - circles[i].r) < EPS)
        {
            Point B_copy = B;
            B_copy.circle_index = i;
            vertices.push_back(B_copy);
        }
    }

    for (int i = 0; i < N; i++)
    {
        vector<Point> tangent_points_A = find_tangent_points(A, i, circles[i]);
        for (size_t j = 0; j < tangent_points_A.size(); j++)
            vertices.push_back(tangent_points_A[j]);

        vector<Point> tangent_points_B = find_tangent_points(B, i, circles[i]);
        for (size_t j = 0; j < tangent_points_B.size(); j++)
            vertices.push_back(tangent_points_B[j]);
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            long double dx = circles[j].x - circles[i].x;
            long double dy = circles[j].y - circles[i].y;
            long double d = hypotl(dx, dy);
            long double r1 = circles[i].r;
            long double r2 = circles[j].r;

            if (d < fabsl(r1 - r2) + EPS)
                continue;

            bool swapped = false;
            Circle C_big = circles[i];
            Circle C_small = circles[j];
            int index_Big = i;
            int index_Small = j;
            long double R = r1 - r2;

            if (r1 < r2)
            {
                swapped = true;
                swap(C_big, C_small);
                swap(index_Big, index_Small);
                R = r2 - r1;
            }

            long double d2 = (C_small.x - C_big.x) * (C_small.x - C_big.x) + (C_small.y - C_big.y) * (C_small.y - C_big.y);
            if (d2 < R * R + EPS)
                continue;

            long double d_ = sqrtl(d2);
            long double beta = acosl(R / d_);
            long double gamma = atan2l(C_small.y - C_big.y, C_small.x - C_big.x);

            int signs[2] = {-1, 1};
            for (int k = 0; k < 2; k++)
            {
                int s = signs[k];
                long double final_angle = gamma + s * beta;
                long double direction_x = cosl(final_angle);
                long double direction_y = sinl(final_angle);

                Point T_big, T_small;
                T_big.x = C_big.x + C_big.r * direction_x;
                T_big.y = C_big.y + C_big.r * direction_y;
                T_big.circle_index = index_Big;
                T_small.x = C_small.x + C_small.r * direction_x;
                T_small.y = C_small.y + C_small.r * direction_y;
                T_small.circle_index = index_Small;

                if (swapped)
                    swap(T_big, T_small);

                vertices.push_back(T_big);
                vertices.push_back(T_small);
            }
        }
    }

    int num_vertices = vertices.size();
    vector<vector<Edge> > graph(num_vertices);
    for (int i = 0; i < num_vertices; i++)
    {
        for (int j = i + 1; j < num_vertices; j++)
        {
            if (can_use_segment(vertices[i], vertices[j], circles))
            {
                long double d = distance(vertices[i], vertices[j]);
                add_edge(i, j, d, graph);
            }
        }
    }

    for (int i = 0; i < num_vertices; i++)
    {
        for (int j = i + 1; j < num_vertices; j++)
        {
            int current_circle_index = vertices[i].circle_index;

            if (current_circle_index < 0 || current_circle_index != vertices[j].circle_index)
                continue;

            const Circle &C = circles[current_circle_index];
            long double Ax = vertices[i].x - C.x;
            long double Ay = vertices[i].y - C.y;
            long double Bx = vertices[j].x - C.x;
            long double By = vertices[j].y - C.y;
            long double angle_A = adjust_angle(atan2l(Ay, Ax));
            long double angle_B = adjust_angle(atan2l(By, Bx));

            long double adjusted_angel_B = angle_B;
            if (angle_B < angle_A)
                adjusted_angel_B += 2.0L * M_PI;
            long double arc_length_CCW = (adjusted_angel_B - angle_A) * C.r;

            long double adjusted_angel_A = angle_A;
            if (angle_A < angle_B)
                adjusted_angel_A += 2.0L * M_PI;
            long double arc_length_CW = (adjusted_angel_A - angle_B) * C.r;

            bool arc_CCW = can_use_arc(vertices[i], C.r, angle_A, angle_B, true, circles);
            bool arc_CW  = can_use_arc(vertices[i], C.r, angle_A, angle_B, false, circles);

            long double shortest_arc_lenght;
            if (arc_CCW && (!arc_CW || arc_length_CCW <= arc_length_CW))
                shortest_arc_lenght = arc_length_CCW;
            else if (arc_CW)
                shortest_arc_lenght = arc_length_CW;
            else
                continue;

            add_edge(i, j, shortest_arc_lenght, graph);
        }
    }

    vector<long double> shortest_distance(num_vertices, INF);
    vector<bool> done(num_vertices, false);
    priority_queue<pair<long double,int>, vector<pair<long double,int> >, greater<pair<long double,int> > > pq;
    shortest_distance[0] = 0.0L;
    pq.push(pair<long double,int>(0.0L, 0));

    while (!pq.empty())
    {
        pair<long double,int> current_pair = pq.top();
        pq.pop();
        long double current_distance = current_pair.first;
        int current_vertex = current_pair.second;

        if (done[current_vertex])
            continue;
        done[current_vertex] = true;

        if (current_vertex == 1)
            break;

        for (size_t i = 0; i < graph[current_vertex].size(); i++)
        {
            int neighbor = graph[current_vertex][i].index;
            long double new_distance = current_distance + graph[current_vertex][i].length;

            if (new_distance + EPS < shortest_distance[neighbor])
            {
                shortest_distance[neighbor] = new_distance;
                pq.push(pair<long double, int>(new_distance, neighbor));
            }
        }
    }

    if (shortest_distance[1] > INF / 2)
        cout << "Path between A and B does not exist.\n";
    else
        cout << "The shortest path between A and B is -➜ " << fixed << setprecision(5) << shortest_distance[1] << "\n";

    return 0;
}

// UNIT TESTY:

void error(const string &msg)
{
    cout << "[ERROR] " << msg << "\n";
}

#define ASSERT_NEAR(val, expected, eps) \
do                                      \
{                           \
    long double diff = fabsl((long double)(val) - (long double)(expected));   \
    if (diff > (eps))                                                  \
        error("Value " + to_string((long double)(val)) + " != " + to_string((long double)(expected)) + " (diff=" + to_string((long double)(diff)) + "), EPS=" + #eps);    \
}                                       \
while(0)

#define ASSERT_TRUE(cond) \
do                        \
{             \
    if (!(cond))                         \
        error("Condition failed: " #cond); \
}                         \
while(0)

#define ASSERT_FALSE(cond) \
do                         \
{            \
    if ((cond))                         \
        error("Condition failed (expected false): " #cond);\
}                          \
while(0)

void test_header(const string &testName)
{
    cout << "⇨" << testName << "\n";
}

void test_distance()
{
    test_header("test_distance");
    Point p1;
    p1.x = 0.0L; p1.y = 0.0L; p1.circle_index = -1;
    Point p2;
    p2.x = 3.0L; p2.y = 4.0L; p2.circle_index = -1;

    long double dist = distance(p1, p2);
    ASSERT_NEAR(dist, 5.0L, 1e-15);

    ASSERT_NEAR(distance(p1, p1), 0.0L, 1e-15);
}

void test_inside_circle()
{
    test_header("test_inside_circle");
    Circle c1;
    c1.x = 0.0L; c1.y = 0.0L; c1.r = 2.0L;

    ASSERT_TRUE( inside_circle(1.0L, 1.0L, c1) );

    ASSERT_FALSE( inside_circle(2.0L, 0.0L, c1) );

    ASSERT_FALSE( inside_circle(3.0L, 0.0L, c1) );
}

void test_segment_intersects_circle()
{
    test_header("test_segment_intersects_circle");
    Circle circle2;
    circle2.x = 0.0L; circle2.y = 0.0L; circle2.r = 1.0L;

    Point a; a.x = 2.0L; a.y = 0.0L; a.circle_index = -1;
    Point b; b.x = -2.0L; b.y = 0.0L; b.circle_index = -1;
    ASSERT_TRUE( segment_intersects_circle(a, b, circle2) );

    Point c; c.x = 2.0L; c.y = 2.0L; c.circle_index = -1;
    Point d; d.x = -2.0L; d.y = 2.0L; d.circle_index = -1;
    ASSERT_FALSE( segment_intersects_circle(c, d, circle2) );

    Point e; e.x = 1.1L; e.y = 1.1L; e.circle_index = -1;
    Point f; f.x = 2.0L; f.y = 2.0L; f.circle_index = -1;
    ASSERT_FALSE( segment_intersects_circle(e, f, circle2) );
}

void test_can_use_segment()
{
    test_header("test_can_use_segment");

    vector<Circle> circles;
    {
        Circle cA;
        cA.x = 0.0L; cA.y = 0.0L; cA.r = 1.0L;
        circles.push_back(cA);
    }
    {
        Circle cB;
        cB.x = 3.0L; cB.y = 3.0L; cB.r = 1.0L;
        circles.push_back(cB);
    }

    Point p1; p1.x = 2.0L; p1.y = 0.0L; p1.circle_index = -1;
    Point p2; p2.x = 4.0L; p2.y = 0.0L; p2.circle_index = -1;
    ASSERT_TRUE( can_use_segment(p1, p2, circles) );

    Point p3; p3.x = 2.0L; p3.y = 2.0L; p3.circle_index = -1;
    Point p4; p4.x = 4.0L; p4.y = 4.0L; p4.circle_index = -1;
    ASSERT_FALSE( can_use_segment(p3, p4, circles) );
}

void test_find_tangent_points()
{
    test_header("test_find_tangent_points");

    Circle c;
    c.x = 0.0L; c.y = 0.0L; c.r = 1.0L;

    Point p;
    p.x = 3.0L; p.y = 0.0L; p.circle_index = -1;

    vector<Point> tangents = find_tangent_points(p, 0, c);
    ASSERT_TRUE(tangents.size() == 2);

    for (size_t i = 0; i < tangents.size(); i++)
    {
        long double dx = tangents[i].x - c.x;
        long double dy = tangents[i].y - c.y;
        long double dist2 = dx * dx + dy * dy;
        ASSERT_NEAR(dist2, 1.0L, 1e-10);
    }
}

void test_can_use_arc()
{
    test_header("test_can_use_arc");

    Circle mainC;
    mainC.x = 0.0L; mainC.y = 0.0L; mainC.r = 2.0L;

    vector<Circle> circles;
    {
        Circle cX;
        cX.x = 3.0L; cX.y = 0.0L; cX.r = 1.0L;
        circles.push_back(cX);
    }

    Point center;
    center.x = mainC.x;
    center.y = mainC.y;
    center.circle_index = -1;

    long double R = mainC.r;

    bool ok_arc_ccw = can_use_arc(center, R, 0.0L, (long double)M_PI/2.0L, true, circles);
    ASSERT_TRUE(ok_arc_ccw);

    bool ok_arc_ccw2 = can_use_arc(center, R, 0.0L, 0.3L, true, circles);
    ASSERT_TRUE(ok_arc_ccw2);
}

void run_unit_tests()
{
    cout << "\nUNIT TESTY:\n";
    test_distance();
    test_inside_circle();
    test_segment_intersects_circle();
    test_can_use_segment();
    test_find_tangent_points();
    test_can_use_arc();
    cout << "Vsetky testovacie funkcie boli vykonane. :)\n\n";
}

/*
Examples of inputs and outputs:

1)
0 0 10 0
0

10.00000

2)
0 0 10 0
1
5 0 2

10.81122

3)
0 0 10 0
1
0 0 5

Path between A and B does not exist.

4)
0 0 10 0
2
3 0 2
7 0 2

11.39105

5)
0 0 10 0
2
3 0 3
7 0 3

13.42478

6)
0 0 10 10
1
5 5 3

15.43514

7)
0 0 10 0
2
2 2 2
8 2 2

10.00000

8)
0 0 10 0
8
3 3 2
-3 -3 2
-3 3 2
3 -3 2
0 3 2
3 0 2
0 -3 2
-3 0 2

Path between A and B does not exist.
 */


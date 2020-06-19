#include <iostream>
#include <queue>
#include <set>
#include <math.h>

#include "morph/tools.h"
#include "morph/HdfData.h"
#include "morph/Config.h"

using morph::Config;
using morph::Tools;
using morph::HdfData;

using namespace std;

// Notation for working with points
typedef pair<double, double> point;
#define x first
#define y second

// Arc, event, and segment datatypes
struct arc;
struct seg;

struct event {
    double x;
    point p;
    arc *a;
    bool valid;

    event(double xx, point pp, arc *aa)
    : x(xx), p(pp), a(aa), valid(true) {}
};

struct arc {
    point p;
    arc *prev, *next;
    event *e;

    seg *s0, *s1;

    arc(point pp, arc *a=0, arc *b=0)
    : p(pp), prev(a), next(b), e(0), s0(0), s1(0) {}
};

vector<seg*> output;  // Array of output segments.

struct seg {
    point start, end;
    bool done;

    seg(point p)
    : start(p), end(0,0), done(false)
    { output.push_back(this); }

    // Set the end point and mark as "done."
    void finish(point p) { if (done) return; end = p; done = true; }
};




class Voronoi {

public:

    int N;
    arc *root = 0; // First item in the parabolic front linked list.
    vector<double> x1, x2, y1, y2;
    vector<vector<int> > G;
    vector<vector<int> > S;


    Voronoi(int);

    // Function declarations
    void process_point();
    void process_event();
    void front_insert(point  p);

    bool circle(point a, point b, point c, double *x, point *o);
    void check_circle_event(arc *i, double x0);

    bool intersect(point p, arc *i, point *result = 0);
    point intersection(point p0, point p1, double l);

    void finish_edges();
    void print_output();
    void construct_graph();
    void save_points();

    // "Greater than" comparison, for reverse sorting in priority queue.
    struct gt {
        bool operator()(point a, point b) {return a.x==b.x ? a.y>b.y : a.x>b.x;}
        bool operator()(event *a, event *b) {return a->x>b->x;}
    };

    // Bounding box coordinates.
    double X0 = 0, X1 = 0, Y0 = 0, Y1 = 0;

    vector<double> px, py;

    priority_queue<point,  vector<point>,  gt> points; // site events
    priority_queue<event*, vector<event*>, gt> events; // circle events


};

Voronoi::Voronoi(int N){

    this->N = N;
    G.resize(N);
    S.resize(N);

    for(int i=0;i<N;i++){
        point p;
        p.x = morph::Tools::randDouble();
        p.y = morph::Tools::randDouble();
        points.push(p);

        px.push_back(p.x);
        py.push_back(p.y);

        // Keep track of bounding box size.
        if (p.x < X0) X0 = p.x;
        if (p.y < Y0) Y0 = p.y;
        if (p.x > X1) X1 = p.x;
        if (p.y > Y1) Y1 = p.y;
    }

    // Add margins to the bounding box.
    double dx = (X1-X0+1)/5.0, dy = (Y1-Y0+1)/5.0;
    X0 -= dx;  X1 += dx;  Y0 -= dy;  Y1 += dy;

    // Process the queues; select the top element with smaller x coordinate.
    while (!points.empty()){
        if (!events.empty() && events.top()->x <= points.top().x){ process_event(); }
        else { process_point();}
    }
    // After all points are processed, do the remaining circle events.
    while (!events.empty()){ process_event(); }

    finish_edges(); // Clean up dangling edges.

    construct_graph(); // populate G and S with graph information

}


void Voronoi::process_point()
{
    // Get the next point from the queue.
    point p = points.top();
    points.pop();

    // Add a new arc to the parabolic front.
    front_insert(p);
}

void Voronoi::process_event()
{
    // Get the next event from the queue.
    event *e = events.top();
    events.pop();

    if (e->valid) {
        // Start a new edge.
        seg *s = new seg(e->p);

        // Remove the associated arc from the front.
        arc *a = e->a;
        if (a->prev) {
            a->prev->next = a->next;
            a->prev->s1 = s;
        }
        if (a->next) {
            a->next->prev = a->prev;
            a->next->s0 = s;
        }

        // Finish the edges before and after a.
        if (a->s0) a->s0->finish(e->p);
        if (a->s1) a->s1->finish(e->p);

        // Recheck circle events on either side of p:
        if (a->prev) check_circle_event(a->prev, e->x);
        if (a->next) check_circle_event(a->next, e->x);
    }
    delete e;
}

void Voronoi::front_insert(point p)
{
    if (!root) {
        root = new arc(p);
        return;
    }

    // Find the current arc(s) at height p.y (if there are any).
    for (arc *i = root; i; i = i->next) {
        point z, zz;
        if (intersect(p,i,&z)) {
            // New parabola intersects arc i.  If necessary, duplicate i.
            if (i->next && !intersect(p,i->next, &zz)) {
                i->next->prev = new arc(i->p,i,i->next);
                i->next = i->next->prev;
            }
            else i->next = new arc(i->p,i);
            i->next->s1 = i->s1;

            // Add p between i and i->next.
            i->next->prev = new arc(p,i,i->next);
            i->next = i->next->prev;

            i = i->next; // Now i points to the new arc.

            // Add new half-edges connected to i's endpoints.
            i->prev->s1 = i->s0 = new seg(z);
            i->next->s0 = i->s1 = new seg(z);

            // Check for new circle events around the new arc:
            check_circle_event(i, p.x);
            check_circle_event(i->prev, p.x);
            check_circle_event(i->next, p.x);

            return;
        }
    }

    // Special case: If p never intersects an arc, append it to the list.
    arc *i;
    for (i = root; i->next; i=i->next) ; // Find the last node.

    i->next = new arc(p,i);
    // Insert segment between p and i
    point start;
    start.x = X0;
    start.y = (i->next->p.y + i->p.y) / 2;
    i->s1 = i->next->s0 = new seg(start);
}

// Look for a new circle event for arc i.
void Voronoi::check_circle_event(arc *i, double x0)
{
    // Invalidate any old event.
    if (i->e && i->e->x != x0)
    i->e->valid = false;
    i->e = NULL;

    if (!i->prev || !i->next)
    return;

    double x;
    point o;

    if (circle(i->prev->p, i->p, i->next->p, &x,&o) && x > x0) {
        // Create new event.
        i->e = new event(x, o, i);
        events.push(i->e);
    }
}

// Find the rightmost point on the circle through a,b,c.
bool Voronoi::circle(point a, point b, point c, double *x, point *o)
{
    // Check that bc is a "right turn" from ab.
    if ((b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y) > 0)
    return false;

    // Algorithm from O'Rourke 2ed p. 189.
    double A = b.x - a.x,  B = b.y - a.y,
    C = c.x - a.x,  D = c.y - a.y,
    E = A*(a.x+b.x) + B*(a.y+b.y),
    F = C*(a.x+c.x) + D*(a.y+c.y),
    G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

    if (G == 0) return false;  // Points are co-linear.

    // Point o is the center of the circle.
    o->x = (D*E-B*F)/G;
    o->y = (A*F-C*E)/G;

    // o.x plus radius equals max x coordinate.
    *x = o->x + sqrt( pow(a.x - o->x, 2) + pow(a.y - o->y, 2) );
    return true;
}

// Will a new parabola at point p intersect with arc i?
bool Voronoi::intersect(point p, arc *i, point *res)
{
    if (i->p.x == p.x) return false;

    double a,b;
    if (i->prev) // Get the intersection of i->prev, i.
    a = intersection(i->prev->p, i->p, p.x).y;
    if (i->next) // Get the intersection of i->next, i.
    b = intersection(i->p, i->next->p, p.x).y;

    if ((!i->prev || a <= p.y) && (!i->next || p.y <= b)) {
        res->y = p.y;

        // Plug it back into the parabola equation.
        res->x = (i->p.x*i->p.x + (i->p.y-res->y)*(i->p.y-res->y) - p.x*p.x)
        / (2*i->p.x - 2*p.x);

        return true;
    }
    return false;
}

// Where do two parabolas intersect?
point Voronoi::intersection(point p0, point p1, double l)
{
    point res, p = p0;

    if (p0.x == p1.x)
    res.y = (p0.y + p1.y) / 2;
    else if (p1.x == l)
    res.y = p1.y;
    else if (p0.x == l) {
        res.y = p0.y;
        p = p1;
    } else {
        // Use the quadratic formula.
        double z0 = 2*(p0.x - l);
        double z1 = 2*(p1.x - l);

        double a = 1/z0 - 1/z1;
        double b = -2*(p0.y/z0 - p1.y/z1);
        double c = (p0.y*p0.y + p0.x*p0.x - l*l)/z0
        - (p1.y*p1.y + p1.x*p1.x - l*l)/z1;

        res.y = ( -b - sqrt(b*b - 4*a*c) ) / (2*a);
    }
    // Plug back into one of the parabola equations.
    res.x = (p.x*p.x + (p.y-res.y)*(p.y-res.y) - l*l)/(2*p.x-2*l);
    return res;
}

void Voronoi::finish_edges()
{
    // Advance the sweep line so no parabolas can cross the bounding box.
    double l = X1 + (X1-X0) + (Y1-Y0);

    // Extend each remaining segment to the new parabola intersections.
    for (arc *i = root; i->next; i = i->next)
    if (i->s1)
    i->s1->finish(intersection(i->p, i->next->p, l*2));
}

void Voronoi::print_output()
{
    // Bounding box coordinates.
    cout << X0 << " "<< X1 << " " << Y0 << " " << Y1 << endl;

    // Each output segment in four-column format.
    vector<seg*>::iterator i;
    for (i = output.begin(); i != output.end(); i++) {
        point p0 = (*i)->start;
        point p1 = (*i)->end;
        cout << p0.x << " " << p0.y << " " << p1.x << " " << p1.y << endl;
    }
}

void Voronoi::construct_graph() {

    vector<double> xa, xb, ya, yb;

    vector<seg*>::iterator i;
    for (i = output.begin(); i != output.end(); i++) {
        point p0 = (*i)->start;
        point p1 = (*i)->end;
        xa.push_back(p0.x);
        ya.push_back(p0.y);
        xb.push_back(p1.x);
        yb.push_back(p1.y);
    }

    // Identify colinear segments (to be removed)
    vector<bool> flag(xa.size(),false);
    vector<vector<int> > colinear;
    for(int i=0;i<xa.size();i++){
        for(int j=0;j<xa.size();j++){
            if(i!=j){
                double ang1, ang2;
                if((xa[i]==xa[j]) && (ya[i]==ya[j])){
                    ang1 = atan2(yb[i]-ya[i],xb[i]-xa[i]);
                    ang2 = atan2(ya[j]-yb[j],xa[j]-xb[j]);
                }
                if((xa[i]==xb[j]) && (ya[i]==yb[j])){
                    ang1 = atan2(yb[i]-ya[i],xb[i]-xa[i]);
                    ang2 = atan2(yb[j]-ya[j],xb[j]-xa[j]);
                }
                if((xb[i]==xa[j]) && (yb[i]==ya[j])){
                    ang1 = atan2(ya[i]-yb[i],xa[i]-xb[i]);
                    ang2 = atan2(ya[j]-yb[j],xa[j]-xb[j]);
                }
                if((xb[i]==xb[j]) && (yb[i]==yb[j])){
                    ang1 = atan2(ya[i]-yb[i],xa[i]-xb[i]);
                    ang2 = atan2(yb[j]-ya[j],xb[j]-xa[j]);
                }
                if(fabs(ang1-ang2)<1e-9){
                    flag[i] = true;
                    vector<int> c(2,i); c[1] = j;
                    colinear.push_back(c);
                }
            }
        }
    }

    // append segments that will replace colinear segments
    for(int k=0;k<colinear.size();k++){
        int i = colinear[k][0];
        int j = colinear[k][1];
        if((xa[i]==xa[j]) && (ya[i]==ya[j])){
            xa.push_back(xb[i]);
            ya.push_back(yb[i]);
            xb.push_back(xb[j]);
            yb.push_back(yb[j]);
        }
        if((xa[i]==xb[j]) && (ya[i]==yb[j])){
            xa.push_back(xb[i]);
            ya.push_back(yb[i]);
            xb.push_back(xa[j]);
            yb.push_back(ya[j]);
        }
        if((xb[i]==xa[j]) && (yb[i]==ya[j])){
            xa.push_back(xa[i]);
            ya.push_back(ya[i]);
            xb.push_back(xb[j]);
            yb.push_back(yb[j]);
        }
        if((xb[i]==xb[j]) && (yb[i]==yb[j])){
            xa.push_back(xa[i]);
            ya.push_back(ya[i]);
            xb.push_back(xa[j]);
            yb.push_back(ya[j]);
        }
        flag.push_back(false);
    }

    // remove colinear segments
    for(int i=0;i<xa.size();i++){
        if(!flag[i]){
            x1.push_back(xa[i]);
            x2.push_back(xb[i]);
            y1.push_back(ya[i]);
            y2.push_back(yb[i]);
        }
    }



    // Obtain graph information in terms of connected points (equivalent to Delauny graph?) and the identity of the segments that partition each point from its neighbour
    vector<vector<int> > Gall(px.size());
    vector<vector<int> > Sall(px.size());

    int NN = x1.size();
    for(int i=0;i<NN;i++){

        // Identify nearest 3 seed points to end of segment
        vector<double> D(px.size());
        for(int j=0;j<px.size();j++){
            D[j] = (x2[i]-px[j])*(x2[i]-px[j]) + (y2[i]-py[j])*(y2[i]-py[j]);
        }
        vector<int> I(3);
        for(int k=0;k<3;k++){
            double minV = 1e9;
            int minI = 0;
            for(int j=0;j<px.size();j++){
                if(minV>D[j]){
                    minV = D[j];
                    minI = j;
                }
            }
            I[k] = minI;
            D[minI] = 1e9;
        }

        // rotate each seed 3 points about the angle of the segment
        double ang = atan2(y2[i]-y1[i],x2[i]-x1[i]);
        vector<double> r(3);
        for(int k=0;k<3;k++){
            r[k] = (px[I[k]]-x1[i])*sin(-ang) + (py[I[k]]-y1[i])*cos(-ang);
        }

        // each point in the pair connected by this segment are equidistant from y=0 after rotation
        for(int a=0;a<3;a++){
            for(int b=0;b<3;b++){
                if(fabs(r[a]+r[b])<1e-9){
                    Gall[I[a]].push_back(I[b]);
                    Gall[I[b]].push_back(I[a]);
                    Sall[I[a]].push_back(i);   // index into segments for each pair
                    Sall[I[b]].push_back(i);   // index into segments for each pair
                }
            }

        }
    }

    // clean up to get rid of duplicates
    for(int g=0;g<G.size();g++){
        for(int i=0;i<Gall[g].size();i++){
            bool uni = true;
            for(int k=0;k<G[g].size();k++){
                if(Gall[g][i]==G[g][k]){
                    uni = false; break;
                }
            } if(uni){
                G[g].push_back(Gall[g][i]);
                S[g].push_back(Sall[g][i]);
            }
        }
    }

}

void Voronoi::save_points() {

    HdfData outdata("points.h5");
    outdata.add_contained_vals ("x", px);
    outdata.add_contained_vals ("y", py);
    outdata.add_contained_vals ("x1", x1);
    outdata.add_contained_vals ("x2", x2);
    outdata.add_contained_vals ("y1", y1);
    outdata.add_contained_vals ("y2", y2);
}



struct cell {

public:

    int N;              // number of neighbours
    double x, y;        // coords of the center point
    vector<int> G, S;   // info relating to original Voronoi structure
    vector<double> X, Y, segL, segA, segX, segY; // coords of neighbours, and length and angle of adjoining segment

    double perim;
    vector<double> segLscaled;

    cell(Voronoi V, int index){
        x = V.px[index];
        y = V.py[index];
        G = V.G[index];
        S = V.S[index];
        N = G.size();
        segL.resize(N);
        segLscaled.resize(N);
        X.resize(N);
        Y.resize(N);
        segX.resize(N);
        segY.resize(N);
        segA.resize(N);

        for(int i=0;i<N;i++){
            X[i] = V.px[G[i]];
            Y[i] = V.py[G[i]];
            int s = S[i];
            double dx = V.x2[s]-V.x1[s];
            double dy = V.y2[s]-V.y1[s];
            segX[i] = V.x1[s];
            segY[i] = V.y1[s];
            segA[i] = atan2(dy,dx);
            segL[i] = pow(dx*dx+dy*dy,0.5);
        }


        for(int i=0;i<N;i++){
            if(segL[i]>pow(2.0,0.5)){
                segL[i] = 0.0;
            }
        }

        // PLAYING AROUND WITH SCALING TO USE FOR COUPLING - NEEDS THOUGHT!!
        perim = 0.;
        for(int i=0;i<N;i++){
            perim += segL[i];
        }
        for(int i=0;i<N;i++){
            segLscaled[i] = segL[i];//0.01*segL[i]/((double)N*perim);
        }

    }

    void save(){
        HdfData outdata("cell.h5");
        outdata.add_contained_vals ("X", X);
        outdata.add_contained_vals ("Y", Y);
        outdata.add_contained_vals ("segL", segL);
        outdata.add_contained_vals ("segX", segX);
        outdata.add_contained_vals ("segY", segY);
        outdata.add_contained_vals ("segA", segA);
    }
};

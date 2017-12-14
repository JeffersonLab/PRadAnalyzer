#ifndef CANA_GEOMETRY_H
#define CANA_GEOMETRY_H

#include <cmath>

namespace cana
{
    // get the intersection
    // of line (x1, y1) (x2, y2) and line (x3, y3) (x4, y4)
    // return status
    // -1 two lines are parallel, no intersection point found
    // 0 intersection point found, within both line segments
    // 1 intersection point found, out of the former line segment (x1, y1)&(x2, y2)
    // 2 intersection point found, out of the latter line segment (x3, y3)&(x4, y4)
    // 3 intersection point found, out of both line segments
    template<typename T>
    inline int intersection(T x1, T y1, T x2, T y2, T x3, T y3, T x4, T y4,
                            T &xc, T &yc, T inf = 0.001)
    {
        // denominator
        double denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);

        // parallel
        if(std::abs(denom) < inf) return -1;

        // find the cross point
        xc = ((x1*y2 - y1*x2)*(x3 - x4) - (x3*y4 - y3*x4)*(x1 - x2))/denom;
        yc = ((x1*y2 - y1*x2)*(y3 - y4) - (x3*y4 - y3*x4)*(y1 - y2))/denom;

        // check if the cross point is on the boundary line
        bool out_of_line1 = ((x1 - xc)*(xc - x2) < -inf) || ((y1 - yc)*(yc - y2) < -inf);
        bool out_of_line2 = ((x3 - xc)*(xc - x4) < -inf) || ((y3 - yc)*(yc - y4) < -inf);

        int status = (out_of_line1) ? 1 : 0;
        status += (out_of_line2) ? 2 : 0;
        return status;
    }

    // Check if a point is in a polygon
    // Implementation of the method described by Dan Sunday
    // Check http://geomalgorithms.com/a03-_inclusion.html for details
    // helper function to determine if the point is at the left side of a line
    // > 0 : (x,y) is at left of the line by (x0, y0) and (x1, y1)
    // = 0 : on the line
    // < 0 : right of the line
    template<typename T>
    inline T is_left(const T &x, const T &y,
                     const T &x0, const T &y0, const T &x1, const T &y1)
    {
        return (x1 - x0) * (y - y0) - (x - x0) * (y1 - y0);
    }

    // input: A 2d point and the iterators of polygon vertices
    //        Vertices order determine how the polygon is reconstructed
    //        The polygon is assumed to be closed, the last boundary connects
    //        its first and last vertices
    // output: winding number, 0 means outside
    template<class Iter, typename T>
    int inside_polygon_2d(T x, T y, Iter pv_beg, Iter pv_end)
    {
        int wn = 0;    // the  winding number counter

        // loop through all edges of the polygon
        Iter it_next = pv_beg;
        it_next++;
        for(Iter it = pv_beg; it != pv_end; it++, it_next++)
        {
            // last vertex, the boundary is formed between it and the first one
            if(it_next == pv_end)
                it_next = pv_beg;

            // start y <= point.y
            if(it->y <= y)
            {
                // upward crossing
                if(it_next->y > y)
                    // left of this boundary
                    if(is_left(x, y, it->x, it->y, it_next->x, it_next->y) > 0)
                        ++wn;
            }
            else
            {
                // downward crossing
                if(it_next->y <= y)
                    // right of this boundary
                    if(is_left(x, y, it->x, it->y, it_next->x, it_next->y) < 0)
                        --wn;
            }
        }

        // wn is positivie for counter clockwise and negative for clockwise
        return abs(wn);
    }
} // namespace cana

#endif // CANA_GEOMETRY_H

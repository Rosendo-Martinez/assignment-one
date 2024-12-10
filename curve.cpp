#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    
}
    

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> (" << P[i].x() << "," << P[i].y() << "," << P[i].z() << ")" << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;

    // NOTE: still not handling setting T, N, or B!

    // Bezier curves count
    const int curveCount = (P.size() - 4)/3  + 1;
    // Number of points returned curve should have
    const int curvePointsCount = steps * (curveCount) - (curveCount - 1);

    cout << "Curve count: " << curveCount << '\n';
    cout << "Sample Points Count: " << curvePointsCount << '\n';

    Curve cur(curvePointsCount);

    // Bernstein basis (B)
    const Matrix4f B(
        1, -3,  3, -1,
        0,  3, -6,  3,
        0,  0,  3, -3,
        0,  0,  0,  1
    );

    // Derivate of Bernstein basis (D)
    const Matrix4f D(
        -3,  3,  0, 0,
         0, -6,  6, 0,
         0,  0, -3, 3,
         0,  0,  0, 0
    );

    // sample t = 0 for first curve
    {
        // Control points matrix (G)
        const Matrix4f G(
            P[0].x(), P[1].x(), P[2].x(), P[3].x(),
            P[0].y(), P[1].y(), P[2].y(), P[3].y(),
            P[0].z(), P[1].z(), P[2].z(), P[3].z(),
            0,        0,        0,        0
        ); 

        // G * B
        const Matrix4f GB = G * B;
        // Monomial basis (T(t))
        const Vector4f T(1, 0, 0, 0);
        // G * B * T(t) 
        const Vector4f GBT = GB * T;
        CurvePoint cp;
        cp.V = Vector3f(GBT.x(), GBT.y(), GBT.z());
        cur[0] = cp;
    }

    // sample points from each curve (except t = 0)
    // two curves share 1 point, so don't sample it twice
    // index for current point to set in Curve
    int curvePointIndex = 1;
    for (int i = 0; i < curveCount; i++)
    {
        // index of first control point of current curve
        const int indexP0 = i * 3;

        // Control Point matrix (G)
        const Matrix4f G(
            P[indexP0].x(), P[indexP0 + 1].x(), P[indexP0 + 2].x(), P[indexP0 + 3].x(),
            P[indexP0].y(), P[indexP0 + 1].y(), P[indexP0 + 2].y(), P[indexP0 + 3].y(),
            P[indexP0].z(), P[indexP0 + 1].z(), P[indexP0 + 2].z(), P[indexP0 + 3].z(),
            0,              0,                  0,                  0
        );

        // G * B
        const Matrix4f GB = G * B;       

        // for each curve sample (steps - 1) points, exclude t = 0
        // each curve will sample t = 1, so the curve afterwards won't need to sample its t = 0
        for (int k = 0; k < steps - 1; k++)
        {
            const float t = ((float) (k + 1)) / (steps - 1);

            // Monomial basis (T(t))
            const Vector4f T(1, t, (t * t), (t * t * t));

            // G * B * T(t) 
            const Vector4f GBT = GB * T;

            CurvePoint cp;
            cp.V = Vector3f(GBT.x(), GBT.y(), GBT.z());
            cur[curvePointIndex] = cp;
            curvePointIndex++;
        }
    }

    return cur;
}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> (" << P[i].x() << "," << P[i].y() << "," << P[i].z() << ")" << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning empty curve." << endl;
    
    // For now assume 1 bezier curve (i.e 4 points)

    // B-spline basis (B_Spline)
    Matrix4f B_Spline(
        1/6.f, -3/6.f,  3/6.f, -1/6.f,
        4/6.f,  0/6.f, -6/6.f,  3/6.f,
        1/6.f,  3/6.f,  3/6.f, -3/6.f,
        0/6.f,  0/6.f,  0/6.f,  1/6.f
    );

    // Bernstein basis (B_Bezier)
    const Matrix4f B_Bezier(
        1, -3,  3, -1,
        0,  3, -6,  3,
        0,  0,  3, -3,
        0,  0,  0,  1
    );

    // B3 = B_spline * B_bezier^(-1)
    const Matrix4f B3 = B_Spline * B_Bezier.inverse();

    // B-Spline curve count
    const int numberOfBSplineCurves = P.size() - 3;
    // Count of control points for Bezier curve
    const int numberOfBezierControlPoints = (numberOfBSplineCurves - 1) * 3 + 4;
    vector<Vector3f> P_Bezier(numberOfBezierControlPoints);

    // First, set P_Bezier[0]
    {
        // Control points matrix (G_spline)
        const Matrix4f G_Spline(
            P[0].x(), P[1].x(), P[2].x(), P[3].x(),
            P[0].y(), P[1].y(), P[2].y(), P[3].y(),
            P[0].z(), P[1].z(), P[2].z(), P[3].z(),
            0,        0,        0,        0
        );
        // Change of basis from b-spline to bezier
        const Matrix4f G_Bezier = G_Spline * B3;
        
        P_Bezier[0] = G_Bezier.getCol(0).xyz();
    }

    // Then for each curve set only the last 3 of its respective control points
    // index of next control point in P_Bezier to set
    int pBezierIndex = 1;
    for (int i = 0; i < P.size() - 3; i++)
    {
        // Control points matrix (G_spline)
        const Matrix4f G_Spline(
            P[i].x(), P[i + 1].x(), P[i + 2].x(), P[i + 3].x(),
            P[i].y(), P[i + 1].y(), P[i + 2].y(), P[i + 3].y(),
            P[i].z(), P[i + 1].z(), P[i + 2].z(), P[i + 3].z(),
            0,        0,        0,        0
        );

        // Change of basis from b-spline to bezier
        const Matrix4f G_Bezier = G_Spline * B3;

        // only set last 3 control points of current curve
        for (int k = 0; k < 3; k++)
        {
            P_Bezier[pBezierIndex] = G_Bezier.getCol(k + 1).xyz();
            pBezierIndex++;
        }
    }

    return evalBezier(P_Bezier, steps);
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M.getElements() );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}


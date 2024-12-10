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

    // Print/log input information
    cerr << "\t>>> evalBezier has been called with the following input:" << endl;
    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> (" << P[i].x() << "," << P[i].y() << "," << P[i].z() << ")" << endl;
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


    // NOTE: still not handling setting T, N, or B!

    // Bezier curves count
    const unsigned curveCount = (P.size() - 4)/3  + 1;

    // Number of points returned curve should have
    const unsigned curvePointsCount = steps * (curveCount) - (curveCount - 1);

    // Curve / sample points to return
    Curve curve(curvePointsCount);

    // Bernstein basis (B)
    const Matrix4f B(
        1, -3,  3, -1,
        0,  3, -6,  3,
        0,  0,  3, -3,
        0,  0,  0,  1
    );

    // Derivate of Bernstein basis (D)
    const Matrix4f D(
        -3,   6, -3, 0,
         3, -12,  9, 0,
         0,   6, -9, 0,
         0,   0,  3, 0
    );

    // First, sample first point of curve
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

        // G * D
        const Matrix4f GD = G * D;

        // Monomial basis (T(t))
        const Vector4f T(1, 0, 0, 0);

        // Sample point
        // G * B * T(t) 
        const Vector4f GBT = GB * T;
        
        // Tangent
        // G * D * T(t)
        const Vector4f GDT = GD * T;

        CurvePoint cp;
        cp.V = Vector3f(GBT.x(), GBT.y(), GBT.z());
        cp.T = GDT.xyz();
        cp.T.normalize();

        // Arbitrary binormal
        Vector3f B0(1,0,0); 

        // Check that B0 and T aren't parallel
        if (approx(cp.T, B0))
        {
            B0 = Vector3f(0,1,0);
        }

        cp.N = Vector3f::cross(B0,cp.T).normalized();
        cp.B = Vector3f::cross(cp.T,cp.N).normalized();

        // set first sample point
        curve[0] = cp;
    }

    // Finally, sample (steps - 1) points from each curve

    // But don't sample t = 0 from any curve. Each curve
    // will be sampled at t = 1, which is same sample point for next 
    // curve at t = 0, thats why we skip t = 0.

    // Index of next sample point to set
    unsigned u = 1;

    for (unsigned i = 0; i < curveCount; i++)
    {
        // index of first control point of current curve
        const unsigned p_0 = i * 3;

        // Control Point matrix (G)
        const Matrix4f G(
            P[p_0].x(), P[p_0 + 1].x(), P[p_0 + 2].x(), P[p_0 + 3].x(),
            P[p_0].y(), P[p_0 + 1].y(), P[p_0 + 2].y(), P[p_0 + 3].y(),
            P[p_0].z(), P[p_0 + 1].z(), P[p_0 + 2].z(), P[p_0 + 3].z(),
            0,              0,                  0,                  0
        );

        // G * B
        const Matrix4f GB = G * B;

        // G * D
        const Matrix4f GD = G * D;

        // Sample points from the curve
        for (unsigned k = 0; k < steps - 1; k++)
        {
            // Current value of parameter t
            const float t = ((float) (k + 1)) / (steps - 1);

            // Monomial basis (T(t))
            const Vector4f T(1, t, (t * t), (t * t * t));

            // Sample point
            // Q(t) = G * B * T(t)
            const Vector4f GBT = GB * T;

            // Tangent
            // G * D * T(t)
            const Vector4f GDT = GD * T;

            CurvePoint cp;
            cp.V = GBT.xyz();
            cp.T = GDT.xyz();
            cp.T.normalize();

            cp.N = Vector3f::cross(curve[u-1].B,cp.T).normalized();
            cp.B = Vector3f::cross(cp.T,cp.N).normalized();

            // add sample point to Curve
            curve[u] = cp;

            // increment to next sample point index
            u++;
        }
    }

    return curve;
}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // Print input information
    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;
    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> (" << P[i].x() << "," << P[i].y() << "," << P[i].z() << ")" << endl;
    }

    // This function simply does of change of basis from B-Spline to Bezier, then calls
    // evalBezier to do rest of work.
    
    // B-spline basis
    const Matrix4f B1(
        1/6.f, -3/6.f,  3/6.f, -1/6.f,
        4/6.f,  0/6.f, -6/6.f,  3/6.f,
        1/6.f,  3/6.f,  3/6.f, -3/6.f,
        0/6.f,  0/6.f,  0/6.f,  1/6.f
    );

    // Bernstein/Bezier basis
    const Matrix4f B2(
        1, -3,  3, -1,
        0,  3, -6,  3,
        0,  0,  3, -3,
        0,  0,  0,  1
    );

    // B1 * B2^(-1)
    const Matrix4f B3 = B1 * B2.inverse();

    // B-Spline curve count
    const unsigned bSplineCurvesCount = P.size() - 3;

    // Count of control points for Bezier curves
    const unsigned bezierControlPointsCount = (bSplineCurvesCount - 1) * 3 + 4;

    // Control points in Bezier basis
    vector<Vector3f> P2(bezierControlPointsCount);

    // First, do change of basis to Bezier for first control point only
    {
        // Control points matrix in b-spline basis
        const Matrix4f G1(
            P[0].x(), P[1].x(), P[2].x(), P[3].x(),
            P[0].y(), P[1].y(), P[2].y(), P[3].y(),
            P[0].z(), P[1].z(), P[2].z(), P[3].z(),
            0,        0,        0,        0
        );

        // Control points matrix in Bezier basis
        const Matrix4f G2 = G1 * B3;
        
        // set first control point
        P2[0] = G2.getCol(0).xyz();
    }

    // Then, for each curve set only the last 3 of its respective control points

    // index of next Bezier control point to set
    unsigned u = 1;

    for (unsigned i = 0; i < P.size() - 3; i++)
    {
        // Control points matrix in B-Spline basis
        const Matrix4f G1(
            P[i].x(), P[i + 1].x(), P[i + 2].x(), P[i + 3].x(),
            P[i].y(), P[i + 1].y(), P[i + 2].y(), P[i + 3].y(),
            P[i].z(), P[i + 1].z(), P[i + 2].z(), P[i + 3].z(),
            0,        0,        0,        0
        );

        // Control points matrix in Bezier basis
        const Matrix4f G2 = G1 * B3;

        // only set last 3 control points of current curve
        for (unsigned k = 0; k < 3; k++)
        {
            P2[u] = G2.getCol(k + 1).xyz();
            u++;
        }
    }

    return evalBezier(P2, steps);
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


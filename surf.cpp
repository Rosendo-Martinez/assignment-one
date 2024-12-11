#include "surf.h"
#include "extra.h"
#include <cmath>
using namespace std;

namespace
{
    
    // We're only implenting swept surfaces where the profile curve is
    // flat on the xy-plane.  This is a check function.
    static bool checkFlat(const Curve &profile)
    {
        for (unsigned i=0; i<profile.size(); i++)
            if (profile[i].V[2] != 0.0 ||
                profile[i].T[2] != 0.0 ||
                profile[i].N[2] != 0.0)
                return false;
    
        return true;
    }
}

Surface makeSurfRev(const Curve &profile, unsigned steps)
{
    // Some assumptions are made, if they are  not kept,
    // then likely will produce incorrect results.

    Surface surface;
    
    if (!checkFlat(profile))
    {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // Two main steps:
    // 1) create rotated curves
    // 2) triangulate the curves

    // Create curves of surface
    for (unsigned i = 0; i < steps; i++)
    {
        // rotation angle
        const double a = (2*M_PI)/(steps) * i;

        // Rotation matrix (R):
        // counter-clock wise rotation around +y axis
        const Matrix3f R(
            cos(a),  0, sin(a),
            0,       1, 0,
            -sin(a), 0, cos(a)
        );

        // rotate profile curve
        for (unsigned j = 0; j < profile.size(); j++)
        {
            // Rotate vertex (V)
            const Vector3f V = R * profile[j].V;
            // Rotate normal (N)
            const Vector3f N = R * (-1 * profile[j].N);

            // add to surface
            surface.VV.push_back(V);
            surface.VN.push_back(N);
        }
    }

    // generate 'upper' triangles
    for (unsigned i = 0; i < surface.VV.size(); i++)
    {
        // Vertex is last point on its curve
        // so it can't form a upper triangle with it as top left
        if ((i + 1) % profile.size() == 0)
        {
            continue;
        }

        // Upper Triangle
        // A-----C
        // |   /
        // | /
        // B

        // Note, A and B are on same curve, and C is on adjacent curve.

        // A is i
        // B is i + 1
        // C is i + profile.size()

        surface.VF.push_back(Tup3u(i, i+1, (i+profile.size()) % surface.VV.size()));
    }

    // generate 'lower' triangles
    for (unsigned i = 0; i < surface.VV.size(); i++)
    {
        // Vertex is first point on its curve
        // so it can't form a lower triangle with it as bottom left
        if (i % profile.size() == 0)
        {
            continue;
        }

        /**
         * Lower Triangle
         *       C
         *      /|
         *    /  |
         *  /    |
         * A-----B
         * 
         * Note, B and C are on same curve, and A is on adjacent curve.
         */

        // A is i
        // B is i + profile.size()
        // C is i + profile.size() - 1

        surface.VF.push_back(Tup3u(i, (i + profile.size()) % surface.VV.size(), (i+profile.size() - 1) % surface.VV.size()));   
    }
 
    return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep )
{
    Surface surface;

    if (!checkFlat(profile))
    {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    // TODO: Here you should build the surface.  See surf.h for details.

    cerr << "\t>>> makeGenCyl called (but not implemented).\n\t>>> Returning empty surface." <<endl;

    return surface;
}

void drawSurface(const Surface &surface, bool shaded)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    if (shaded)
    {
        // This will use the current material color and light
        // positions.  Just set these in drawScene();
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // This tells openGL to *not* draw backwards-facing triangles.
        // This is more efficient, and in addition it will help you
        // make sure that your triangles are drawn in the right order.
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    else
    {        
        glDisable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        glColor4f(0.4f,0.4f,0.4f,1.f);
        glLineWidth(1);
    }

    glBegin(GL_TRIANGLES);
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        glNormal(surface.VN[surface.VF[i][0]]);
        glVertex(surface.VV[surface.VF[i][0]]);
        glNormal(surface.VN[surface.VF[i][1]]);
        glVertex(surface.VV[surface.VF[i][1]]);
        glNormal(surface.VN[surface.VF[i][2]]);
        glVertex(surface.VV[surface.VF[i][2]]);
    }
    glEnd();

    glPopAttrib();
}

void drawNormals(const Surface &surface, float len)
{
    // Save current state of OpenGL
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_LIGHTING);
    glColor4f(0,1,1,1);
    glLineWidth(1);

    glBegin(GL_LINES);
    for (unsigned i=0; i<surface.VV.size(); i++)
    {
        glVertex(surface.VV[i]);
        glVertex(surface.VV[i] + surface.VN[i] * len);
    }
    glEnd();

    glPopAttrib();
}

void outputObjFile(ostream &out, const Surface &surface)
{
    
    for (unsigned i=0; i<surface.VV.size(); i++)
        out << "v  "
            << surface.VV[i][0] << " "
            << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (unsigned i=0; i<surface.VN.size(); i++)
        out << "vn "
            << surface.VN[i][0] << " "
            << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;
    
    for (unsigned i=0; i<surface.VF.size(); i++)
    {
        out << "f  ";
        for (unsigned j=0; j<3; j++)
        {
            unsigned a = surface.VF[i][j]+1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "searchableCone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableCone, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableCone,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableCone,
        dict,
        cone
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::searchableCone::coordinates() const
{
    return tmp<pointField>::New(1, 0.5*(point1_ + point2_));
}


void Foam::searchableCone::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = 0.5*(point1_ + point2_);

    radiusSqr.setSize(1);
    if (radius1_ > radius2_)
    {
        radiusSqr[0] = Foam::magSqr(point1_-centres[0]) + Foam::sqr(radius1_);
    }
    else
    {
        radiusSqr[0] = Foam::magSqr(point2_-centres[0]) + Foam::sqr(radius2_);
    }

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


Foam::tmp<Foam::pointField> Foam::searchableCone::points() const
{
    auto tpts = tmp<pointField>::New(2);
    auto& pts = tpts.ref();

    pts[0] = point1_;
    pts[1] = point2_;

    return tpts;
}


void Foam::searchableCone::findNearestAndNormal
(
    const point& sample,
    const scalar nearestDistSqr,
    pointIndexHit& info,
    vector& nearNormal
) const
{
    vector v(sample - point1_);

    // Decompose sample-point1 into normal and parallel component
    const scalar parallel = (v & unitDir_);

    // Remove the parallel component and normalise
    v -= parallel*unitDir_;

    const scalar magV = mag(v);
    v.normalise();

    // Nearest and normal on disk at point1
    point disk1Point(point1_ + min(max(magV, innerRadius1_), radius1_)*v);
    vector disk1Normal(-unitDir_);

    // Nearest and normal on disk at point2
    point disk2Point(point2_ + min(max(magV, innerRadius2_), radius2_)*v);
    vector disk2Normal(unitDir_);

    // Nearest and normal on cone. Initialise to far-away point so if not
    // set it picks one of the disk points
    point nearCone(point::uniform(GREAT));
    vector normalCone(1, 0, 0);

    // Nearest and normal on inner cone. Initialise to far-away point
    point iCnearCone(point::uniform(GREAT));
    vector iCnormalCone(1, 0, 0);

    point projPt1 = point1_+ radius1_*v;
    point projPt2 = point2_+ radius2_*v;

    vector p1 = (projPt1 - point1_);
    if (mag(p1) > ROOTVSMALL)
    {
        p1 /= mag(p1);

        // Find vector along the two end of cone
        const vector b = normalised(projPt2 - projPt1);

        // Find the vector along sample pt and pt at one end of cone
        vector a = (sample - projPt1);

        if (mag(a) <= ROOTVSMALL)
        {
            // Exception: sample on disk1. Redo with projPt2.
            a = (sample - projPt2);

            // Find normal unitvector
            nearCone = (a & b)*b + projPt2;

            vector b1 = (p1 & b)*b;
            normalCone = normalised(p1 - b1);
        }
        else
        {
            // Find nearest point on cone surface
            nearCone = (a & b)*b + projPt1;

            // Find projection along surface of cone
            vector b1 = (p1 & b)*b;
            normalCone = normalised(p1 - b1);
        }

        if (innerRadius1_ > 0 || innerRadius2_ > 0)
        {
            // Same for inner radius
            point iCprojPt1 = point1_+ innerRadius1_*v;
            point iCprojPt2 = point2_+ innerRadius2_*v;

            const vector iCp1 = normalised(iCprojPt1 - point1_);

            // Find vector along the two end of cone
            const vector iCb = normalised(iCprojPt2 - iCprojPt1);


            // Find the vector along sample pt and pt at one end of conde
            vector iCa(sample - iCprojPt1);

            if (mag(iCa) <= ROOTVSMALL)
            {
                iCa = (sample - iCprojPt2);

                // Find normal unitvector
                iCnearCone = (iCa & iCb)*iCb+iCprojPt2;

                vector b1 = (iCp1 & iCb)*iCb;
                iCnormalCone = normalised(iCp1 - b1);
            }
            else
            {
                // Find nearest point on cone surface
                iCnearCone = (iCa & iCb)*iCb+iCprojPt1;

                // Find projection along surface of cone
                vector b1 = (iCp1 & iCb)*iCb;
                iCnormalCone = normalised(iCp1 - b1);
            }
        }
    }


    // Select nearest out of the 4 points (outer cone, disk1, disk2, inner
    // cone)

    FixedList<scalar, 4> dist;
    dist[0] = magSqr(nearCone-sample);
    dist[1] = magSqr(disk1Point-sample);
    dist[2] = magSqr(disk2Point-sample);
    dist[3] = magSqr(iCnearCone-sample);

    const label minI = findMin(dist);


    // Snap the point to the corresponding surface

    if (minI == 0)  // Outer cone
    {
        // Closest to (infinite) outer cone. See if needs clipping to end disks

        {
            vector v1(nearCone-point1_);
            scalar para = (v1 & unitDir_);
            // Remove the parallel component and normalise
            v1 -= para*unitDir_;
            const scalar magV1 = mag(v1);
            v1 = v1/magV1;

            if (para < 0.0 && magV1 >= radius1_)
            {
                // Near point 1. Set point to intersection of disk and cone.
                // Keep normal from cone.
                nearCone = disk1Point;
            }
            else if (para < 0.0 && magV1 < radius1_)
            {
                // On disk1
                nearCone = disk1Point;
                normalCone = disk1Normal;
            }
            else if (para > magDir_ && magV1 >= radius2_)
            {
                // Near point 2. Set point to intersection of disk and cone.
                // Keep normal from cone.
                nearCone = disk2Point;
            }
            else if (para > magDir_ && magV1 < radius2_)
            {
                // On disk2
                nearCone = disk2Point;
                normalCone = disk2Normal;
            }
        }
        info.setPoint(nearCone);
        nearNormal = normalCone;
    }
    else if (minI == 1) // Near to disk1
    {
        info.setPoint(disk1Point);
        nearNormal = disk1Normal;
    }
    else if (minI == 2) // Near to disk2
    {
        info.setPoint(disk2Point);
        nearNormal = disk2Normal;
    }
    else                // Near to inner cone
    {
        {
            vector v1(iCnearCone-point1_);
            scalar para = (v1 & unitDir_);
            // Remove the parallel component and normalise
            v1 -= para*unitDir_;

            const scalar magV1 = mag(v1);
            v1 = v1/magV1;

            if (para < 0.0 && magV1 >= innerRadius1_)
            {
               iCnearCone = disk1Point;
            }
            else if (para < 0.0 && magV1 < innerRadius1_)
            {
               iCnearCone = disk1Point;
               iCnormalCone = disk1Normal;
            }
            else if (para > magDir_ && magV1 >= innerRadius2_)
            {
               iCnearCone = disk2Point;
            }
            else if (para > magDir_ && magV1 < innerRadius2_)
            {
               iCnearCone = disk2Point;
               iCnormalCone = disk2Normal;
            }
        }
        info.setPoint(iCnearCone);
        nearNormal = iCnormalCone;
    }


    if (magSqr(sample - info.rawPoint()) < nearestDistSqr)
    {
        info.setHit();
        info.setIndex(0);
    }
}


Foam::scalar Foam::searchableCone::radius2
(
    const searchableCone& cone,
    const point& pt
)
{
    const vector x = (pt-cone.point1_) ^ cone.unitDir_;
    return x&x;
}


// From mrl.nyu.edu/~dzorin/rend05/lecture2.pdf,
// Ray Tracing II, Infinite cone ray intersection.
void Foam::searchableCone::findLineAll
(
    const searchableCone& cone,
    const scalar innerRadius1,
    const scalar innerRadius2,
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    vector point1Start(start-cone.point1_);
    vector point2Start(start-cone.point2_);
    vector point1End(end-cone.point1_);

    // Quick rejection of complete vector outside endcaps
    scalar s1 = point1Start & (cone.unitDir_);
    scalar s2 = point1End & (cone.unitDir_);

    if ((s1 < 0.0 && s2 < 0.0) || (s1 > cone.magDir_ && s2 > cone.magDir_))
    {
        return;
    }

    // Line as P = start+t*V where V is unit vector and t=[0..mag(end-start)]
    vector V(end-start);
    scalar magV = mag(V);
    if (magV < ROOTVSMALL)
    {
        return;
    }
    V /= magV;


    // We now get the nearest intersections to start. This can either be
    // the intersection with the end plane or with the cylinder side.

    // Get the two points (expressed in t) on the end planes. This is to
    // clip any cylinder intersection against.
    scalar tPoint1;
    scalar tPoint2;

    // Maintain the two intersections with the endcaps
    scalar tNear = VGREAT;
    scalar tFar = VGREAT;

    scalar radius_sec = cone.radius1_;

    {
        // Find dot product: mag(s)>VSMALL suggest that it is greater
        scalar s = (V & unitDir_);
        if (mag(s) > VSMALL)
        {
            tPoint1 = -s1/s;
            tPoint2 = -(point2Start&(cone.unitDir_))/s;

            if (tPoint2 < tPoint1)
            {
                std::swap(tPoint1, tPoint2);
            }
            if (tPoint1 > magV || tPoint2 < 0)
            {
                return;
            }
            // See if the points on the endcaps are actually inside the cylinder
            if (tPoint1 >= 0 && tPoint1 <= magV)
            {
                scalar r2 = radius2(cone, start+tPoint1*V);
                vector p1 = (start+tPoint1*V-point1_);
                vector p2 = (start+tPoint1*V-point2_);
                radius_sec = cone.radius1_;
                scalar inC_radius_sec = innerRadius1_;

                if (mag(p2&(cone.unitDir_)) < mag(p1&(cone.unitDir_)))
                {
                    radius_sec = cone.radius2_;
                    inC_radius_sec = innerRadius2_;
                }

                if (r2 <= sqr(radius_sec) && r2 >= sqr(inC_radius_sec))
                {
                    tNear = tPoint1;
                }
            }
            if (tPoint2 >= 0 && tPoint2 <= magV)
            {
                // Decompose sample-point1 into normal and parallel component
                vector p1 = (start+tPoint2*V-cone.point1_);
                vector p2 = (start+tPoint2*V-cone.point2_);
                radius_sec = cone.radius1_;
                scalar inC_radius_sec = innerRadius1_;
                if (mag(p2&(cone.unitDir_)) < mag(p1&(cone.unitDir_)))
                {
                    radius_sec = cone.radius2_;
                    inC_radius_sec = innerRadius2_;
                }
                scalar r2 = radius2(cone, start+tPoint2*V);
                if (r2 <= sqr(radius_sec) && r2 >= sqr(inC_radius_sec))
                {
                    // Check if already have a near hit from point1
                    if (tNear <= 1)
                    {
                        tFar = tPoint2;
                    }
                    else
                    {
                        tNear = tPoint2;
                    }
                }
            }
        }
        else
        {
            // Vector perpendicular to cylinder. Check for outside already done
            // above so just set tpoint to allow all.
            tPoint1 = -VGREAT;
            tPoint2 = VGREAT;
        }
    }


    // Second order equation of the form a*t^2 + b*t + c
    scalar a, b, c;

    scalar deltaRadius = cone.radius2_-cone.radius1_;
    if (mag(deltaRadius) <= ROOTVSMALL)
    {
        vector point1Start(start-cone.point1_);
        const vector x = point1Start ^ cone.unitDir_;
        const vector y = V ^ cone.unitDir_;
        const scalar d = sqr(0.5*(cone.radius1_ + cone.radius2_));

        a = (y&y);
        b = 2*(x&y);
        c = (x&x)-d;
    }
    else
    {
        vector va = cone.unitDir_;
        vector v1 = normalised(end-start);
        scalar p  = (va&v1);
        vector a1 = (v1-p*va);

        // Determine the end point of the cone
        point pa =
            cone.unitDir_*cone.radius1_*mag(cone.point2_-cone.point1_)
           /(-deltaRadius)
           +cone.point1_;

        scalar l2 = sqr(deltaRadius)+sqr(cone.magDir_);
        scalar sqrCosAlpha = sqr(cone.magDir_)/l2;
        scalar sqrSinAlpha = sqr(deltaRadius)/l2;


        vector delP(start-pa);
        vector p1 = (delP-(delP&va)*va);

        a = sqrCosAlpha*((v1-p*va)&(v1-p*va))-sqrSinAlpha*p*p;
        b =
            2.0*sqrCosAlpha*(a1&p1)
           -2.0*sqrSinAlpha*(v1&va)*(delP&va);
        c =
            sqrCosAlpha
           *((delP-(delP&va)*va)&(delP-(delP&va)*va))
           -sqrSinAlpha*sqr(delP&va);
    }

    const scalar disc = b*b-4.0*a*c;

    scalar t1 = -VGREAT;
    scalar t2 = VGREAT;

    if (disc < 0)
    {
        // Fully outside
        return;
    }
    else if (disc < ROOTVSMALL)
    {
        // Single solution
        if (mag(a) > ROOTVSMALL)
        {
            t1 = -b/(2.0*a);

            if (t1 >= 0.0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // Valid. Insert sorted.
                if (t1 < tNear)
                {
                    tFar = tNear;
                    tNear = t1;
                }
                else if (t1 < tFar)
                {
                    tFar = t1;
                }
            }
            else
            {
                return;
            }
        }
        else
        {
            // Aligned with axis. Check if outside radius
            if (c > 0.0)
            {
                return;
            }
        }
    }
    else
    {
        if (mag(a) > ROOTVSMALL)
        {
            scalar sqrtDisc = sqrt(disc);

            t1 = (-b - sqrtDisc)/(2.0*a);
            t2 = (-b + sqrtDisc)/(2.0*a);
            if (t2 < t1)
            {
                std::swap(t1, t2);
            }

            if (t1 >= 0.0 && t1 <= magV && t1 >= tPoint1 && t1 <= tPoint2)
            {
                // Valid. Insert sorted.
                if (t1 < tNear)
                {
                    tFar = tNear;
                    tNear = t1;
                }
                else if (t1 < tFar)
                {
                    tFar = t1;
                }
            }
            if (t2>=0 && t2 <= magV && t2 >= tPoint1 && t2 <= tPoint2)
            {
                // Valid. Insert sorted.
                if (t2 < tNear)
                {
                    tFar = tNear;
                    tNear = t2;
                }
                else if (t2 < tFar)
                {
                    tFar = t2;
                }
            }
        }
        else
        {
            // Aligned with axis. Check if outside radius
            if (c > 0.0)
            {
                return;
            }
        }
    }

    // Check tNear, tFar
    if (tNear>=0 && tNear <= magV)
    {
        near.setPoint(start+tNear*V);
        near.setHit();
        near.setIndex(0);
        if (tFar <= magV)
        {
            far.setPoint(start+tFar*V);
            far.setHit();
            far.setIndex(0);
        }
    }
    else if (tFar>=0 && tFar <= magV)
    {
        near.setPoint(start+tFar*V);
        near.setHit();
        near.setIndex(0);
    }
}


void Foam::searchableCone::insertHit
(
    const point& start,
    const point& end,
    List<pointIndexHit>& info,
    const pointIndexHit& hit
) const
{
    scalar smallDistSqr = SMALL*magSqr(end-start);

    scalar hitMagSqr = magSqr(hit.hitPoint()-start);

    forAll(info, i)
    {
        scalar d2 = magSqr(info[i].hitPoint()-start);

        if (d2 > hitMagSqr+smallDistSqr)
        {
            // Insert at i.
            label sz = info.size();
            info.setSize(sz+1);
            for (label j = sz; j > i; --j)
            {
                info[j] = info[j-1];
            }
            info[i] = hit;
            return;
        }
        else if (d2 > hitMagSqr-smallDistSqr)
        {
            // hit is same point as info[i].
            return;
        }
    }
    // Append
    label sz = info.size();
    info.setSize(sz+1);
    info[sz] = hit;
}


Foam::boundBox Foam::searchableCone::calcBounds() const
{
    // Adapted from
    // http://www.gamedev.net/community/forums
    //       /topic.asp?topic_id=338522&forum_id=20&gforum_id=0

    // Let cylinder have end points A,B and radius r,

    // Bounds in direction X (same for Y and Z) can be found as:
    // Let A.X<B.X (otherwise swap points)
    // Good approximate lowest bound is A.X-r and highest is B.X+r (precise for
    // capsule). At worst, in one direction it can be larger than needed by 2*r.

    // Accurate bounds for cylinder is
    // A.X-kx*r, B.X+kx*r
    // where
    // kx=sqrt(((A.Y-B.Y)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))

    // similar thing for Y and Z
    // (i.e.
    // ky=sqrt(((A.X-B.X)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // kz=sqrt(((A.X-B.X)^2+(A.Y-B.Y)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // )

    // How derived: geometric reasoning. Bounds of cylinder is same as for 2
    // circles centered on A and B. This sqrt thingy gives sine of angle between
    // axis and direction, used to find projection of radius.

    vector kr
    (
        sqrt(sqr(unitDir_.y()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.y()))
    );

    if (radius2_ >= radius1_)
    {
        kr *= radius2_;
    }
    else
    {
        kr *= radius1_;
    }

    point min = point1_ - kr;
    point max = point1_ + kr;

    min = ::Foam::min(min, point2_ - kr);
    max = ::Foam::max(max, point2_ + kr);

    return boundBox(min, max);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableCone::searchableCone
(
    const IOobject& io,
    const point& point1,
    const scalar radius1,
    const scalar innerRadius1,
    const point& point2,
    const scalar radius2,
    const scalar innerRadius2
)
:
    searchableSurface(io),
    point1_(point1),
    radius1_(radius1),
    innerRadius1_(innerRadius1),
    point2_(point2),
    radius2_(radius2),
    innerRadius2_(innerRadius2),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_)
{
    bounds() = calcBounds();
}


Foam::searchableCone::searchableCone
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    point1_(dict.get<point>("point1")),
    radius1_(dict.get<scalar>("radius1")),
    innerRadius1_(dict.getOrDefault<scalar>("innerRadius1", 0)),
    point2_(dict.get<point>("point2")),
    radius2_(dict.get<scalar>("radius2")),
    innerRadius2_(dict.getOrDefault<scalar>("innerRadius2", 0)),
    magDir_(mag(point2_-point1_)),
    unitDir_((point2_-point1_)/magDir_)
{
    bounds() = calcBounds();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableCone::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchableCone::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        vector normal;
        findNearestAndNormal(samples[i], nearestDistSqr[i], info[i], normal);
    }
}


void Foam::searchableCone::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        pointIndexHit b;
        findLineAll
        (
            *this,
            innerRadius1_,
            innerRadius2_,
            start[i],
            end[i],
            info[i],
            b
        );
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }


    if (innerRadius1_ > 0.0 || innerRadius2_ > 0.0)
    {
        IOobject io(*this);
        io.rename(name()+"Inner");
        searchableCone innerCone
        (
            io,
            point1_,
            innerRadius1_,
            0.0,
            point2_,
            innerRadius2_,
            0.0
        );

        forAll(info, i)
        {
            point newEnd;
            if (info[i].hit())
            {
                newEnd = info[i].hitPoint();
            }
            else
            {
                newEnd = end[i];
            }
            pointIndexHit near;
            pointIndexHit far;
            findLineAll
            (
                innerCone,
                innerRadius1_,
                innerRadius2_,
                start[i],
                newEnd,
                near,
                far
            );

            if (near.hit())
            {
                info[i] = near;
            }
            else if (far.hit())
            {
                info[i] = far;
            }
        }
    }
}


void Foam::searchableCone::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Discard far intersection
        pointIndexHit b;
        findLineAll
        (
            *this,
            innerRadius1_,
            innerRadius2_,
            start[i],
            end[i],
            info[i],
            b
        );
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
    if (innerRadius1_ > 0.0 || innerRadius2_ > 0.0)
    {
        IOobject io(*this);
        io.rename(name()+"Inner");
        searchableCone cone
        (
            io,
            point1_,
            innerRadius1_,
            0.0,
            point2_,
            innerRadius2_,
            0.0
        );

        forAll(info, i)
        {
            if (!info[i].hit())
            {
                pointIndexHit b;
                findLineAll
                (
                    cone,
                    innerRadius1_,
                    innerRadius2_,
                    start[i],
                    end[i],
                    info[i],
                    b
                );
                if (!info[i].hit() && b.hit())
                {
                    info[i] = b;
                }
            }
        }
    }
}


void Foam::searchableCone::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll
        (
            *this,
            innerRadius1_,
            innerRadius2_,
            start[i],
            end[i],
            near,
            far
        );

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }

    if (innerRadius1_ > 0.0 || innerRadius2_ > 0.0)
    {
        IOobject io(*this);
        io.rename(name()+"Inner");
        searchableCone cone
        (
            io,
            point1_,
            innerRadius1_,
            0.0,
            point2_,
            innerRadius2_,
            0.0
        );

        forAll(info, i)
        {
            pointIndexHit near;
            pointIndexHit far;
            findLineAll
            (
                cone,
                innerRadius1_,
                innerRadius2_,
                start[i],
                end[i],
                near,
                far
            );

            if (near.hit())
            {
                insertHit(start[i], end[i], info[i], near);
            }
            if (far.hit())
            {
                insertHit(start[i], end[i], info[i], far);
            }
        }
    }
}


void Foam::searchableCone::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableCone::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            pointIndexHit nearInfo;
            findNearestAndNormal
            (
                info[i].hitPoint(),
                Foam::sqr(GREAT),
                nearInfo,
                normal[i]
            );
        }
    }
}


void Foam::searchableCone::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        volType[pointi] = volumeType::OUTSIDE;

        vector v(pt - point1_);

        // Decompose sample-point1 into normal and parallel component
        const scalar parallel = (v & unitDir_);

        // Quick rejection. Left of point1 endcap, or right of point2 endcap
        if (parallel < 0 || parallel > magDir_)
        {
            continue;
        }

        const scalar radius_sec =
            radius1_ + parallel * (radius2_-radius1_)/magDir_;

        const scalar radius_sec_inner =
            innerRadius1_ + parallel * (innerRadius2_-innerRadius1_)/magDir_;

        // Remove the parallel component
        v -= parallel*unitDir_;

        if (mag(v) >= radius_sec_inner && mag(v) <= radius_sec)
        {
            volType[pointi] = volumeType::INSIDE;
        }
    }
}


// ************************************************************************* //

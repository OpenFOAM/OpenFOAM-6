/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

Class
    Foam::searchableCone

Description
    Searching on (optionally hollow) cone.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required | Default
        type        | cone                              | selector |
        point1      | coordinate of endpoint            | yes   |
        radius1     | radius at point1                  | yes   |
        innerRadius1| inner radius at point1            | no    | 0
        point2      | coordinate of endpoint            | yes   |
        radius2     | radius at point2                  | yes   |
        innerRadius2| inner radius at point2            | no    | 0
    \endtable

Note
    Initial implementation, might suffer from robustness (e.g. radius1==radius2)

Note
    Longer type name : \c searchableCone

SourceFiles
    searchableCone.C

\*---------------------------------------------------------------------------*/

#ifndef searchableCone_H
#define searchableCone_H

#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class searchableCone Declaration
\*---------------------------------------------------------------------------*/

class searchableCone
:
    public searchableSurface
{
    // Private Member Data

        //- 'Left' point
        const point point1_;

        //- Outer radius at point1
        const scalar radius1_;

        //- Inner radius at point1
        const scalar innerRadius1_;


        //- 'Right' point
        const point point2_;

        //- Outer radius at point2
        const scalar radius2_;

        //- Inner radius at point2
        const scalar innerRadius2_;


        //- Length of vector point2-point1
        const scalar magDir_;

        //- Normalised vector point2-point1
        const vector unitDir_;


        //- Names of regions
        mutable wordList regions_;


   // Private Member Functions

        //- Find nearest point on cylinder.
        void findNearestAndNormal
        (
            const point& sample,
            const scalar nearestDistSqr,
            pointIndexHit & nearInfo,
            vector& normal
        ) const;

        //- Determine radial coordinate (squared)
        static scalar radius2(const searchableCone& cone, const point& pt);

        //- Find both intersections with cone. innerRadii supplied externally.
        void findLineAll
        (
            const searchableCone& cone,
            const scalar innerRadius1,
            const scalar innerRadius2,
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;

        //- Insert a hit if it differs (by a tolerance) from the existing ones
        void insertHit
        (
            const point& start,
            const point& end,
            List<pointIndexHit>& info,
            const pointIndexHit& hit
        ) const;

        //- Return the boundBox of the cylinder
        boundBox calcBounds() const;

        //- No copy construct
        searchableCone(const searchableCone&) = delete;

        //- No copy assignment
        void operator=(const searchableCone&) = delete;


public:

    //- Runtime type information
    TypeName("searchableCone");


    // Constructors

        //- Construct from components
        searchableCone
        (
            const IOobject& io,
            const point& point1,
            const scalar radius1,
            const scalar innerRadius1,
            const point& point2,
            const scalar radius2,
            const scalar innerRadius2
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableCone
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableCone() = default;


    // Member Functions

        //- Names of regions
        virtual const wordList& regions() const;

        //- Whether supports volume type (below)
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- What is type of points outside bounds
        virtual volumeType outsideVolumeType() const
        {
            return volumeType::OUTSIDE;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            NotImplemented;
            return false;
        }


        // Multiple point queries.

            //- Find nearest point on cylinder
            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            //- Find nearest intersection on line from start to end
            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Find all intersections in order from start to end
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- Find any intersection on line from start to end
            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point.
            //  Unknown if cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField& points,
                List<volumeType>& volType
            ) const;


        // regIOobject implementation

            virtual bool writeData(Ostream&) const
            {
                NotImplemented;
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

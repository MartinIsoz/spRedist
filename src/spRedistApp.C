/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

Application
    spRedistApp

Description
    Generates sootVolFrac and distributes it
    Note: soot ~ solid phase (sp)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// -- necessary includes
#include "ListOps.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// !!! THIS YOU NEED TO PLUG INTO YOUR CODE //
#include "fluxSootRedistribution.H"
// !!! THIS YOU NEED TO PLUG INTO YOUR CODE (end) //
// Note (MI): this needs to come BEFORE your main function

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);
    
    // - dictionary for app settings
    IOdictionary sootTestAppDict
    (
        IOobject
        (
            "sootTestAppDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    const List<labelList> addCells = sootTestAppDict.lookup("addCells");
    
    scalar addSoot  = readScalar(sootTestAppDict.lookup("addSoot"));
    scalar sootSurplusTol = readScalar(sootTestAppDict.lookup("sootSurplusTol"));
    
    label maxIter   = readLabel(sootTestAppDict.lookup("maxIter"));
    
    // - fields used for computation
    // the main field -> supplied by your program
    volScalarField sootVolFrac
    (
        IOobject
        (
            "sootVolFrac",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sootVolFrac",dimless,0)
    );
    
    // auxiliary volScalarField, it is faster to declare it outside of
    // the main compuational loop. whatever values are supplied, the
    // code does not used them
    volScalarField surPlusSoot               
    (
        IOobject
        (
            "surPlusSoot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("surPlusSoot",dimless,0)
    );
    //~ volScalarField surPlusSoot(sootVolFrac);                        //alternative yet usable definition
    
    // auxiliary surfaceScalarField, it is faster to declare it outside of
    // the main compuational loop. whatever values are supplied, the
    // code does not used them
    surfaceScalarField phiS
    (
        IOobject
        (
            "phiS",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phiS",dimless,0)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // -- prepare objects for soot redistribution
    
    // access preparation (to make the function code better readable)
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbor = mesh.faceNeighbour();
    
    // random number generator (to surpress the soot redistribution
    // directionality
    Random randObj(clock::getTime());
    
    // memory allocation
    // Note (MI): these are passed as non-const references to the
    //            function, which is faster than always allocating the
    //            memory
    // Note (MI): the memory is assumed to be 10% of the corresponding
    //            mesh object. this might have to be changed based
    //            on your actual needs
    DynamicList<label> cellsToTest;
    cellsToTest.reserve(ceil(mesh.nCells()*0.1));
    DynamicList<label> problemCells;
    problemCells.reserve(ceil(mesh.nCells()*0.1));
    DynamicList<label> noFluxFaces;
    noFluxFaces.reserve(ceil(mesh.nFaces()*0.1));
    
    DynamicList<label> cellFluxFaces;                                   //faces in individual cells that are used for flux
    cellFluxFaces.reserve(20);                                          //no cell with more than 20 faces
    
    // -- prepare the algorithm controls
    List<scalar> addedSoot(Pstream::nProcs(),0.0);                      //counter for added soot
    scalar  initialSoot(fvc::domainIntegrate(sootVolFrac).value());
    scalar  discardedSoot(0.0);                                         //counter for discarded soot
    // Note (MI): there are two objects because addedSoot is updated
    //            by adding to each processor separately and then
    //            summing over the processors. In contrast, discardedSoot
    //            is updated via fvc::domainIntegrate, which does the
    //            sync accross the processors
    
    // -- dummies usefull just for this particular test
    scalar meshVol(gSum(mesh.V()));                                     //total mesh volume
    

    Info << "\nStarting the soot re-distribution loop\n" << endl;
    
    scalarList cellV;
    
    
    forAll (addCells,procI)
    {
        if (Pstream::myProcNo() == procI)
        {
            forAll (addCells[procI],cellI)
            {
                label addCell(addCells[procI][cellI]);
                cellV.append(mesh.V()[addCell]);
            }
        }
    }
    
    scalar averageAddV(gAverage(cellV));
    
        
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // prepare list of the cells to check
        cellsToTest.clear();
        
        // add soot to the wanted cells
        
        forAll (addCells,procI)
        {
            if (Pstream::myProcNo() == procI)
            {
                forAll (addCells[procI],cellI)
                {
                    label addCell(addCells[procI][cellI]);
                    sootVolFrac[addCell] += addSoot*averageAddV/mesh.V()[addCell];
                    addedSoot[Pstream::myProcNo()] += addSoot*averageAddV;
                    cellsToTest.append(addCell);
                }
            }
        }

        
        // !!! THIS YOU NEED TO PLUG INTO YOUR CODE //
        discardedSoot += fluxRedistributeSoot(
            sootVolFrac,
            surPlusSoot,
            phiS,
            sootSurplusTol,
            maxIter,
            cellsToTest,
            problemCells,
            noFluxFaces,
            cellFluxFaces,
            randObj,
            mesh,
            faceOwner,
            faceNeighbor
        );
        // !!! THIS YOU NEED TO PLUG INTO YOUR CODE (end) //
        
        // -- potentially usefull control outputs
        scalar sumSootVolFrac(fvc::domainIntegrate(sootVolFrac).value());
        sumSootVolFrac -= initialSoot;     
        
        Info << "added   soot mass  : " << gSum(addedSoot) << endl;
        Info << "current soot mass  : " << sumSootVolFrac << endl;
        Info << "current mass loss  : " << (gSum(addedSoot) - sumSootVolFrac)/(gSum(addedSoot))*100.0 << "%" << endl;
        Info << "max sootVolFrac    : " << gMax(sootVolFrac) << endl;
        Info << "total discard. soot: " << discardedSoot << " ~ " << discardedSoot/gSum(addedSoot)*100.0 << " % of the added mass" << endl;
        
        // -- dummy output usefull for this particular test
        if (gSum(addedSoot) > meshVol)
        {
            Info << "!! no more soot can be redistributed !!" << endl;
            maxIter = -1;
        }
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

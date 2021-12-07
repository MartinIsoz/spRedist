scalar fluxRedistributeSoot
(
    volScalarField& sootVolFrac,
    volScalarField& surPlusSoot,
    surfaceScalarField& phiS,
    scalar& sootSurplusTol,
    label& maxIter,
    DynamicList<label>& cellsToTest,
    DynamicList<label>& problemCells,
    DynamicList<label>& noFluxFaces,
    DynamicList<label>& cellFluxFaces,
    Random& randObj,
    const fvMesh& mesh,
    const labelList& faceOwner,
    const labelList& faceNeighbor
)
{
    // -- prepare inputs  
    bool iterate(true);label itCnt(0);
    label cellI;
    
    surPlusSoot *= 0.0;
    phiS *= 0.0;
    
    // -- initialize the tick-tock face blocking
    bool tick(true);
    
    // -- descend into the iteration loop
    while (iterate)
    {
        
        // - just auxiliaries for output
        label nCellsToTest(cellsToTest.size());
        reduce(nCellsToTest, sumOp<label>());
        Info << "iteration " << itCnt << " cells to test: "  << nCellsToTest << endl;
        
        // - clearing the lists
        noFluxFaces.clear();
        problemCells.clear();
        
        // - localize the soot surPlus
        if (tick)
        {
            forAll(cellsToTest,tCellI)                                  //go only through the cells to be tested
            {
                cellI = cellsToTest[tCellI];
                surPlusSoot[cellI] = max(0.0,sootVolFrac[cellI]-1.0);
                if (surPlusSoot[cellI] > sootSurplusTol) 
                {
                    problemCells.append(cellI);
                    // go through the cell neighbors and check if they are full
                    const cell& cellFaces = mesh.cells()[cellI];
                    forAll (cellFaces,faceI)
                    {
                        if (mesh.isInternalFace(cellFaces[faceI]))
                        {
                            if
                            (
                                (sootVolFrac[faceNeighbor[cellFaces[faceI]]] - (1.0 - SMALL) > 0.0)
                                and
                                (sootVolFrac[faceOwner[cellFaces[faceI]]]    - (1.0 - SMALL) > 0.0)
                            )
                            {
                                noFluxFaces.append(cellFaces[faceI]);
                            }
                        }
                        //~ else
                        //~ {
                            //~ label facePatchI = mesh.boundaryMesh().whichPatch(cellFaces[faceI]);
                            //~ const polyPatch& patch = mesh.boundaryMesh()[facePatchI];
                            //~ if (patch.coupled())
                            //~ {
                                //~ const scalarField& sootVolFracRed
                                //~ (
                                    //~ sootVolFrac.boundaryFieldRef()[facePatchI].patchNeighbourField()
                                //~ );
                                //~ if (sootVolFracRed[cellFaces[faceI]-patch.start()] - (1.0 - SMALL) > 0.0)
                                //~ {
                                    //~ noFluxFaces.append(cellFaces[faceI]);
                                    //~ Pout << "HIT BOUNDARY" << endl;
                                //~ }
                            //~ }
                        //~ }
                        //~ // THERE IS SOMETHING ROTTEN IN THE STATE OF DENMARK
                    }
                }
            }
            tick = false;
        }
        else
        {
            forAll(cellsToTest,tCellI)                                  //go only through the cells to be tested
            {
                cellI = cellsToTest[tCellI];
                surPlusSoot[cellI] = max(0.0,sootVolFrac[cellI]-1.0);
                if (surPlusSoot[cellI] > sootSurplusTol) 
                {
                    problemCells.append(cellI);
                }
            }
            tick = true;
        }
        
        
        // - shuffle to surpress directionality
        for (label posi = problemCells.size()-1; posi > 0; --posi)
        {
            const label i = randObj.sampleAB<label>(0, posi);
            Foam::Swap(problemCells[i], problemCells[posi]);
        }
        
        // - create soot fluxes accross the cell faces
        forAll(problemCells,pCellI)
        {
            cellI = problemCells[pCellI];
            scalar totalTransferArea(0.0);                          //total area available to transfer
            
            label nTransferFaces(0);
            const cell& cellFaces = mesh.cells()[cellI];

            // - go through cellFaces and see which are available
            //   for soot transfer
            // Note (MI): transfer soot only via proc-proc boundaries
            cellFluxFaces.clear();
                
            forAll(cellFaces,faceI)
            {
                labelList minIndLs = Foam::findIndices(noFluxFaces,cellFaces[faceI]);
                if (minIndLs.size() == 0)
                {
                    if (mesh.isInternalFace(cellFaces[faceI])) 
                    {
                            nTransferFaces++;
                            totalTransferArea += mesh.magSf()[cellFaces[faceI]];
                            cellFluxFaces.append(cellFaces[faceI]);
                    }
                    else
                    {
                        label facePatchI = mesh.boundaryMesh().whichPatch(cellFaces[faceI]);
                        const polyPatch& patch = mesh.boundaryMesh()[facePatchI];
                        if (patch.coupled())
                        {
                            nTransferFaces++;
                            totalTransferArea += mesh.magSf().boundaryField()[facePatchI][cellFaces[faceI]-patch.start()];
                            cellFluxFaces.append(cellFaces[faceI]);
                        }
                    }
                }
            }
            
            // - modify the current cells by the soot outflux
            // Note (MI): the soot fluxes are always positive
            //            thus, 1.0*surPlusSoot[cellI] will be
            //            added back through the fluxes
            if (nTransferFaces > label(0))
            {
                sootVolFrac[cellI] -= 2.0*surPlusSoot[cellI];
                surPlusSoot[cellI] *= mesh.V()[cellI]/(totalTransferArea);
            }
            
            // - go through the faces and assign the soot fluxes
            //   to the faces
            forAll(cellFluxFaces,cFFaceI)
            {
                label faceI(cellFluxFaces[cFFaceI]);
                // - nothing special for internal cells
                if (mesh.isInternalFace(faceI))
                {
                    phiS[faceI] = surPlusSoot[cellI]*mesh.magSf()[faceI];
                    noFluxFaces.append(faceI);
                }
                // a bit more fun for boundaries
                else
                {
                    label facePatchId = mesh.boundaryMesh().whichPatch(faceI);
                    const polyPatch& patch = mesh.boundaryMesh()[facePatchId];
                    label facePatchI = faceI - patch.start();
                    if (patch.coupled())                        //work only with coupled patches
                    {
                        phiS.boundaryFieldRef()[facePatchId][facePatchI] = surPlusSoot[cellI]*mesh.magSf().boundaryField()[facePatchId][facePatchI];
                        noFluxFaces.append(faceI);
                    }
                }
            }
        }
        
        // Note (MI): just a sync dummy
        surPlusSoot.correctBoundaryConditions();
        sootVolFrac.correctBoundaryConditions();
        
        // Note (MI): now, you should have
        // - soot fluxes across all the internal faces
        // - soot fluxes across all BUT the wall boundary faces
        // - no coupled boundaries are synced at the moment !!!
        
        // - sync the processor boundaries
        // Note (MI): I need to do SUM operation and rewrite the
        //            boundary field at the given patches
        // Note (MI): the HARD step: identify the current coupled
        //            patch neighbor and sync the data - CIRCUMVENTED
        
        forAll(mesh.boundaryMesh(),patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];
            if (patch.coupled())
            {
                const scalarField& phiSRed
                (
                    surPlusSoot.boundaryFieldRef()[patchI].patchNeighbourField()
                );
                forAll(patch,faceI)
                {
                    phiS.boundaryFieldRef()[patchI][faceI] += phiSRed[faceI]*mesh.magSf().boundaryField()[patchI][faceI];
                }
            }
        }
        
        // - apply the fluxes
        // Note (MI): the idea is that at this time, all should be
        //            synced accross the coupled faces
        // Note (MI): if something happens on the face, add the
        //            adjacent cells to the list of cells to be
        //            tested (should communicate over proc-proc
        //            boundaries)
        forAll(phiS.internalField(),faceI)
        {
            if (phiS.internalField()[faceI]/mesh.magSf()[faceI] > VSMALL)//assume magSf > 0 for all the mesh faces
            {
                cellsToTest.append(faceOwner[faceI]);
                cellsToTest.append(faceNeighbor[faceI]);
                sootVolFrac[faceOwner[faceI]]    += phiS.internalField()[faceI]/mesh.V()[faceOwner[faceI]];
                sootVolFrac[faceNeighbor[faceI]] += phiS.internalField()[faceI]/mesh.V()[faceNeighbor[faceI]];
            }
        }
        forAll(phiS.boundaryField(),patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];
            forAll(patch,faceI)
            {
                if (phiS.boundaryFieldRef()[patchI][faceI]/mesh.magSf().boundaryField()[patchI][faceI] > VSMALL)//assume magSf > 0 for all the mesh faces
                {
                    cellsToTest.append(faceOwner[faceI+patch.start()]);
                    sootVolFrac[faceOwner[faceI+patch.start()]] += phiS.boundaryFieldRef()[patchI][faceI]/mesh.V()[faceOwner[faceI+patch.start()]];
                }
            }
        }
        // Note (MI): the division by face area is quite necessary
        //            for the computations on catalysts
        //            (extremely small cells)
        
        phiS *= 0.0;
        surPlusSoot.correctBoundaryConditions();
        sootVolFrac.correctBoundaryConditions();
        
        // - make the list cellsToTest unique
        List<label> cellsToTestOrder;
        Foam::uniqueOrder(cellsToTest,cellsToTestOrder);
        const label lenCellsToTestOrder = cellsToTestOrder.size();
        List<label> auxLst(lenCellsToTestOrder);
        
        for (label i=0; i < lenCellsToTestOrder; ++i)
        {
            auxLst[i] = std::move(cellsToTest[cellsToTestOrder[i]]);
        }
        cellsToTest.transfer(auxLst);
    
        // - increase iterator and check for the current problem status
        itCnt++;
        if (itCnt > maxIter or problemCells.size() < 1){iterate=false;}
        reduce(iterate, orOp<bool>());
        
        // Note (MI): just a sync dummy
        surPlusSoot.correctBoundaryConditions();
        sootVolFrac.correctBoundaryConditions();
        
    }
    
    // - finish up computations
    surPlusSoot = max(0.0,sootVolFrac-1.0);
    surPlusSoot.correctBoundaryConditions();
    scalar stepDiscardedSoot = fvc::domainIntegrate(surPlusSoot).value();
    
    // - discard the soot surplus
    sootVolFrac -= surPlusSoot;
    
    return stepDiscardedSoot;
}

/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTMAKEANDCLOSEWOUND_HPP_
#define TESTMAKEANDCLOSEWOUND_HPP_
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */
#include "PetscSetupAndFinalize.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "Cell.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "FarhadifarForce.hpp"
#include "SmartPointers.hpp"
#include "SimpleTargetAreaModifier.hpp"

class TestWoundHealing : public AbstractCellBasedWithTimingsTestSuite
{
public:
    void TestMakeAndCloseWound()
    {
        // Create a simple 2D MutableVertexMesh
        // MutableVertexMesh<2,2>* p_mesh;
        VoronoiVertexMeshGenerator mesh_generator(10,10,5,1.0); // number of cells, number of cells, number of LLoyds steps, elementArea

        // Get the centroid of the mesh
//        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();
        Toroidal2dVertexMesh* p_mesh = mesh_generator.GetToroidalMesh();
        c_vector<double, 2> mesh_centroid = zero_vector<double>(2);
        for (unsigned element_index = 0; element_index < p_mesh->GetNumAllElements(); element_index++)
        {
            mesh_centroid += p_mesh->GetCentroidOfElement(element_index);
        }
        mesh_centroid /= p_mesh->GetNumAllElements();

        // delete all elements in the mesh that fall within a circle of radius 5
        c_vector<double, 2> this_centroid = zero_vector<double>(2);

        for (unsigned element_index = 0; element_index < p_mesh->GetNumAllElements(); element_index++)
        {
            this_centroid = p_mesh->GetCentroidOfElement(element_index);
            if (norm_2(this_centroid - mesh_centroid) < 1.5)
            {
                p_mesh->DeleteElementPriorToReMesh(element_index);
            }
        }

        p_mesh->ReMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-(double)i -19.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestWoundHealing");
        simulator.SetEndTime(20.0);
        simulator.SetSamplingTimestepMultiple(10);
        // Create a force law and pass it to the simulation
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetBoundaryLineTensionParameter(1.0);
        simulator.AddForce(p_force);

        // A FarhadifarForce has to be used together with an AbstractTargetAreaModifier#2488
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.SetDt(0.01);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTHELLO_HPP_*/

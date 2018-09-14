/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "WoundHealingForce.hpp"

template<unsigned DIM>
WoundHealingForce<DIM>::WoundHealingForce()
   : AbstractForce<DIM>(),
     mWoundTensionParameter(0.12) // this parameter as such does not exist in Farhadifar's model.
{
}

template<unsigned DIM>
WoundHealingForce<DIM>::~WoundHealingForce()
{
}

template<unsigned DIM>
void WoundHealingForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    // First, find the first boundary node
    unsigned first_boundary_node_start = 0;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        if ( p_this_node->IsBoundaryNode() )
        {
            first_boundary_node_start = node_index;
            break;
        }
    }

    // Then a while loop to collect the nodes along the boundary
    std::vector<unsigned> first_boundary_nodes;
    first_boundary_nodes.push_back(first_boundary_node_start);
    bool boundary_incomplete = true;
    unsigned current_boundary_node = first_boundary_node_start;
    while (boundary_incomplete)
    {
        // find adjacent cells A and B
        std::set<unsigned> containing_elem_indices = p_cell_population->
                GetNode(current_boundary_node)->rGetContainingElementIndices();
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned number_of_nodes_in_this_element = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(current_boundary_node);
            unsigned next_node_local_index = (local_index+1)%number_of_nodes_in_this_element;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);
            unsigned next_node_index = p_next_node->GetIndex();
            if (next_node_index == first_boundary_node_start)
            {
                boundary_incomplete = false;
                break;
            }
            else if (p_next_node->IsBoundaryNode())
            {
                first_boundary_nodes.push_back( p_next_node->GetIndex() );
                current_boundary_node = p_next_node->GetIndex();
                break;
            }
        } // element for loop
    }// while loop

    // find a node on the next boundary
    unsigned second_boundary_node_start = 0;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        if ( p_this_node->IsBoundaryNode() )
        {
            if(std::find(first_boundary_nodes.begin(),
                    first_boundary_nodes.end(), node_index) == first_boundary_nodes.end())
            {
                second_boundary_node_start = node_index;
                break;
            }
        }
    }

    // Then a while loop to collect the nodes along the boundary
    std::vector<unsigned> second_boundary_nodes;
    boundary_incomplete = true;
    second_boundary_nodes.push_back(second_boundary_node_start);
    current_boundary_node = second_boundary_node_start;
    while (boundary_incomplete)
    {
        // find adjacent cells A and B
        std::set<unsigned> containing_elem_indices = p_cell_population->
                GetNode(current_boundary_node)->rGetContainingElementIndices();
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned number_of_nodes_in_this_element = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(current_boundary_node);
            unsigned next_node_local_index = (local_index+1)%number_of_nodes_in_this_element;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);
            unsigned next_node_index = p_next_node->GetIndex();
            if (next_node_index == second_boundary_node_start)
            {
                boundary_incomplete = false;
                break;
            }
            else if (p_next_node->IsBoundaryNode())
            {
                second_boundary_nodes.push_back(p_next_node->GetIndex());
                current_boundary_node = p_next_node->GetIndex();
                break;
            }
        } // element for loop
    }// while loop

    // Then, identify inner boundary

    std::vector<unsigned> inner_boundary_nodes;
    if ( first_boundary_nodes.size() < second_boundary_nodes.size() )
    {
        inner_boundary_nodes = first_boundary_nodes;
    }
    else
    {
        inner_boundary_nodes = second_boundary_nodes;
    }

    // Then, add forces to that boundary a for loop

    for ( auto &node_index : inner_boundary_nodes )
    {
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        c_vector<double, DIM> line_tension_contribution = zero_vector<double>(DIM);
        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's area elasticity (note the minus sign)

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            if ( p_previous_node->IsBoundaryNode() )
            {
                c_vector<double, DIM> previous_edge_gradient =
                        -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element,
                                previous_node_local_index);
                line_tension_contribution -= mWoundTensionParameter*previous_edge_gradient;
            }

            if ( p_next_node->IsBoundaryNode() )
            {
                c_vector<double, DIM> next_edge_gradient = p_cell_population->
                        rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);
                line_tension_contribution -= mWoundTensionParameter*next_edge_gradient;
            }

        }
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(line_tension_contribution);
    }

}

template<unsigned DIM>
double WoundHealingForce<DIM>::GetWoundTensionParameter()
{
    return mWoundTensionParameter;
}

template<unsigned DIM>
void WoundHealingForce<DIM>::SetWoundTensionParameter(double woundTension)
{
    mWoundTensionParameter = woundTension;
}

template<unsigned DIM>
void WoundHealingForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<WoundTensionParameter>" << mWoundTensionParameter << "</WoundTensionParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class WoundHealingForce<1>;
template class WoundHealingForce<2>;
template class WoundHealingForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WoundHealingForce)

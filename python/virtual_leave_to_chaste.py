import untangle

def convert_xml_to_chaste(xml_path):
    xml_data = untangle.parse(xml_path)

    nodes = xml_data.leaf.nodes.node

    x_coordinates = [ float(node['x']) for node in nodes]
    y_coordinates = [ float(node['y']) for node in nodes]
    boundary_list = [ node['boundary'] == 'true' for node in nodes]
    
    cells = xml_data.leaf.cells.cell
    
    cell_entries = []
    for cell in cells:
        cell_nodes = cell.node
        these_indices = [ int(node['n']) for node in cell_nodes ]
        cell_entries.append(these_indices)
        
    filtered_node_indices = []
    
    for node_index in range(len(x_coordinates)):
        number_of_cells_found = 0
        for cell in cell_entries:
            if node_index in cell:
                number_of_cells_found +=1
        if number_of_cells_found > 2 or boundary_list[node_index]:
            filtered_node_indices.append(node_index)
            
    print cell_entries

    new_cell_entries = []
    for cell in cell_entries:
        new_cell = []
        for node_index in cell:
            if node_index in filtered_node_indices:
                new_cell.append(node_index)
        new_cell_entries.append(new_cell)
        
    print new_cell_entries
    
    nodes_file = open('virtual_leaf.node','w')
    
    nodes_file.write(str(len(filtered_node_indices)) + '\t2\t0\t1\n')
    for node_index in filtered_node_indices:
        nodes_file.write(str(node_index) + '\t' + 
                   str(x_coordinates[node_index]) + '\t' +
                   str(y_coordinates[node_index]) + '\t' +
                   str(int(boundary_list[node_index])) + '\n' )

    nodes_file.close()

    elements_file = open('virtual_leaf.cell','w')
    elements_file.write(str(len(new_cell_entries)) + '\t1\n')
    for cell_index, cell_entries in enumerate(new_cell_entries):
        elements_file.write(str(cell_index) + '\t' + 
                   str(len(cell_entries)))
        for node_index in cell_entries:
            elements_file.write( '\t' + str(node_index) )
        elements_file.write( '\t0\n' )

    elements_file.close()

if __name__ == "__main__":
    path = 'data/virtual_leaf_example.xml'
    convert_xml_to_chaste(path)
#     with open(path) as fd:
#     doc = xmltodict.parse(fd.read())

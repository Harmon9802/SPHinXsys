/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
 /**
  * @file 	MeshHelper.h
  * @brief 	Here, we define the functions used for reading the mesh data.
  * @author	Yash Mandaokar Zhentong Wang and Xiangyu Hu
  */
#ifndef MESHTYPE_H
#define MESHTYPE_H

#include "unstructured_mesh.h"
using namespace std;
namespace SPH
{

    class MeshFileHelpers 
    {
        public:

        MeshFileHelpers()
        {
     
        };
        virtual ~MeshFileHelpers(){};

        static void meshDimension(ifstream& mesh_file, size_t& dimension, string& text_line);
        static void numberofNodes(ifstream& mesh_file, size_t& number_of_points, string& text_line);
        static void nodeCoordinates(ifstream& mesh_file, StdLargeVec<Vecd>& node_coordinates_, string& text_line,size_t& dimension);
        static void numberofElements(ifstream& mesh_file, size_t& number_of_elements, string& text_line);
        static void dataStruct(vector<vector<vector<size_t>>>& mesh_topology_, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_,
                                size_t number_of_elements, size_t mesh_type, size_t dimension);
        static size_t findBoundaryType(string& text_line, size_t boundary_type);
        static Vecd nodeIndex(string& text_line);
        static Vec2d cellIndex(string& text_line);
        static void updateElementsNodesConnection(StdLargeVec<StdVec<size_t>>& elements_nodes_connection_, Vecd nodes, Vec2d cells,
                                                    bool& check_neighbor_cell1, bool& check_neighbor_cell2);
        static void updateCellLists(vector<vector<vector<size_t>>>& mesh_topology_, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_, Vecd nodes,
                                    Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type);
        static void updateBoundaryCellLists(vector<vector<vector<size_t>>>& mesh_topology_, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_,
                                            Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type);
        static void cellCenterCoordinates(StdLargeVec<StdVec<size_t>>& elements_nodes_connection_, std::size_t& element,
                                        StdLargeVec<Vecd>& node_coordinates_, StdLargeVec<Vecd>& elements_center_coordinates_, Vecd& center_coordinate);
        static void elementVolume(StdLargeVec<StdVec<size_t>>& elements_nodes_connection_, std::size_t& element,
                                    StdLargeVec<Vecd>& node_coordinates_, StdLargeVec<Real>& elements_volumes_);
        static void minimumdistance(vector<Real>& all_data_of_distance_between_nodes, StdLargeVec<Real>& elements_volumes_, 
                                    vector<vector<vector<size_t>>>& mesh_topology_, StdLargeVec<Vecd>& node_coordinates_);
        static void vtuFileHeader(std::ofstream& out_file);
        static void vtuFileNodeCoordinates(std::ofstream& out_file, StdLargeVec<Vecd>& nodes_coordinates_, 
                                    StdLargeVec<StdVec<size_t>>& elements_nodes_connection_, SPHBody& bounds_,Real& rangemax);
        static void vtuFileInformationKey(std::ofstream& out_file, Real& rangemax);
        static void vtuFileCellConnectivity(std::ofstream& out_file, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_, StdLargeVec<Vecd>& node_coordinates_);
        static void vtuFileOffsets(std::ofstream& out_file, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_);
        static void vtuFileTypeOfCell(std::ofstream& out_file, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_);

        /*Functions for .msh file from FLUENT*/
        static void numberofNodesFluent(ifstream& mesh_file, size_t& number_of_points, string& text_line);
        static void numberofElementsFluent(ifstream& mesh_file, size_t& number_of_elements, string& text_line);
        static void nodeCoordinatesFluent(ifstream& mesh_file, StdLargeVec<Vecd>& node_coordinates_, string& text_line,
                                            size_t& dimension);
        static void updateBoundaryCellListsFluent(vector<vector<vector<size_t>>>& mesh_topology_, StdLargeVec<StdVec<size_t>>& elements_nodes_connection_,
                                            Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type);
    };


}
#endif // MESHTYPE_H
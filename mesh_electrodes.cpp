/*
 * This program takes in electrode data generated from the svg file containing the HOA trap geometry
 * and converts it to a coordinate mesh for use in finite element analysis. An electrode, for our purposes,
 * is one or more regions of the trap that we can control with a single voltage source.  Electrodes with
 * multiple "patches" are typically distributed in various symmetric regions about the trap. See RS1096.svg
 * for an image of where the electrodes are located and how they are distributed.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Delaunay_mesh_face_base_2<Kernel> MeshFaceBase;
typedef CGAL::Triangulation_data_structure_2<VertexBase, MeshFaceBase> TriangulationData;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationData> CDTriangulation;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDTriangulation> MeshSizeCriteria;

typedef CGAL::Exact_predicates_exact_constructions_kernel IntKernel;
typedef IntKernel::Point_2 Point_2;
typedef CGAL::Polygon_2<IntKernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<IntKernel> Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2> Pwh_list_2;

typedef CGAL::Cartesian_converter<IntKernel, Kernel> Converter;


// Triple is a bundle of three elements of data type T.
template<typename T>
struct Triple
{
    T a, b, c;
};

// A Vertex is a point with three double coordinates.
typedef Triple<double> Vertex;
// A Triangle is a set of three vertices.
typedef Triple<size_t> Triangle;

// Each electrode has a name, a set of triangles that defines its surface,
// and the vertices corresponding to the coordinates of the points making
// up the triangles.
struct Electrode
{
    std::string name;
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
};

void generate_mesh(std::ifstream &electrode_data, std::ofstream &mesh_data);
Polygon_2 create_region_of_interest(std::string region_name);
CDTriangulation triangulate_region_with_boundary(Polygon_2 region_boundary);

int main(int argc, char **argv)
{
    // Notify the user of proper program usage if they fail to provide the correct arguments.
    if (argc != 3)
    {
        std::cout << "Usage " << argv[0] << " [electrode_input_data_file] [mesh_output_data_file] " << std::endl;
        return -1;
    }

    // Attempt to open the input data file.
    std::ifstream electrode_data_file(argv[1]);
    if(!electrode_data_file)
    {
        std::cout << "Invalid input file." << std::endl;
        return -1;
    }

    // Attempt to open the output data file.
    std::ofstream mesh_data_file(argv[2]);
    if(!mesh_data_file)
    {
        std::cout << "Invalid output file." << std::endl;
        return -1;
    }

    // If everything is peachy with the files, generate the mesh data.
    generate_mesh(electrode_data_file, mesh_data_file);

    return 0;
}

void generate_mesh(std::ifstream &electrode_data, std::ofstream &mesh_data)
{
    // Some dynamic storage for the electrodes.
    std::vector<Electrode> electrodes;

    // Get the region of interest over which we would like to generate a mesh.
    Polygon_2 region_of_interest = create_region_of_interest("Region 1");

    // Placeholder for the electrode names.
    std::string electrode_name;

    // As long as there is more data
    while (electrode_data >> electrode_name)
    {
        // Create a new electrode struct at the end of the vector.
        electrodes.push_back(Electrode());
        // Grab a reference to the electrode.
        Electrode &electrode = electrodes.back();
        // Set the electrode's name
        electrode.name = electrode_name;

        // Allocate some space for the number of patches that make up the electrode.
        size_t number_of_patches;
        // Create a variable to hold the z values of the patches. z is defined separately from x and y because an
        // electrode patch is at a uniform z level, so the value is repeated for each corner in the patch.
        double z = 0.0;
        // Grab the number of patches from the file.
        electrode_data >> number_of_patches;
        // Let the user know which electrode we're processing.
        std::cout << "Processing electrode : " << electrode_name << std::endl;

        // Loop over all the patches that make up the electrode.
        for (size_t patch_number = 0; patch_number < number_of_patches; ++patch_number)
        {
            // Let the user know which patch we're currently processing.
            std::cout << "Processing  electrode patch : " << patch_number << std::endl;

            // Grab the number of corners that make up the patch from the data.
            size_t number_of_corners;
            electrode_data >> number_of_corners;
            // We now want to generate a 2D CGAL polygon from our patch data.
            Polygon_2 patch_polygon;
            // Since our patches are closed paths, the last corner and the first corner are the same point,
            // so we loop over the all corners except the last.
            for (size_t corner_number = 0; corner_number < number_of_corners - 1; ++corner_number)
            {
                // Grab the coordinate data from the file.
                double x, y;
                electrode_data >> x >> y >> z;
                // We then ignore the z data for now, since it's the same for all corners in the patch.
                patch_polygon.push_back(Point_2(x, y));
            }
            // Pull the repeated corner data from the file.
            double repeat_x, repeat_y;
            electrode_data >> repeat_x >> repeat_y >> z;

            // We now want a list of the patch regions that intersect the region of interest
            Pwh_list_2 intersecting_regions;
            // This method checks for the intersection of the two polygons and adds them to the list.
            std::cout << patch_polygon << std::endl;
            std::cout << region_of_interest << std::endl;
            CGAL::intersection(region_of_interest, patch_polygon, std::back_inserter(intersecting_regions));
            std::cout << "Intersections checked" << std::endl;
            // If the patch doesn't overlap the region of interest, we don't generate a mesh for it.
            if (!intersecting_regions.size()) continue;

            // Let the user know if the electrode patch does overlap the region of interest.
            std::cout << "Overlap with region of interest found on "
            << electrode_name << ":" << patch_number << std::endl;
            std::cout << intersecting_regions.size() << " polygons" << std::endl;

            // Now loop over every region in the list.
            for (Pwh_list_2::iterator region = intersecting_regions.begin();
                 region != intersecting_regions.end(); ++region)
            {
                // Get the polygon representation of the overlapping region.
                const Polygon_2 &region_boundary = region->outer_boundary();
                // Let the user know how many boundary points (polygon corners) there are.
                std::cout << region_boundary.size() << " boundary points" << std::endl;

                // Generate the triangulated mesh from the boundary of the polygon.
                CDTriangulation triangulation = triangulate_region_with_boundary(region_boundary);

                // Let the user know how many vertices are in the mesh.
                std::cout << triangulation.number_of_vertices() << " mesh vertices" << std::endl;

                // Create a map to index the vertices is the mesh.
                std::map<CDTriangulation::Vertex_handle, size_t> indices;

                // Loop over all the vertices in the mesh.
                for (CDTriangulation::Vertex_iterator mesh_vertex = triangulation.vertices_begin();
                     mesh_vertex != triangulation.vertices_end(); ++mesh_vertex)
                {
                    // Convert from the CGAL vertex handle to our own vertex data structure.
                    Vertex vert;
                    vert.a = mesh_vertex->point().x();
                    vert.b = mesh_vertex->point().y();
                    vert.c = z;
                    // Add the vertex to the current electrode.
                    electrode.vertices.push_back(vert);
                    // Assign the next index to the processed vertex.
                    indices[mesh_vertex] = electrode.vertices.size() - 1;
                }

                // Grab the first face in the generated mesh.
                CDTriangulation::Finite_faces_iterator mesh_face = triangulation.finite_faces_begin();

                // Loop over all faces in the mesh.
                for (mesh_face; mesh_face != triangulation.finite_faces_end(); ++mesh_face)
                {
                    // Sanity check, the mesh should be in our computational domain, otherwise skip it.
                    if (!mesh_face->is_in_domain()) continue;

                    // Second sanity check, make sure that at least one of the vertices of
                    // the mesh face showed up as a separate vertex of the triangulation.
                    if (!indices.count(mesh_face->vertex(0))
                        || !indices.count(mesh_face->vertex(1))
                        || !indices.count(mesh_face->vertex(2)))
                    {
                        // If none of the vertices of the mesh face showed up, tell user know there was a bad triangle.
                        std::cout << "BAD TRI: Bad " << electrode_name << ", "
                        << indices[mesh_face->vertex(0)] << ":"
                        << indices[mesh_face->vertex(1)] << ":"
                        << indices[mesh_face->vertex(2)] << std::endl;
                        // Then move on to the next mesh face.
                        continue;
                    }

                    // The triangles of our mesh consist of the indexes associated with the vertices of the face.
                    // Take the combined information about the mesh faces and the vertices and coerce it into
                    // our own triangle data structure.
                    Triangle triangle;
                    triangle.a = indices[mesh_face->vertex(0)];
                    triangle.b = indices[mesh_face->vertex(1)];
                    triangle.c = indices[mesh_face->vertex(2)];
                    // Add the triangle to the current electrode.
                    electrode.triangles.push_back(triangle);
                }
            }
        }

        // Make sure the current electrode has a mesh associated with it.
        if (!electrode.vertices.size() || !electrode.triangles.size())
        {
            // If not, let the user know.
            std::cout << "No mesh associated with electrode : " << electrode_name << std::endl;
            continue;
        }

        // Write the electrode name and number of vertices to the data file.
        mesh_data << "ELEC " << electrode.name << std::endl;
        mesh_data << "VERTS " << electrode.vertices.size() << std::endl;
        // For every vertex in the electrode mesh
        for (size_t vertex_index = 0; vertex_index < electrode.vertices.size(); ++vertex_index)
        {
            // Grab a reference to the vertex associated with the index.
            Vertex &vertex = electrode.vertices[vertex_index];
            // Write the vertex coordinate information as a line in the mesh file.
            mesh_data << vertex.a << " " << vertex.b << " " << vertex.c << std::endl;
        }

        // Write the number of faces associated with the mesh to the file
        mesh_data << "FACES " << electrode.triangles.size() << std::endl;
        // For every triangle in the electrode mesh
        for (size_t triangle_index = 0; triangle_index < electrode.triangles.size(); ++triangle_index)
        {
            // Grab a reference to the triangle associated with the triangle index.
            Triangle &triangle = electrode.triangles[triangle_index];
            // Write the triangle coordinate information as a line in the mesh file.
            mesh_data << triangle.a << " " << triangle.b << " " << triangle.c << std::endl;
        }
    }
}

Polygon_2 create_region_of_interest(std::string region_name)
{
    // Scaling is still pretty magical
    double scaling = 100 / 9.;
    std::cout << (region_name == "Region 1") << std::endl;
    if (region_name == "Region 1")
    {
        double x1 = scaling * 80;  // 888.8
        double x2 = scaling * 160; // 1777.7
        double y1 = scaling * 140; // 1555.5
        double y2 = scaling * 75;  // 833.3
        Polygon_2 region_of_interest;
        region_of_interest.push_back(Point_2(x1, y1));
        region_of_interest.push_back(Point_2(x1, y2));
        region_of_interest.push_back(Point_2(x2, y2));
        region_of_interest.push_back(Point_2(x2, y1));

        return region_of_interest;
    }
}

CDTriangulation triangulate_region_with_boundary(Polygon_2 region_boundary)
{
    // The converter allows us to switch to a kernel where doubles are used for the geometric
    // constructions (the meshes to be generated from our coordinate data).  This is faster but
    // suffers the normal pitfalls of computational arithmetic.
    Converter converter;

    // We use a Delaunay triangulation to avoid sliver triangles in the mesh.
    CDTriangulation triangulation;

    // Grab an iterator on the boundary vertices of the intersection region.
    Polygon_2::Vertex_iterator boundary_vertex = region_boundary.vertices_begin();
    // Create some vertex handles for insterting constraints into the triangulation (ie making sure
    // the boundary is included as part of the mesh).
    CDTriangulation::Vertex_handle start_vertex, current_vertex, previous_vertex;

    // Keep a handle to the first vertex in the boundary, converted and inserted into the triangulation.
    start_vertex = triangulation.insert(converter(*boundary_vertex));

    // Fencepost assignment to start_vertex the boundary definition in the triangulation.
    previous_vertex = start_vertex;
    // Skip to the second vertex and loop over the remaining vertices in the boundary.
    for (++boundary_vertex; boundary_vertex != region_boundary.vertices_end(); ++boundary_vertex)
    {
        // Convert the vertex, insert it into the triangulation, and grab a handle to it.
        current_vertex = triangulation.insert(converter(*boundary_vertex));
        // Set the segment between the previous and current points
        // as a constraint (boundary) of the triangulation.
        triangulation.insert_constraint(previous_vertex, current_vertex);
        // Update the previous vertex handle.
        previous_vertex = current_vertex;
    }
    // Provide the final boundary segment to the triangulation.
    triangulation.insert_constraint(previous_vertex, start_vertex);

    // Specify the criteria for the mesh generation.  The shape criteria is the default (corresponding
    // a minimum angle in a given triangle of ~20.7 degrees and guaranteeing convergence of the
    // algorithm).  The size criteria bounds the maximum length of a triangle edge.
    // TODO: Still, with the scaling, figure out why the value 15 was chosen, and what it corresponds to.
    double shape_criteria = 0.125;
    double size_criteria = 15.0;

    // Now generate the mesh for the given criteria.
    CGAL::refine_Delaunay_mesh_2(triangulation, MeshSizeCriteria(shape_criteria, size_criteria));

    return triangulation;
}

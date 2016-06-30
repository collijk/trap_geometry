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
typedef std::list<Kernel::Segment_2> ConstraintList;

typedef CGAL::Exact_predicates_exact_constructions_kernel IntKernel;
typedef IntKernel::Point_2 Point_2;
typedef CGAL::Polygon_2<IntKernel> Polygon_2;
typedef CGAL::Polygon_with_holes_2<IntKernel> Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2> Pwh_list_2;

typedef CGAL::Cartesian_converter<IntKernel, Kernel> Converter;


// Triple is a bundle of three elements of data type T.
template<typename T>
struct Triple {
    T a, b, c;
};

// A Vertex is a point with three double coordinates.
typedef Triple<double> Vertex;
// A Triangle is a set of three vertices.
typedef Triple<size_t> Triangle;

// Each electrode has a name, a set of triangles that defines its surface,
// and the vertices corresponding to the coordinates of the points making
// up the triangles.
struct Electrode {
    std::string name;
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
};


int main(int argc, char **argv) {



    std::vector<Electrode> electrodes;

    Polygon_2 roi;
    roi.push_back(Point_2(80 * 100 / 9., 140 * 100 / 9.));
    roi.push_back(Point_2(80 * 100 / 9., 75 * 100 / 9.));
    roi.push_back(Point_2(160 * 100 / 9., 75 * 100 / 9.));
    roi.push_back(Point_2(160 * 100 / 9., 140 * 100 / 9.));


    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    std::string name;

    while (input >> name) {
        electrodes.push_back(Electrode());
        Electrode &electrode = electrodes.back();
        electrode.name = name;

        size_t nelectrodes;
        double z = 0.0;
        input >> nelectrodes;
        std::cout << "Processing " << name << std::endl;

        for (size_t e = 0; e < nelectrodes; ++e) {
            std::cout << " Electrode: " << e << std::endl;
            size_t npoints;
            input >> npoints;

            Polygon_2 whole_polygon;
            for (size_t p = 0; p < npoints - 1; ++p) {
                double x, y;
                input >> x >> y >> z;
                whole_polygon.push_back(Point_2(x, y));
            }
            double repeat_x, repeat_y;
            input >> repeat_x >> repeat_y >> z;

            Pwh_list_2 intersections;
            CGAL::intersection(roi, whole_polygon,
                               std::back_inserter(intersections));
            if (!intersections.size()) {
                continue;
            }

            std::cout << "Overlap found on " << name << ":" << e << std::endl;
            std::cout << intersections.size() << " polygons" << std::endl;

            for (Pwh_list_2::iterator pwh = intersections.begin();
                 pwh != intersections.end(); ++pwh) {

                const Polygon_2 &polygon = pwh->outer_boundary();
                std::cout << polygon.size() << " points" << std::endl;

                Converter converter;

                CDTriangulation cdt;
                //ConstraintList constraints;
                //std::transform( polygon.edges_begin(), polygon.edges_end(),
                //std::back_inserter( constraints  ), converter );
                //cdt.insert_constraints(
                //constraints.begin(), constraints.end() );

                CDTriangulation::Vertex_handle begin, prev;
                Polygon_2::Vertex_iterator v = polygon.vertices_begin();
                begin = prev = cdt.insert(converter(*v));
                for (++v; v != polygon.vertices_end(); ++v) {
                    CDTriangulation::Vertex_handle cur = cdt.insert(converter(*v));
                    cdt.insert_constraint(prev, cur);
                    prev = cur;
                }
                cdt.insert_constraint(prev, begin);

                double shape_criteria = 0.125;
                double size_criteria = 15.0;
                CGAL::refine_Delaunay_mesh_2(cdt,
                                             MeshSizeCriteria(shape_criteria, size_criteria));

                size_t mesh_nvert = cdt.number_of_vertices();
                std::cout << mesh_nvert << " vertices" << std::endl;

                std::map<CDTriangulation::Vertex_handle, size_t> indices;
                size_t index = 0;

                for (CDTriangulation::Vertex_iterator v = cdt.vertices_begin();
                     v != cdt.vertices_end(); ++v, ++index) {

                    Vertex vert;
                    vert.a = v->point().x();
                    vert.b = v->point().y();
                    vert.c = z;
                    electrode.vertices.push_back(vert);
                    indices[v] = electrode.vertices.size() - 1;
                }

                for (CDTriangulation::Finite_faces_iterator f = cdt.finite_faces_begin();
                     f != cdt.finite_faces_end(); ++f) {
                    if (!f->is_in_domain()) continue;

                    if (!indices.count(f->vertex(0)) ||
                        !indices.count(f->vertex(1)) ||
                        !indices.count(f->vertex(2))
                            ) {
                        std::cout << "BAD TRI: Bad " << name << ", "
                        << indices[f->vertex(0)] << ":"
                        << indices[f->vertex(1)] << ":"
                        << indices[f->vertex(2)] << std::endl;
                        continue;
                    }

                    Triangle tri;
                    tri.a = indices[f->vertex(0)];
                    tri.b = indices[f->vertex(1)];
                    tri.c = indices[f->vertex(2)];
                    electrode.triangles.push_back(tri);
                }
            }
        }

        if (!electrode.vertices.size() || !electrode.triangles.size())
            continue;

        output << "ELEC " << electrode.name << std::endl;
        output << "VERTS " << electrode.vertices.size() << std::endl;
        for (size_t v = 0; v < electrode.vertices.size(); ++v) {
            Vertex &vert = electrode.vertices[v];
            output << vert.a << " " << vert.b << " " << vert.c << std::endl;
        }

        output << "FACES " << electrode.triangles.size() << std::endl;
        for (size_t t = 0; t < electrode.triangles.size(); ++t) {
            Triangle &tri = electrode.triangles[t];
            output << tri.a << " " << tri.b << " " << tri.c << std::endl;
        }
    }


    output.close();
    return 0;
}
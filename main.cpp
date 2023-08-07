/*
This file is part of ``offset_surface'', a library for point based surface mesh offset.
Copyright (C) 2019 Bill He <github.com/easterngarden>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
// Author : Bill He
//
#include "C2t3_type.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <iostream>
#include <map>
#include <assert.h> 
#include "kdprocessor.h"

namespace CGAL {

	template <class TriangleMesh, class GeomTraits>
	class Offset_function
	{
		typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
		typedef AABB_traits<GeomTraits, Primitive> Traits;
		typedef AABB_tree<Traits> Tree;
		typedef Side_of_triangle_mesh<TriangleMesh, GeomTraits> Side_of;

	public:

		Offset_function(TriangleMesh& tm, double offset_distance, kdprocessor* &pProc)
			: m_tree_ptr(new Tree(boost::begin(faces(tm)),
				boost::end(faces(tm)),
				tm))
			, m_side_of_ptr(new Side_of(*m_tree_ptr))
			, m_offset_distance(offset_distance)
			, m_is_closed(is_closed(tm))
			, proc(pProc)
		{
			m_tree_ptr->accelerate_distance_queries();
		}

		double operator()(const typename GeomTraits::Point_3& p) const
		{
			using CGAL::sqrt;
			double factor = m_offset_distance > 0.0 ? 1.0 : -1.0;
			double offset_distance = proc ? factor * proc->findOffset(p.x(), p.y(), p.z()) : m_offset_distance;
			Bounded_side side = m_is_closed ? m_side_of_ptr->operator()(p) : ON_UNBOUNDED_SIDE;
			if (side == ON_BOUNDARY) return offset_distance;

			typename GeomTraits::Point_3 closest_point = m_tree_ptr->closest_point(p);
			double distance = sqrt(squared_distance(p, closest_point));

			return (side == ON_UNBOUNDED_SIDE ? -distance : distance) + offset_distance;
		}

	private:
		boost::shared_ptr<Tree> m_tree_ptr;
		boost::shared_ptr<Side_of> m_side_of_ptr;
		double m_offset_distance;
		bool m_is_closed;
		kdprocessor* proc;
	};

} //end of CGAL namespace

// declare the CGAL function
template<class Mesh>
int cgal_off_meshing(
	Mesh* tm_ptr,
	kdprocessor* pProc,
	const double offset_value,
	int tag)
{
	typedef Tr::Geom_traits GT;
	typedef CGAL::Offset_function<Mesh, GT> Offset_function;
	typedef CGAL::Implicit_surface_3<GT, Offset_function> Surface_3;
	typedef GT::Sphere_3 Sphere_3;

	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(*tm_ptr);

	GT::Point_3 center((bbox.xmax() + bbox.xmin()) / 2,
		(bbox.ymax() + bbox.ymin()) / 2,
		(bbox.zmax() + bbox.zmin()) / 2);
	double sqrad = 0.6 * std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
		CGAL::square(bbox.ymax() - bbox.ymin()) +
		CGAL::square(bbox.zmax() - bbox.zmin()))
		+ offset_value;
	sqrad = CGAL::square(sqrad);

	double dx = bbox.xmax() - bbox.xmin();
	double dy = bbox.ymax() - bbox.ymin();
	double dz = bbox.zmax() - bbox.zmin();
	double diag = std::sqrt(dx*dx + dy * dy + dz * dz);

	const double angle = 25.0;
	const double sizing = diag * 0.05;  //? default value
	const double approx = 0.1 * sizing;

	CGAL::Timer timer;
	timer.start();

	Offset_function offset_function(*tm_ptr, offset_value, pProc);

	Tr tr;
	C2t3 c2t3(tr);

	// defining the surface
	Surface_3 surface(offset_function,
		Sphere_3(center, sqrad)); // bounding sphere

    // defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle, sizing, approx);

	// meshing surface
	switch (tag) {
	case 0:
		CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
		break;
	case 1:
		CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
		break;
	default:
		CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

	}

	timer.stop();
	std::cerr << "done (" << timer.time() << " ms, " << c2t3.triangulation().number_of_vertices() << " vertices)" << std::endl;

	if (c2t3.triangulation().number_of_vertices() > 0)
	{
		// add remesh as new polyhedron
		Mesh *pRemesh = new Mesh;
		CGAL::facets_in_complex_2_to_triangle_mesh<C2t3, Mesh>(c2t3, *pRemesh);
		if (c2t3.number_of_facets() != num_faces(*pRemesh))
		{
			delete pRemesh;
			std::stringstream temp_file;
			if (!CGAL::output_surface_facets_to_off(temp_file, c2t3))
			{
				std::cerr << "Cannot write the mesh to an off file!\n";
				return 0;
			}
			std::ofstream outFile("offset.off");
			outFile << temp_file.rdbuf();
			outFile.close();
			return 1;
		}
		else {
			//pRemesh --> off
			std::ofstream outFile("offset.off");
			CGAL::write_off (outFile, *pRemesh);
			delete pRemesh;
			return 0;
		}
	}
	else
		return 0;
}

int main(int argc, char* argv[])
{
	assert(argc > 2);

	// create and read Polyhedron
	Polyhedron mesh;
	std::ifstream input(argv[1]);
	if (!input || !(input >> mesh) || mesh.empty() || (!CGAL::is_triangle_mesh(mesh))) {
		std::cerr << "Input is not a triangle mesh." << std::endl;
		return EXIT_FAILURE;
	}

	std::string inFile = argv[2];
	assert(!inFile.empty());
	kdprocessor* pProc = new kdprocessor();
	assert(pProc);
	double offset_value = -1.0; //? bounding offset
	pProc->buildtree(inFile, offset_value);

	CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
	const int tag = (argc > 3)? atoi(argv[3]) : 1; //CGAL::Manifold_tag()

	//offset surfurce mesh
	cgal_off_meshing(&mesh, pProc, offset_value, tag);
	
	if (pProc)
		delete pProc;
	return 0;
}


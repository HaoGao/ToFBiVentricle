// The Next Great Finite Element Library.
// Copyright (C) 2003  Benjamin S. Kirk
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Mesh library includes
//#include "libmesh_common.h"
#include "libmesh/libmesh_common.h"


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/point.h"
#include "libmesh/boundary_info.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;





/**
 * This is the exact solution that
 * we are trying to obtain.  We will solve
 *
 * - (u_xx + u_yy) = f
 *
 * and take a finite difference approximation using this
 * function to get f.  This is the well-known "method of
 * manufactured solutions".
 */
Real exact_solution (const Real x,
		     const Real y,
		     const Real z = 0.)
{
  static const Real pi = acos(-1.);

  return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
}


/** trying to assign u to gradux, a test whether it work
**/
void assign_gradient(EquationSystems& es,
                     const std::string& system_name)
{
	
   std::cout << "working on subsystem " << system_name << "\n"; 	

  //mesh is shared across the subsystems	
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  //this is the solution we are going to estimate the gradient
  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem& system_u = es.get_system<LinearImplicitSystem> ("Poisson");
  int u_system_num = system_u.number();
  
  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap& dof_map_u = system_u.get_dof_map();
  std::vector<unsigned int> u_dof_indices;
  system_u.solution->localize(*system_u.current_local_solution);
  NumericVector<double>& u_data = *(system_u.current_local_solution);
  u_data.close();
  
  //extract the gradient sub_system
  LinearImplicitSystem& system_GradUx = es.get_system<LinearImplicitSystem> (system_name);
  int GradUx_system_num = system_GradUx.number();
  const DofMap& dof_map_GradUx = system_GradUx.get_dof_map();
  std::vector<unsigned int> GradUx_dof_indices;
  NumericVector<double>& GradUx_data = *(system_GradUx.solution);
  
  //get constant reference to the finite element type
  FEType fe_type = dof_map_GradUx.variable_type(0);
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  
  
  
        MeshBase::const_element_iterator el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
  //loop all elements, using ++el
  for (; el != end_el;  ++el)
  {
	  const Elem* elem = *el;
	  
	  dof_map_u.dof_indices (elem, u_dof_indices);
	  dof_map_GradUx.dof_indices (elem, GradUx_dof_indices); 
	  
	  //update each element
	  fe->reinit(elem);
	  
	  for (unsigned int i=0; i<u_dof_indices.size(); ++i)
	  {
		  GradUx_data.set(GradUx_dof_indices[i], u_data(u_dof_indices[i]) );
	  }
	  
	  
  }
  
  GradUx_data.close();
  
  return ;	                
}

void assign_gradient_centre(EquationSystems& es,
                     const std::string& system_name)
{
	std::cout << "working on subsystem " << system_name << "\n"; 
	
	//mesh is shared across the subsystems	
    // Get a constant reference to the mesh object.
	const MeshBase& mesh = es.get_mesh();
	// The dimension that we are running
	const unsigned int dim = mesh.mesh_dimension();
  
	//this is the solution we are going to estimate the gradient
	// Get a reference to the LinearImplicitSystem we are solving
	LinearImplicitSystem& system_u = es.get_system<LinearImplicitSystem> ("Poisson");
	int u_system_num = system_u.number();
  
	// A reference to the  DofMap object for this system.  The  DofMap
	// object handles the index translation from node and element numbers
	// to degree of freedom numbers.  We will talk more about the  DofMap
	// in future examples.
	const DofMap& dof_map_u = system_u.get_dof_map();
	std::vector<unsigned int> u_dof_indices;
	system_u.solution->localize(*system_u.current_local_solution);
	NumericVector<double>& u_data = *(system_u.current_local_solution);
	u_data.close();
	
	
	//extract the gradient sub_system
	System& system_GradU = es.get_system<System> (system_name);
	int GradU_system_num = system_GradU.number();
	const DofMap& dof_map_GradU = system_GradU.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > GradU_dof_indices(dim);  
	NumericVector<double>& GradU_data = *(system_GradU.solution);
	
	
	
	//set up the FE type
	AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, CONSTANT);
	AutoPtr<FEBase> u_fe(FEBase::build(dim, dof_map_u.variable_type(0)));
	u_fe->attach_quadrature_rule(qrule.get());
	const std::vector<std::vector<RealGradient> >& dphi_u = u_fe->get_dphi();
	
	
	
	
	
        MeshBase::const_element_iterator el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
	//loop all elements, using ++el
	for (; el != end_el;  ++el)
	{
	  const Elem* elem = *el;
	  
	  dof_map_u.dof_indices (elem, u_dof_indices);
	  for (unsigned int d=0; d<dim; ++d)
	  {
			dof_map_GradU.dof_indices (elem, GradU_dof_indices[d], d); 
	  }
	  
	  //update each element
	  u_fe->reinit(elem);
	  
	  const unsigned int n_qp = qrule->n_points();
	  libmesh_assert(1 == n_qp);
	  const unsigned int qp = 0;
	  
	  std::vector<double> u_node (u_dof_indices.size() );
	  for (unsigned int node_i = 0; node_i < u_dof_indices.size(); ++node_i)
	  {
		  u_node[node_i] = u_data(u_dof_indices[node_i] );
	  }
	  //now calculate the gradient in the element centre
	  double GradUx, GradUy, GradUz;
	  GradUx = 0.0; 
	  GradUy = 0.0;
	  GradUz = 0.0;
	  for (unsigned int node_i = 0; node_i < u_dof_indices.size(); ++node_i)
	  {
		  const VectorValue<double>& dphi_ds = dphi_u[node_i][qp];
		  GradUx = GradUx + dphi_ds(0)*u_node[node_i];
		  GradUy = GradUy + dphi_ds(1)*u_node[node_i];
		  GradUz = GradUz + dphi_ds(2)*u_node[node_i];
	  }
	  
	  VectorValue<double> rad_dir (GradUx, GradUy, GradUz);
	  VectorValue<double> rad_dir_norm = rad_dir.unit();
	  
	  //save to the global solution vector
	  GradU_data.set( GradU_dof_indices[0][0], rad_dir_norm(0) );
	  GradU_data.set( GradU_dof_indices[1][0], rad_dir_norm(1) );
	  GradU_data.set( GradU_dof_indices[2][0], rad_dir_norm(2) );
	}
  
	GradU_data.close();
	
	
	return ;
} 




void assign_circum_centre(EquationSystems& es, 
                             const std::string system_name,
                             const Point global_centre)
{
	std::cout << "working on subsystem  "<< system_name <<"\n";
	//mesh is shared across the subsystems	
    // Get a constant reference to the mesh object.
	const MeshBase& mesh = es.get_mesh();
	// The dimension that we are running
	const unsigned int dim = mesh.mesh_dimension();

	//this is the solution we are going to estimate the gradient
	// Get a reference to the LinearImplicitSystem we are solving
	LinearImplicitSystem& system_u = es.get_system<LinearImplicitSystem> ("Poisson");
	int u_system_num = system_u.number();
  
	// A reference to the  DofMap object for this system.  The  DofMap
	// object handles the index translation from node and element numbers
	// to degree of freedom numbers.  We will talk more about the  DofMap
	// in future examples.
	const DofMap& dof_map_u = system_u.get_dof_map();
	std::vector<unsigned int> u_dof_indices;
	system_u.solution->localize(*system_u.current_local_solution);
	NumericVector<double>& u_data = *(system_u.current_local_solution);
	u_data.close();


    //extract the circum sub_system
	System& system_circum = es.get_system<System> (system_name);
	int circum_system_num = system_circum.number();
	const DofMap& dof_map_circum = system_circum.get_dof_map();
	//there are three components for the circum
	std::vector<std::vector<unsigned int> > circum_dof_indices(dim);  
	NumericVector<double>& circum_data = *(system_circum.solution);
	
	
	//extract the logi_sub_system
	System& system_logi = es.get_system<System> ("LogiCentre");
	int logi_system_num = system_logi.number();
	const DofMap& dof_map_logi = system_logi.get_dof_map();
	//there are three components for the circum
	std::vector<std::vector<unsigned int> > logi_dof_indices(dim);  
	NumericVector<double>& logi_data = *(system_logi.solution);
	
	
	//extract the rad_sub_system
	System& system_GradU = es.get_system<System> ("GradUCentre");
	int GradU_system_num = system_GradU.number();
	const DofMap& dof_map_GradU = system_GradU.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > GradU_dof_indices(dim);  
	NumericVector<double>& GradU_data = *(system_GradU.solution);
	
	
	
	//set up the FE type
	AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, CONSTANT);
	AutoPtr<FEBase> u_fe(FEBase::build(dim, dof_map_u.variable_type(0)));
	u_fe->attach_quadrature_rule(qrule.get());
	const std::vector<std::vector<RealGradient> >& dphi_u = u_fe->get_dphi();
	const std::vector<std::vector<Real> >& phi_u = u_fe->get_phi();
	const std::vector<Point>& q_point = u_fe->get_xyz();
	
	MeshBase::const_element_iterator el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
	//loop all elements, using ++el
	for (; el != end_el;  ++el)
	{
	  const Elem* elem = *el;
	  
	  dof_map_u.dof_indices (elem, u_dof_indices);
	  for (unsigned int d=0; d<dim; ++d)
	  {
			dof_map_circum.dof_indices (elem, circum_dof_indices[d], d);
			dof_map_logi.dof_indices   (elem,   logi_dof_indices[d], d);
			dof_map_GradU.dof_indices  (elem,  GradU_dof_indices[d], d); 
	  }
	  
	  //update each element
	  u_fe->reinit(elem);
	  
	  const unsigned int n_qp = qrule->n_points();
	  libmesh_assert(1 == n_qp);
	  const unsigned int qp = 0;
	  
	  
	  //now calculate the gradient in the element centre
	  double circum_x, circum_y;
	  Point p_centre = q_point[qp];
	  
	  //assume the centre location is [0,0]
	  //which is passed from global_centre
	  double xvec = p_centre(0) - global_centre(0);
	  double yvec = p_centre(1) - global_centre(1);
	  double zvec = 0.0;
	  double mag_vec = sqrt( xvec*xvec + yvec*yvec );
	  if (mag_vec<0.5)
	  {//the signual point in the apical region 
		  circum_x = 1.0;
		  circum_y = 0.0;	  
	  }
	  circum_x = yvec/mag_vec ;
	  circum_y = -xvec/mag_vec;	  
	  
	  
	  VectorValue<double> circum_dir(circum_x,circum_y,0);
	  VectorValue<double> circum_dir_norm = circum_dir.unit();
	  
	  double rad_x = GradU_data(  GradU_dof_indices[0][0]);
	  double rad_y = GradU_data(  GradU_dof_indices[1][0]);
	  double rad_z = GradU_data(  GradU_dof_indices[2][0]);
	  VectorValue<double> rad_dir(rad_x, rad_y, rad_z);
	  VectorValue<double> rad_dir_norm = rad_dir.unit();
	 
	  //here we will trust the radial direction is fixed, should not be updated
	  VectorValue<double> logi_dir = circum_dir_norm.cross(rad_dir_norm);
	  VectorValue<double> logi_dir_norm = logi_dir.unit();
	  
	  
	  //here we can add some further check if logi is not point the z positive
	  if (logi_dir_norm(3)<0)
	  {
		  libmesh_error();
	  }
	  
	  
	  //update the circumferential direction 
	  circum_dir = rad_dir_norm.cross(logi_dir_norm);
	  circum_dir_norm = circum_dir.unit();
	  
	  //save to the global solution vector for cir direction
	  circum_data.set( circum_dof_indices[0][0], circum_dir_norm(0));
	  circum_data.set( circum_dof_indices[1][0], circum_dir_norm(1));
	  circum_data.set( circum_dof_indices[2][0], circum_dir_norm(2));
	  
	  
	  //save to the global solution vector for logi_direction
	  logi_data.set( logi_dof_indices[0][0], logi_dir_norm(0));
	  logi_data.set( logi_dof_indices[1][0], logi_dir_norm(1));
	  logi_data.set( logi_dof_indices[2][0], logi_dir_norm(2));
	}
  
	circum_data.close();
	logi_data.close();
	
	
	return ;



}




















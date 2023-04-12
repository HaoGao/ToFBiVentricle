/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */


 // <h1>Introduction Example 3 - Solving a Poisson Problem</h1>
 //
 // This is the third example program.  It builds on
 // the second example program by showing how to solve a simple
 // Poisson system.  This example also introduces the notion
 // of customized matrix assembly functions, working with an
 // exact solution, and using element iterators.
 // We will not comment on things that
 // were already explained in the second example.

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
#include "libmesh/cell_tet4.h"

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

// read the options from a input file
#include "libmesh/getpot.h"

#include "tetVolumeCalculation.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;
using namespace std;

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems& es,
                      const std::string& system_name);

void assign_gradient(EquationSystems& es,
                     const std::string& system_name);
                     
void assign_gradient_centre(EquationSystems& es,
                     const std::string& system_name);  
                     
void assign_circum_centre(EquationSystems& es, 
                             const std::string system_name,
                             const Point global_centre);  
                             
void assign_fsn_centre(EquationSystems& es,
                       const std::string& system_name,
                       const std::string& circum_system_name,
                       const std::string& rad_system_name,
                       const std::string& logi_system_name, 
                       const unsigned int Num_Var,
                       const double alpha_fibre_endo, 
                       const double alpha_fibre_epi,
                       const double alpha_sheet_endo,
                       const double alpha_sheet_epi);  
                       
                       
void assign_fsn_centre(EquationSystems& es,
                       const std::string& system_name,
                       const std::string& circum_system_name,
                       const std::string& rad_system_name,
                       const std::string& logi_system_name,
                       const std::string& AHA_system_name, 
                       const unsigned int Num_Var,
                       std::map<int,  VectorValue< double> >& fibre_angle_aha);   
                       
                     
void output_fsn(EquationSystems& es,
					   const std::string& system_name,
                       const std::string& fibreDirName,
                       const std::string& sheetDirName,
                       const unsigned int Num_Var);
                     
void volumeCalculation(EquationSystems& es, double& volTotal);
              

// Function prototype for the exact solution.
Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.);

int main (int argc, char** argv)
{
  // Initialize libraries, like in example 2.
  LibMeshInit init (argc, argv);
  
  
  std::string node_file_name;
  std::string elem_file_name;
  std::string exodus_output_name;
  std::string xda_output_name;
  std::string result_output_name;
  std::string fibreDir_output_name;
  std::string sheetDir_output_name;
  
  std::string AHA_segments_file_name;
  std::string AHA_degree_file_name;
  
  double fibre_angle = 0.0;
  double sheet_angle = 0.0; 
  double using_AHASegments = 0.0;
  
  double alpha_fibre_endo = 0.0;
  double alpha_fibre_epi = 0.0;
  double alpha_sheet_endo = 0.0;
  double alpha_sheet_epi = 0.0;
  
  
  
  std::map<int, VectorValue<double> > fibre_angle_aha;
  
  
  
  std::string input_file_name = "mesh.in";
  if (argc >= 2)
  {
     input_file_name = argv[1];
  }
  else
  {
    std::cerr << "Missing input file, usage: exhg mesh.in\n";
  }
  GetPot input_file(input_file_name);
  
  if (input_file.have_variable("node_file_name") )
  {
     node_file_name = input_file("node_file_name","");
     std::cout << "node file name is "<< node_file_name << std::endl;
  }
  else std::cerr<< "Missing varialbe node_file_name, check your input file.\n";
  
  
  if (input_file.have_variable("elem_file_name") )
  {
     elem_file_name = input_file("elem_file_name","");
     std::cout << "elem file name is "<< elem_file_name << std::endl;
  }
  else std::cerr<< "Missing varialbe elem_file_name, check your input file.\n";
  
  if (input_file.have_variable("exodus_output_name") )
  {
     exodus_output_name = input_file("exodus_output_name","LVtrunk_heart_exodusII.e");
  }
  else std::cerr<< "Missing varialbe exodus_output_name, check your input file.\n";
  
  
   if (input_file.have_variable("xda_output_name") )
  {
     xda_output_name = input_file("xda_output_name","LVtrunk_heart_real.xda");
  }
  else std::cerr<< "Missing varialbe xda_output_name, check your input file.\n";
  
  if (input_file.have_variable("result_output_name") )
  {
     result_output_name = input_file("result_output_name","result.e");
  }
  else std::cerr<< "Missing varialbe result_output_name, check your input file.\n";
  
  if (input_file.have_variable("fibreDir_output_name") )
  {
     fibreDir_output_name = input_file("fibreDir_output_name","fibreDir.txt");
  }
  else std::cerr<< "Missing varialbe fibreDir_output_name, check your input file.\n";
  
  if (input_file.have_variable("sheetDir_output_name") )
  {
     sheetDir_output_name = input_file("sheetDir_output_name","sheetDir.txt");
  }
  else std::cerr<< "Missing varialbe sheetDir_output_name, check your input file.\n";
  
  if (input_file.have_variable("using_AHASegments"))
  {
	  using_AHASegments = input_file("using_AHASegments", 0.0);
	  std::cout << "the mesh files will be divided into 17 AHA segments and each segment can have different angle\n";
	  
	  if (input_file.have_variable("AHA_segments_file_name") )
	  {
		  AHA_segments_file_name = input_file("AHA_segments_file_name", "AHASeg.dat");
	  }
	  else std::cerr << "Missing variable AHA_segments_file_name, check your input file \n";
	  
	  
	  if (input_file.have_variable("AHA_degree_file_name") )
	  {
		  AHA_degree_file_name = input_file("AHA_degree_file_name", "AHADegree.dat");
	  }
	  else std::cerr << "Missing variable AHA_degree_file_name, check your input file \n";
	  
  }
  else   
  {
	  if (input_file.have_variable("fibre_angle")  )
	  {
		 fibre_angle = input_file("fibre_angle", 0.0);
		 std::cout << "reading fibre angle: "<<fibre_angle<<" in degree \n";
	  } 
	  else std::cerr <<"missing fibre angle, check your input file.\n";

	  if (input_file.have_variable("sheet_angle")  )
	  {
		 sheet_angle = input_file("sheet_angle", 0.0);
		 std::cout << "reading sheet angle: "<<sheet_angle<<" in degree \n";
	  } 
	  else std::cerr <<"missing sheet angle, check your input file.\n";
	  
	  
	  alpha_fibre_endo = -fibre_angle/180.0*M_PI;  //-M_PI/3;
	  alpha_fibre_epi  =  fibre_angle/180.0*M_PI;  //M_PI/3;
		
	  alpha_sheet_endo = -sheet_angle/180.0*M_PI;  //-M_PI/4;
	  alpha_sheet_epi  = sheet_angle/180.0*M_PI;   //M_PI/4;
  }  
  
  
  
  
  // Brief message to the user regarding the program name
  // and command line arguments.
  std::cout << "Running " << argv[0];
  
  for (int i=1; i<argc; i++)
    std::cout << " " << argv[i];
  
  std::cout << std::endl << std::endl;
  
  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");
  
  
  
  
  
    unsigned int dim = 3;
    ///reading tetgen mesh manually
	Mesh meshTet(init.comm(), dim);
	meshTet.allow_renumbering(false);
	
	ifstream ifsNode(node_file_name);
	//ifstream ifsEle("heart_real.1.ele");
   	ifstream ifsEle(elem_file_name);
	int NoOfEle, NoOfNode, tempv1, tempv2, tempv3;
	int IDNo, Node1, Node2, Node3, Node4;
	int side1B,side2B,side3B,side4B;
	int side1ID,side2ID,side3ID,side4ID;
	double CoorX, CoorY, CoorZ;
	int domainID;

    ifsNode>>NoOfNode>>tempv1>>tempv2>>tempv3;
	ifsEle>>NoOfEle>>tempv1>>tempv2;

	meshTet.reserve_nodes(NoOfNode);
	meshTet.reserve_elem(NoOfEle);

	///read node inside mesh
    int IDtemp=1;
	while (!ifsNode.eof()& IDtemp<=NoOfNode){
		ifsNode>>IDNo>>CoorX>>CoorY>>CoorZ;
		IDtemp++;
		Point p(CoorX,CoorY,CoorZ);
		meshTet.add_point(p,IDNo-1);
	}
	cout<<IDNo<<'\t'<<CoorX<<'\t'<<CoorY<<'\t'<<CoorZ<<endl;

	///read element inside mesh
	IDtemp = 1;
	while(!ifsEle.eof() & IDtemp<=NoOfEle){
		ifsEle>>IDNo>>Node1>>Node2>>Node3>>Node4>>side1B>>side2B>>side3B>>side4B>>side1ID>>side2ID>>side3ID>>side4ID>>domainID;
		Elem* elem = meshTet.add_elem(new Tet4);
		elem->set_node(0) = meshTet.node_ptr(Node1-1);
		elem->set_node(1) = meshTet.node_ptr(Node2-1);
		elem->set_node(2) = meshTet.node_ptr(Node3-1);
		elem->set_node(3) = meshTet.node_ptr(Node4-1);
		elem->set_id()=IDNo-1;
		elem->subdomain_id()=domainID;
     		
		//set up the boundary information 
		if(side1B == 1)
		{   	meshTet.boundary_info->add_side(elem,0,side1ID);
		}
		if(side2B == 1)
		{	meshTet.boundary_info->add_side(elem,1,side2ID);
		}
		if(side3B == 1)
		{	meshTet.boundary_info->add_side(elem,2,side3ID);
		}
		if(side4B == 1)
		{	meshTet.boundary_info->add_side(elem,3,side4ID);
		}
	
		IDtemp++;
        }
 	
	meshTet.prepare_for_use();
	meshTet.print_info();	
	
///output the mesh file
	#ifdef LIBMESH_HAVE_EXODUS_API
  		ExodusII_IO meshTet_writer(meshTet);
 		meshTet_writer.write(exodus_output_name);
  	//	//mesh_writer.write(argv[2]);
	#endif
//	TetGenIO meshTet_writer(meshTet);
//	meshTet.prepare_for_use();
//	meshTet_writer.write("heart_new");

  meshTet.write(xda_output_name);
//this is the end of reading mesh files	
	
  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());
  mesh.allow_renumbering(false);
  //Mesh mesh;
	
  mesh.read(xda_output_name, NULL);
  mesh.prepare_for_use();
  mesh.print_info();
  ///output the mesh file
  #ifdef LIBMESH_HAVE_EXODUS_API
  		ExodusII_IO mesh_writer(mesh);
 		mesh_writer.write("HeartGeometry.e");
  	//	//mesh_writer.write(argv[2]);
  #endif
  
  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  
  // Declare the Poisson system and its variables.
  // The Poisson system is another example of a steady system.
  equation_systems.add_system<LinearImplicitSystem> ("Poisson");
  
  
  equation_systems.add_system<LinearImplicitSystem>("GradUx");
  equation_systems.add_system<LinearImplicitSystem>("GradUy");
  equation_systems.add_system<LinearImplicitSystem>("GradUz");
  
  
  // Adds the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  //equation_systems.get_system("Poisson").add_variable("u", SECOND);
  equation_systems.get_system("Poisson").add_variable("u", FIRST);
  
  //now need to add the gradient system for calcualte ux, uy and uz in 
  //the nodes, not in the Gaussian Quadrature points 
  
  //equation_systems.get_system("GradUx").add_variable("grad_ux", FIRST);
  //equation_systems.get_system("GradUy").add_variable("grad_uy", FIRST);
  //equation_systems.get_system("GradUz").add_variable("grad_uz", FIRST);
  //equation_systems.get_system("GradientUy").add_variable("uy", FIRST);
  //equation_systems.get_system("GradientUz").add_variable("uz", FIRST);
 
 
  System& GradUCentre_system = equation_systems.add_system<System>("GradUCentre");
  GradUCentre_system.add_variable("grad_ux_centre", CONSTANT, MONOMIAL);
  GradUCentre_system.add_variable("grad_uy_centre", CONSTANT, MONOMIAL);
  GradUCentre_system.add_variable("grad_uz_centre", CONSTANT, MONOMIAL);
  
  System& CircumCentre_system = equation_systems.add_system<System>("CircumCentre");
  CircumCentre_system.add_variable("cir_ux_centre", CONSTANT, MONOMIAL);
  CircumCentre_system.add_variable("cir_uy_centre", CONSTANT, MONOMIAL);
  CircumCentre_system.add_variable("cir_uz_centre", CONSTANT, MONOMIAL);
  
  
  System& LogiCentre_system = equation_systems.add_system<System>("LogiCentre");
  LogiCentre_system.add_variable("log_ux_centre", CONSTANT, MONOMIAL);
  LogiCentre_system.add_variable("log_uy_centre", CONSTANT, MONOMIAL);
  LogiCentre_system.add_variable("log_uz_centre", CONSTANT, MONOMIAL);
  
  
  System& Myofibre_system = equation_systems.add_system<System>("Myofibre");
  Myofibre_system.add_variable("f0_x_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("f0_y_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("f0_z_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("s0_x_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("s0_y_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("s0_z_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("n0_x_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("n0_y_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("n0_z_centre", CONSTANT, MONOMIAL);
  Myofibre_system.add_variable("alpha",CONSTANT, MONOMIAL);
  const unsigned int Num_Var = 10;
  
  
  //add the AHA seg inforamtion into the system, element based
  System& AHA_system = equation_systems.add_system<System>("AHA");
  AHA_system.add_variable("AHA_ID", CONSTANT, MONOMIAL); 
  
  //add the angle system, element based
  System& angle_system = equation_systems.add_system<System>("Angle");
  angle_system.add_variable("fibre_angle", CONSTANT, MONOMIAL);
  angle_system.add_variable("sheet_angle", CONSTANT, MONOMIAL); 
  
  
	
  // Give the system a pointer to the matrix assembly
  // function.  This will be called when needed by the
  // library.
  
  equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
  
  // Initialize the data structures for the equation system.
  equation_systems.init();
  
  //read in the data file to assigan the AHAID
  if (using_AHASegments > 0.5)
  {
	  std::cout << "reading the file "<< AHA_segments_file_name<<"\n";
	  ifstream aha_stream(AHA_segments_file_name.c_str());
	  std::map<int, VectorValue<double> > aha_data;
	  for (unsigned int e=0; e < mesh.n_elem(); ++e)
	  {
		  int e_idx;
		  aha_stream >> e_idx;
		  VectorValue<double>& aha = aha_data[e];
		  aha_stream >> aha(0);
	  }
	  aha_stream.close();
	  
	  NumericVector<double> & aha_vec = *AHA_system.solution;
	  const int aha_system_num = AHA_system.number();
	  for ( MeshBase::element_iterator it = mesh.elements_begin(); it!=mesh.elements_end(); ++it)
	  {
		  Elem* elem = *it;
		  const int elem_id = elem->id();
		  const VectorValue<double> aha = aha_data[elem_id];
		  const int dof_index = elem->dof_number(aha_system_num, 0, 0);
		  aha_vec.set(dof_index, aha(0) );
	  } 
     aha_vec.close();
     aha_vec.localize(*AHA_system.current_local_solution); 
    
    //then we will need to read the degree for each into a standard map fibre_angle_aha
    std::cout <<"reading the degree from "<<AHA_degree_file_name<<"\n";
    ifstream aha_degree_stream(AHA_degree_file_name.c_str());
    int number_of_segments;
    aha_degree_stream >> number_of_segments;
    std::cout << "total segments is "<< number_of_segments<<"\n";
    for (unsigned int seg_id = 0; seg_id < number_of_segments; ++seg_id)
    {
		int aha_id;
		aha_degree_stream >> aha_id;
		VectorValue<double>& aha_degree = fibre_angle_aha[aha_id];
		aha_degree_stream >> aha_degree(0); 
		aha_degree_stream >> aha_degree(1);
		aha_degree_stream >> aha_degree(2);
		aha_degree_stream >> aha_degree(3);
	}
    
    aha_degree_stream.close();
    
    //the following code to check wether fibre_angle_aha is read properly
    std::cout << "fibre_angle_aha size is "<< fibre_angle_aha.size()<<"\n";
    for(std::map<int, VectorValue<double> >::iterator it=fibre_angle_aha.begin(); it!=fibre_angle_aha.end(); ++it)
    {
        int aha_id = it->first;
        VectorValue<double>& aha_degree = fibre_angle_aha[aha_id];
        std::cout << aha_id<<"\t"<<aha_degree(0)<<"\t"<<aha_degree(1)
                  <<"\t"<<aha_degree(2)<<"\t"<<aha_degree(3)<<"\n";
	}
    
    
  }
  
  
  
  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Solve the system "Poisson".  Note that calling this
  // member will assemble the linear system and invoke
  // the default numerical solver.  With PETSc the solver can be
  // controlled from the command line.  For example,
  // you can invoke conjugate gradient with:
  //
  // ./ex3 -ksp_type cg
  //
  // You can also get a nice X-window that monitors the solver
  // convergence with:
  //
  // ./ex3 -ksp_xmonitor
  //
  // if you linked against the appropriate X libraries when you
  // built PETSc.
  equation_systems.get_system("Poisson").solve();
  
  
  
  //solve other systems here 
  //assign_gradient(equation_systems, "GradUx");
  

  //extract gradient at centre
  assign_gradient_centre(equation_systems, "GradUCentre");
  
  //extract the circum at centre, and the longitidual direction
  Point global_centre(0,0,0);
  assign_circum_centre(equation_systems, "CircumCentre", global_centre);
  
 
  // assing the fsn direction 
  if (using_AHASegments > 0.5)
  {
	  assign_fsn_centre(equation_systems,
                        "Myofibre",
                       "CircumCentre",
                       "GradUCentre",
                       "LogiCentre",
                       "AHA", 
                       Num_Var,
                       fibre_angle_aha);
  }
  else
  {
	assign_fsn_centre(equation_systems,
                       "Myofibre",
                       "CircumCentre",
                       "GradUCentre",
                       "LogiCentre",
                       Num_Var,
                       alpha_fibre_endo, alpha_fibre_epi,
                       alpha_sheet_endo, alpha_sheet_epi);
 }
 
 
 
 //out put the fibre and sheet directions
  output_fsn(equation_systems,
						"Myofibre",
                       fibreDir_output_name,
                       sheetDir_output_name,
                       Num_Var);
  
/* comment out this since it will need to read a endo list 
 //calculate the LV cavity volume
 double volTotal = 0;
 if (0 == equation_systems.processor_id())
     {volumeCalculation(equation_systems, volTotal);}
  MPI_Bcast(&volTotal, 1, MPI_DOUBLE, 0, equation_systems.comm().get());
  
  MPI_Barrier(equation_systems.comm().get());
  printf("From processor %d, volTotal is %f\n", equation_systems.comm().rank(), volTotal);

//std::cout<<"from processor " <<init.comm().rank()<<"  ,the total volume is "<<volTotal<<"\n";
  */

#ifdef LIBMESH_HAVE_EXODUS_API
	cout<<"write solution"<<endl;
	ExodusII_IO (mesh).write_equation_systems(result_output_name, equation_systems);
#endif

//output fibreDir.dat and sheetDir.dat


  // All done.  
  return 0;
}



// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_poisson(EquationSystems& es,
                      const std::string& system_name)
{
  
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert (system_name == "Poisson");

  //std::cout << "assemble the system\n"; 
 
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem> ("Poisson");

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap& dof_map = system.get_dof_map();
  
  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);
  
  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as an AutoPtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.  Example 4
  // describes some advantages of  AutoPtr's in the context of
  // quadrature rules.
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);
  
  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  
  // Declare a special finite element object for
  // boundary integration.
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  
  // Boundary integration requires one quadraure rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);
  
  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point>& q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".  These datatypes are templated on
  //  Number, which allows the same code to work for real
  // or complex numbers.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;


  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // iterator el will iterate from the first to the last element on
  // the local processor.  The iterator end_el tells us when to stop.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
 
  // Loop over the elements.  Note that  ++el is preferred to
  // el++ since the latter requires an unnecessary temporary
  // object.
  for ( ; el != end_el ; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);


      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).

      // The  DenseMatrix::resize() and the  DenseVector::resize()
      // members will automatically zero out the matrix  and vector.
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          // Now we will build the element matrix.  This involves
          // a double loop to integrate the test funcions (i) against
          // the trial functions (j).
          for (unsigned int i=0; i<phi.size(); i++)
            for (unsigned int j=0; j<phi.size(); j++)
              {
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              }
          
          // This is the end of the matrix summation loop
          // Now we build the element right-hand-side contribution.
          // This involves a single loop in which we integrate the
          // "forcing function" in the PDE against the test functions.
          {
            const Real x = q_point[qp](0);
            const Real y = q_point[qp](1);
            const Real eps = 1.e-3;
            

            // "fxy" is the forcing function for the Poisson equation.
            // In this case we set fxy to be a finite difference
            // Laplacian approximation to the (known) exact solution.
            //
            // We will use the second-order accurate FD Laplacian
            // approximation, which in 2D is
            //
            // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
            //                u(i-1,j) + u(i+1,j) +
            //                -4*u(i,j))/h^2
            //
            // Since the value of the forcing function depends only
            // on the location of the quadrature point (q_point[qp])
            // we will compute it here, outside of the i-loop
            //const Real fxy = -(exact_solution(x,y-eps) +
            //                   exact_solution(x,y+eps) +
            //                   exact_solution(x-eps,y) +
            //                   exact_solution(x+eps,y) -
            //                   4.*exact_solution(x,y))/eps/eps;
            const Real fxy = 0.0;
			
            for (unsigned int i=0; i<phi.size(); i++)
              Fe(i) += JxW[qp]*fxy*phi[i][qp];
          } 
        } 
      
      // We have now reached the end of the RHS summation,
      // and the end of quadrature point loop, so
      // the interior element integration has
      // been completed.  However, we have not yet addressed
      // boundary conditions.  For this example we will only
      // consider simple Dirichlet boundary conditions.
      //
      // There are several ways Dirichlet boundary conditions
      // can be imposed.  A simple approach, which works for
      // interpolary bases like the standard Lagrange polynomials,
      // is to assign function values to the
      // degrees of freedom living on the domain boundary. This
      // works well for interpolary bases, but is more difficult
      // when non-interpolary (e.g Legendre or Hierarchic) bases
      // are used.
      //
      // Dirichlet boundary conditions can also be imposed with a
      // "penalty" method.  In this case essentially the L2 projection
      // of the boundary values are added to the matrix. The
      // projection is multiplied by some large factor so that, in
      // floating point arithmetic, the existing (smaller) entries
      // in the matrix and right-hand-side are effectively ignored.
      //
      // This amounts to adding a term of the form (in latex notation)
      //
      // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
      //
      // where
      //
      // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
      {

        // The following loop is over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int side=0; side<elem->n_sides(); side++)
        {	
			//const 
			vector<short int> bdry_ids; 
			mesh.boundary_info->boundary_ids(elem, side, bdry_ids);
			// The penalty value.  \frac{1}{\epsilon}
                  // in the discussion above.
            const Real penalty = 1.e10;
			
			short int bdry_id_temp = 0;
			if (bdry_ids.size() >=1)
			  bdry_id_temp = bdry_ids[0];
			
			if ( (elem->neighbor(side) == libmesh_nullptr)  &&  (bdry_id_temp == 4096  || bdry_id_temp == 4097 
			                                                  || bdry_id_temp == 5090  || bdry_id_temp ==5091) )
            {
              // The boundary value.
			  //const Real value = exact_solution(xf, yf);
			  Real value = 0.0;
			  if (bdry_id_temp == 4096) //epi
			  { value = 1.0;}
			  else if (bdry_id_temp == 4097) //endo_LV
			  { value = 0.0;}
			  else if (bdry_id_temp == 5091) //endo_RV_septum
			  { value = 1.0; }
			  else if (bdry_id_temp == 5090) //endo_RV_freewall
			  {	value = 0.0; }
			  else
			  { std::cout<<"in a wrong boundary\n"; }
				  
			
			  // The value of the shape functions at the quadrature
              // points.
              const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
              
              // The Jacobian * Quadrature Weight at the quadrature
              // points on the face.
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
              
              // The XYZ locations (in physical space) of the
              // quadrature points on the face.  This is where
              // we will interpolate the boundary value function.
              const std::vector<Point >& qface_point = fe_face->get_xyz();
              
              // Compute the shape function values on the element
              // face.
              fe_face->reinit(elem, side);
              
              // Loop over the face quadrature points for integration.
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {

                  // The location on the boundary of the current
                  // face quadrature point.
                  const Real xf = qface_point[qp](0);
                  const Real yf = qface_point[qp](1);

                  
				  
                  // Matrix contribution of the L2 projection. 
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

                  // Right-hand-side contribution of the L2
                  // projection.
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                } 
            }// for bc 4091
			
			
		} //for side bc application
      }
      
      // We have now finished the quadrature point loop,
      // and have therefore applied all the boundary conditions.

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The  SparseMatrix::add_matrix()
      // and  NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
  
  // All done!
}




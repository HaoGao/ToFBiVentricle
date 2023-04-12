// Mesh library includes
//#include "libmesh_common.h"
#include "libmesh/libmesh_common.h"


// C++ include files that we need
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
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

#include "tetVolumeCalculation.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


void assign_fsn_centre(EquationSystems& es,
                       const std::string& system_name,
                       const std::string& circum_system_name,
                       const std::string& rad_system_name,
                       const std::string& logi_system_name, 
                       const unsigned int Num_Var,
                       const double alpha_fibre_endo, 
                       const double alpha_fibre_epi,
                       const double alpha_sheet_endo,
                       const double alpha_sheet_epi)
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


    //extract the circum sub_system
	System& system_circum = es.get_system<System> (circum_system_name);
	int circum_system_num = system_circum.number();
	const DofMap& dof_map_circum = system_circum.get_dof_map();
	//there are three components for the circum
	std::vector<std::vector<unsigned int> > circum_dof_indices(dim);  
	NumericVector<double>& circum_data = *(system_circum.solution);
	
	
	//extract the logi_sub_system
	System& system_logi = es.get_system<System> (logi_system_name);
	int logi_system_num = system_logi.number();
	const DofMap& dof_map_logi = system_logi.get_dof_map();
	//there are three components for the circum
	std::vector<std::vector<unsigned int> > logi_dof_indices(dim);  
	NumericVector<double>& logi_data = *(system_logi.solution);
	
	
	//extract the rad_sub_system
	System& system_GradU = es.get_system<System> (rad_system_name);
	int GradU_system_num = system_GradU.number();
	const DofMap& dof_map_GradU = system_GradU.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > GradU_dof_indices(dim);  
	NumericVector<double>& GradU_data = *(system_GradU.solution);
	
	//extract the fsn_sub_system
	System& system_fsn = es.get_system<System> (system_name);
	int fsn_system_num = system_fsn.number();
	const DofMap& dof_map_fsn = system_fsn.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > fsn_dof_indices(Num_Var);  
	NumericVector<double>& fsn_data = *(system_fsn.solution);
	
	
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
		for (unsigned int d=0; d< Num_Var; ++d)
		{
			dof_map_fsn.dof_indices (elem, fsn_dof_indices[d],d);
		}
	  
		//update each element
		u_fe->reinit(elem);
		const unsigned int n_qp = qrule->n_points();
		libmesh_assert(1 == n_qp);
		const unsigned int qp = 0;
	  
		//going to assign the data
		std::vector<double> u_node (u_dof_indices.size() );
		for (unsigned int node_i = 0; node_i < u_dof_indices.size(); ++node_i)
		{
			u_node[node_i] = u_data(u_dof_indices[node_i] );
		}
		
		double u_centre = 0;
		for (unsigned int node_i = 0; node_i < u_dof_indices.size(); ++node_i)
		{
			u_centre = u_centre + u_node[node_i]*phi_u[node_i][qp];
		}
		
		
		double rad_x = GradU_data(  GradU_dof_indices[0][0]);
	    double rad_y = GradU_data(  GradU_dof_indices[1][0]);
		double rad_z = GradU_data(  GradU_dof_indices[2][0]);
		VectorValue<double> rad_dir(rad_x, rad_y, rad_z);
		VectorValue<double> rad_dir_norm = rad_dir.unit();
		
		double circum_x = circum_data(  circum_dof_indices[0][0]);
	    double circum_y = circum_data(  circum_dof_indices[1][0]);
		double circum_z = circum_data(  circum_dof_indices[2][0]);
		VectorValue<double> circum_dir(circum_x, circum_y, circum_z);
		VectorValue<double> circum_dir_norm = circum_dir.unit();
	  
	    double logi_x = logi_data(  logi_dof_indices[0][0]);
	    double logi_y = logi_data(  logi_dof_indices[1][0]);
		double logi_z = logi_data(  logi_dof_indices[2][0]);
		VectorValue<double> logi_dir(logi_x, logi_y, logi_z);
		VectorValue<double> logi_dir_norm = logi_dir.unit();
	  
	    
	    
	    
	    VectorValue<double> f0(0,0,0);
	    VectorValue<double> s0(0,0,0);
	    VectorValue<double> n0(0,0,0);
	    
	    //the alpha at each element center 
	    double u_endo = u_centre;
	    double u_epi   = 1 - u_centre;
	    double alpha_fibre_centre = alpha_fibre_endo + (alpha_fibre_epi 
	                                         - alpha_fibre_endo)*u_endo;
	    
	    
	    //fibre direction in the circum_dir_norm and logi_dir_norm plane 
	    f0 = circum_dir_norm*cos(alpha_fibre_centre) + 
	           logi_dir_norm*sin(alpha_fibre_centre);
	    VectorValue<double> f0_norm = f0.unit();
	    
	    //sheet direction in the rad_dir_norm and -logi_dir_norm
	    double alpha_sheet_centre = alpha_sheet_endo + (alpha_sheet_epi - alpha_sheet_endo)*u_endo;
	    s0 = rad_dir_norm * cos(alpha_sheet_centre) + logi_dir_norm*(sin(alpha_sheet_centre));
	    VectorValue<double> s0_norm = s0.unit();
	    
	    //logitudinal direction 
	    n0 = f0_norm.cross(s0_norm);
	    VectorValue<double> n0_norm = n0.unit();
	    
	    //readjust s0 direction
	    s0 = n0_norm.cross(f0_norm);
	    s0_norm = s0.unit();
	    
	    fsn_data.set(fsn_dof_indices[0][0],f0_norm(0) );
	    fsn_data.set(fsn_dof_indices[1][0],f0_norm(1) );
	    fsn_data.set(fsn_dof_indices[2][0],f0_norm(2) );
	    
	    fsn_data.set(fsn_dof_indices[3][0],s0_norm(0) );
	    fsn_data.set(fsn_dof_indices[4][0],s0_norm(1) );
	    fsn_data.set(fsn_dof_indices[5][0],s0_norm(2) );
	     
	    fsn_data.set(fsn_dof_indices[6][0],n0_norm(0) );
	    fsn_data.set(fsn_dof_indices[7][0],n0_norm(1) );
	    fsn_data.set(fsn_dof_indices[8][0],n0_norm(2) );
	    
	    fsn_data.set(fsn_dof_indices[9][0],alpha_fibre_centre );
	    
	    
	    
	  
     } //end for el
	
	fsn_data.close();
		
		
	return ;	
		   
}



void output_fsn(EquationSystems& es,
					   const std::string& system_name,
                       const std::string& fibreDirName,
                       const std::string& sheetDirName,
                       const unsigned int Num_Var)
{
	std::cout << "output fibre and sheet direction "<< "\n"; 
	
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
	
	//extract the fsn_sub_system
	System& system_fsn = es.get_system<System> (system_name);
	int fsn_system_num = system_fsn.number();
	const DofMap& dof_map_fsn = system_fsn.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > fsn_dof_indices(Num_Var); 
	system_fsn.solution->localize(*system_fsn.current_local_solution); 
	NumericVector<double>& fsn_data = *(system_fsn.current_local_solution);
	std::vector<double> fsn_data_vec;
	system_fsn.solution->localize_to_one(fsn_data_vec);
	
	
	if (0 == es.processor_id())
	{
		std::cout << "writing output from processor:  "<< es.processor_id()<<std::endl;
		
		//open the file to write
		std::fstream fs_fibre;
		fs_fibre.open(fibreDirName.c_str(), std::fstream::out);
		
		std::fstream fs_sheet;
		fs_sheet.open(sheetDirName.c_str(), std::fstream::out);
		
		
		//using the global el iterator
		MeshBase::const_element_iterator el     = mesh.active_elements_begin();
		const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
	
		//loop all elements, using ++el
		for (; el != end_el;  ++el)
		{
				const Elem* elem = *el;
				for (unsigned int d=0; d< Num_Var; ++d)
				{
					dof_map_fsn.dof_indices (elem, fsn_dof_indices[d],d);
				}
				
				//double f0_x = fsn_data(  fsn_dof_indices[0][0]);
	            //double f0_y = fsn_data(  fsn_dof_indices[1][0]);
	            //double f0_z = fsn_data(  fsn_dof_indices[2][0]);
	            
	            //double s0_x = fsn_data ( fsn_dof_indices[3][0] );
	            //double s0_y = fsn_data ( fsn_dof_indices[4][0] );
				//double s0_z = fsn_data ( fsn_dof_indices[5][0] );
				
				double f0_x = fsn_data_vec[  fsn_dof_indices[0][0] ];
	            double f0_y = fsn_data_vec[  fsn_dof_indices[1][0] ];
	            double f0_z = fsn_data_vec[  fsn_dof_indices[2][0] ];
	            
	            double s0_x = fsn_data_vec [ fsn_dof_indices[3][0] ];
	            double s0_y = fsn_data_vec [ fsn_dof_indices[4][0] ];
				double s0_z = fsn_data_vec [ fsn_dof_indices[5][0] ];
	  
				fs_fibre << elem->id()+1 <<"    "<< f0_x << "    " <<f0_y <<"    " <<f0_z<<"\n";
				fs_sheet << elem->id()+1 <<"    "<< s0_x << "    " <<s0_y <<"    " <<s0_z<<"\n";

				
		}
		
		fs_fibre.close();
		fs_sheet.close();
		
	}
	
	
	return;
}


void volumeCalculation(EquationSystems& es, double& volTotal)
{
	
	
	const MeshBase& mesh = es.get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	
	LinearImplicitSystem& system_u = es.get_system<LinearImplicitSystem> ("Poisson");
	int u_system_num = system_u.number();
	
	if (0 == es.processor_id())
	{
		
	
	   std::cout << "calculate the LV cavity volume from processor "<< es.processor_id()<<"\n"; 
	
	   std::vector< std::vector<double> > endo_points;
	   std::ifstream ifsendo("endoList.point");
	   int NoOfEndoNode, nodeEndoID;
	   
	   ifsendo>>NoOfEndoNode;
	   endo_points.resize(NoOfEndoNode);
	   for (unsigned int i = 0; i< NoOfEndoNode; i++)
	   {
		   endo_points[i].resize(3);
	   }
	   
	   unsigned int IDtemp=1; //reused from the initial defintion
	   unsigned int pIndex = 0;
	   Point  node_ref;
	   //initialize end_points
	   while (!ifsendo.eof()& IDtemp<=NoOfEndoNode){
			ifsendo>>nodeEndoID;
			IDtemp++;
			nodeEndoID = nodeEndoID - 1 ; //start from 0		
		
			node_ref = mesh.point(nodeEndoID);		
					
			endo_points[pIndex][0] = node_ref(0);
			endo_points[pIndex][1] = node_ref(1);
			endo_points[pIndex][2] = node_ref(2);
			
			pIndex = pIndex + 1;
			
		}
	   	
		 volTotal = tetVolumeCalculation(endo_points,  NoOfEndoNode);	
		 //std::cout<< "total volume from tet function:  "<<volTotal<<"\n";
		 
		 //send to other processors 
		 //MPI_Send(&volTotal, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		 
	}
	
	//MPI_Bcast(&volTotal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	/*For MPI collective communications, everyone has to particpate; 
	 * everyone has to call the Bcast, or the Allreduce, or what have you. (
	 * That's why the Bcast routine has a parameter that specifies the "root",
	 *  or who is doing the sending; if only the sender called bcast, you wouldn't 
	 * need this.) Everyone calls the broadcast, including the receivers; 
	 * the receviers don't just post a receive.*/
	//std::cout<<"from processor " <<es.processor_id()<<"  ,the total volume is "<<volTotal<<"\n";
	
	
	
	return;
}



void assign_fsn_centre(EquationSystems& es,
                       const std::string& system_name,
                       const std::string& circum_system_name,
                       const std::string& rad_system_name,
                       const std::string& logi_system_name,
                       const std::string& AHA_system_name, 
                       const unsigned int Num_Var,
                       std::map< int, VectorValue< double> >& fibre_angle_aha)
{
	std::cout << "working on subsystem " << system_name << "using AHA defined angle\n"; 
	
	
	
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
	System& system_circum = es.get_system<System> (circum_system_name);
	int circum_system_num = system_circum.number();
	const DofMap& dof_map_circum = system_circum.get_dof_map();
	//there are three components for the circum
	std::vector<std::vector<unsigned int> > circum_dof_indices(dim);  
	NumericVector<double>& circum_data = *(system_circum.solution);
	
	
	//extract the logi_sub_system
	System& system_logi = es.get_system<System> (logi_system_name);
	int logi_system_num = system_logi.number();
	const DofMap& dof_map_logi = system_logi.get_dof_map();
	//there are three components for the circum
	std::vector<std::vector<unsigned int> > logi_dof_indices(dim);  
	NumericVector<double>& logi_data = *(system_logi.solution);
	
	
	//extract the rad_sub_system
	System& system_GradU = es.get_system<System> (rad_system_name);
	int GradU_system_num = system_GradU.number();
	const DofMap& dof_map_GradU = system_GradU.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > GradU_dof_indices(dim);  
	NumericVector<double>& GradU_data = *(system_GradU.solution);
	
	//extract the fsn_sub_system
	System& system_fsn = es.get_system<System> (system_name);
	int fsn_system_num = system_fsn.number();
	const DofMap& dof_map_fsn = system_fsn.get_dof_map();
	//there are three components for the gradient
	std::vector<std::vector<unsigned int> > fsn_dof_indices(Num_Var);  
	NumericVector<double>& fsn_data = *(system_fsn.solution);
	
	
	System& system_aha = es.get_system<System> (AHA_system_name);
	int aha_system_num = system_aha.number();
	const DofMap& dof_map_aha = system_aha.get_dof_map();
	std::vector<unsigned int> aha_dof_indices;
	NumericVector<double>& aha_data = *(system_aha.solution);
	
	//extract the angle system and saved the angle
	System& system_angle = es.get_system<System> ("Angle");
	int angle_system_num = system_angle.number();
	const DofMap& dof_map_angle = system_angle.get_dof_map();
	std::vector<std::vector <unsigned int> > angle_dof_indices(2); //angle_fibre and angle_sheet
	NumericVector<double>& angle_data = *(system_angle.solution);
	 
	//set up the FE type
	AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, CONSTANT);
	AutoPtr<FEBase> u_fe(FEBase::build(dim, dof_map_u.variable_type(0)));
	u_fe->attach_quadrature_rule(qrule.get());
	const std::vector<std::vector<RealGradient> >& dphi_u = u_fe->get_dphi();
	const std::vector<std::vector<Real> >& phi_u = u_fe->get_phi();
	const std::vector<Point>& q_point = u_fe->get_xyz();
	
	
	for(std::map<int, VectorValue<double> >::iterator it=fibre_angle_aha.begin(); it!=fibre_angle_aha.end(); ++it)
    {
        int aha_id = it->first;
        VectorValue<double>& aha_degree = fibre_angle_aha[aha_id];
        std::cout << aha_id<<"\t"<<aha_degree(0)<<"\t"<<aha_degree(1)
                  <<"\t"<<aha_degree(2)<<"\t"<<aha_degree(3)<<"\n";
	}
	
	
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
		for (unsigned int d=0; d< Num_Var; ++d)
		{
			dof_map_fsn.dof_indices (elem, fsn_dof_indices[d],d);
		}
		
		for (unsigned int d=0; d<2; ++d)
		{
			dof_map_angle.dof_indices (elem, angle_dof_indices[d], d);
		}
		
		// figure out what the angle should be for this element
		dof_map_aha.dof_indices(elem, aha_dof_indices, 0);
		int aha_id_elem;
		aha_id_elem  = (int)(aha_data(aha_dof_indices[0]));
		VectorValue<double>& aha_degree  = fibre_angle_aha[aha_id_elem];
		double alpha_fibre_endo = aha_degree(0);
		double alpha_fibre_epi = aha_degree(1);
		double alpha_sheet_endo = aha_degree(2);
		double alpha_sheet_epi = aha_degree(3);
		//std::cout << aha_id_elem<<"\t"<<alpha_fibre_endo <<"\t"<<alpha_fibre_epi<<"\t"<<alpha_sheet_endo<<"\t"<<alpha_sheet_epi<<"\n";
		
	  
		//update each element
		u_fe->reinit(elem);
		const unsigned int n_qp = qrule->n_points();
		libmesh_assert(1 == n_qp);
		const unsigned int qp = 0;
	  
		//going to assign the data
		std::vector<double> u_node (u_dof_indices.size() );
		for (unsigned int node_i = 0; node_i < u_dof_indices.size(); ++node_i)
		{
			u_node[node_i] = u_data(u_dof_indices[node_i] );
		}
		
		double u_centre = 0;
		for (unsigned int node_i = 0; node_i < u_dof_indices.size(); ++node_i)
		{
			u_centre = u_centre + u_node[node_i]*phi_u[node_i][qp];
		}
		
		
		double rad_x = GradU_data(  GradU_dof_indices[0][0]);
	    double rad_y = GradU_data(  GradU_dof_indices[1][0]);
		double rad_z = GradU_data(  GradU_dof_indices[2][0]);
		VectorValue<double> rad_dir(rad_x, rad_y, rad_z);
		VectorValue<double> rad_dir_norm = rad_dir.unit();
		
		double circum_x = circum_data(  circum_dof_indices[0][0]);
	    double circum_y = circum_data(  circum_dof_indices[1][0]);
		double circum_z = circum_data(  circum_dof_indices[2][0]);
		VectorValue<double> circum_dir(circum_x, circum_y, circum_z);
		VectorValue<double> circum_dir_norm = circum_dir.unit();
	  
	    double logi_x = logi_data(  logi_dof_indices[0][0]);
	    double logi_y = logi_data(  logi_dof_indices[1][0]);
		double logi_z = logi_data(  logi_dof_indices[2][0]);
		VectorValue<double> logi_dir(logi_x, logi_y, logi_z);
		VectorValue<double> logi_dir_norm = logi_dir.unit();
	  
	    
	    
	    
	    VectorValue<double> f0(0,0,0);
	    VectorValue<double> s0(0,0,0);
	    VectorValue<double> n0(0,0,0);
	    
	    //the alpha at each element center 
	    double u_endo = u_centre;
	    double u_epi   = 1 - u_centre;
	    double alpha_fibre_centre = alpha_fibre_endo + (alpha_fibre_epi 
	                                         - alpha_fibre_endo)*u_endo;
	    angle_data.set(angle_dof_indices[0][0], alpha_fibre_centre);                                     
	    
	    
	    //fibre direction in the circum_dir_norm and logi_dir_norm plane 
	    f0 = circum_dir_norm*cos(alpha_fibre_centre) + 
	           logi_dir_norm*sin(alpha_fibre_centre);
	    VectorValue<double> f0_norm = f0.unit();
	    
	    //sheet direction in the rad_dir_norm and -logi_dir_norm
	    double alpha_sheet_centre = alpha_sheet_endo + (alpha_sheet_epi 
	                                         - alpha_sheet_endo)*u_endo;
	    angle_data.set(angle_dof_indices[1][0], alpha_sheet_centre);
	    //std::cout<< "alpha_sheet_centre is "<< alpha_sheet_centre<<"\t"<< alpha_sheet_epi<<"\n";
	    
	    s0 = rad_dir_norm * cos(alpha_sheet_centre) + logi_dir_norm*(sin(alpha_sheet_centre));
	    VectorValue<double> s0_norm = s0.unit();
	    
	    //logitudinal direction 
	    n0 = f0_norm.cross(s0_norm);
	    VectorValue<double> n0_norm = n0.unit();
	    
	    //readjust s0 direction
	    s0 = n0_norm.cross(f0_norm);
	    s0_norm = s0.unit();
	    
	    fsn_data.set(fsn_dof_indices[0][0],f0_norm(0) );
	    fsn_data.set(fsn_dof_indices[1][0],f0_norm(1) );
	    fsn_data.set(fsn_dof_indices[2][0],f0_norm(2) );
	    
	    fsn_data.set(fsn_dof_indices[3][0],s0_norm(0) );
	    fsn_data.set(fsn_dof_indices[4][0],s0_norm(1) );
	    fsn_data.set(fsn_dof_indices[5][0],s0_norm(2) );
	     
	    fsn_data.set(fsn_dof_indices[6][0],n0_norm(0) );
	    fsn_data.set(fsn_dof_indices[7][0],n0_norm(1) );
	    fsn_data.set(fsn_dof_indices[8][0],n0_norm(2) );
	    
	    fsn_data.set(fsn_dof_indices[9][0],alpha_fibre_centre );
	    
	    
	    
	  
     } //end for el
	
	fsn_data.close();
	angle_data.close();
		
		
	return ;
}

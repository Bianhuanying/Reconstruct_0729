/*================================================================
*   Copyright (C) 2020 BHY. All rights reserved.
*   
*   文件名称：Reconstruct.h
*   创 建 者：HuanyingBIAN
*   创建日期：2020年07月13日
*   描    述：
*
================================================================*/

#ifndef _RECONSTRUCT_H_
#define _RECONSTRUCT_H_

#include <sstream>
#include <fstream>
#include <iostream>
 
#include <AFEPack/EasyMesh.h> 
#include <AFEPack/Functional.h> 
#include <AFEPack/Operator.h> 
#include <AFEPack/FEMSpace.h>
#include <AFEPack/TemplateElement.h> 
#include <AFEPack/Geometry.h>
#include <AFEPack/AMGSolver.h>
 
#include <lac/sparse_mic.h> 
#include <lac/sparse_ilu.h> 
#include <lac/solver_minres.h> 
#include <lac/solver_gmres.h> 
#include <lac/solver_bicgstab.h> 
#include <lac/solver_cg.h> 
#include <lac/precondition.h> 
#include <lac/sparse_ilu.h> 
#include <lac/sparsity_pattern.h> 
#include <lac/sparse_matrix.h> 

const double PI =4.0*atan(1.0);
const u_int DIM=2;
class Reconstruct
{
public:
	Reconstruct(const std::string& file);
	~Reconstruct();
private:
	TemplateGeometry<DIM> triangle_template_geometry;
	CoordTransform<DIM,DIM> triangle_coord_transform;
	TemplateDOF<DIM> triangle_template_dof;
	BasisFunctionAdmin<double,DIM,DIM> triangle_basis_function;

	std::vector<TemplateElement<double,DIM,DIM> > template_element;
	FEMSpace<double,DIM> *fem_space;

//	Mesh<DIM,DIM> fixed_mesh;
	EasyMesh fixed_mesh;
	std::string mesh_file;

	BoundaryFunction<double,DIM> boundary;
	BoundaryConditionAdmin<double,DIM> boundary_admin;

	AMGSolver solver;

	FEMFunction<double,DIM> *u_h;
	Vector<double> rhs;//在不同的函数中用

	SparseMatrix<double> stiff_matrix;

	 double h_ele;

  	 double norm;
  	 double err;

public:
	void run();
	void initialize();
	void buildFEMSpace();
	
	void build_Ax();
	void build_rhs();
//	void NeummanBC();
//	void DirichletBC();

	void norm_u();
	void norm_gradient_u();
	void norm_gradient_u_R();
	void norm_laplace_u();


	void output();

	void barycenter(AFEPack::Point<DIM>& bc, GeometryBM& geo, Mesh<DIM, DIM>& mesh);


  //Reconstruct.cpp
 private:
 	Vector<double> recon_grad;
	Vector<double> recon_grad_cn;

	std::vector<std::vector<u_int> > reconstruct_patch;
	std::vector<double> weight_patch; //weight_patch[i] is the weight of ith_dof, 1/n_i_patch_ele 

						
  public:
	void testReconstruct();
	void buildReconstructPatch();
	
	};


 
 #endif //_RECONSTRUCT_H_

  

/*================================================================
*   Copyright (C) 2020 BHY. All rights reserved.
*   
*   文件名称：Reconstruct.cpp
*   创 建 者：HuanyingBIAN
*   创建日期：2020年07月13日
*   描    述：
*
================================================================*/


#include "Reconstruct.h"

// face size
double getFaceSize(const AFEPack::Point<DIM>& p0,
		   const AFEPack::Point<DIM>& p1,
		   const AFEPack::Point<DIM>& p2)
{
  double a,b,c,p,s,h;

  a = distance(p0,p1);
  b = distance(p2,p1);
  c = distance(p2,p0);

  p = 0.5*(a+b+c);
  s = sqrt(p*(p-a)*(p-b)*(p-c));
  h = 0.5*a*b*c/s;
  
  return h;
  
}




double f(const double * p)
{
//sin(PI*p[0]) * sin(2*PI*p[1])
	double a,b,c;
	a = 8*pow(PI,4)*(pow(sin(PI*p[0]),2)-pow(cos(PI*p[0]),2))*pow(sin(PI*p[1]),2);
	b = 8*pow(PI,4)*(pow(sin(PI*p[1]),2)-pow(cos(PI*p[1]),2))*pow(sin(PI*p[0]),2);
	c = 8*pow(PI,4)*(pow(sin(PI*p[0]),2)-pow(cos(PI*p[0]),2))*(pow(sin(PI*p[1]),2)-pow(cos(PI*p[1]),2));

  return ( a+b+c  );
};

double u(const double * p)
{
	return( pow(sin(PI*p[0]),2)*pow(sin(PI*p[1]),2) );
}

double g_1(const double * p)
{
	return( 0. );//bas_val|BC = 0,D_BC bu chu li
}

double g_2(const double * p)
{
	return( 0. );
}



int main(int argc, char *argv[]){
	Reconstruct myReconstruct(argv[1]);
	myReconstruct.run();

	return 0;
}

Reconstruct::Reconstruct(const std::string & file): mesh_file(file) {};
Reconstruct::~Reconstruct(){};

void Reconstruct::run()
{
	initialize();
	buildFEMSpace();

	build_rhs();
  std::cout << " rhs has passed"<<std::endl;
	build_Ax();
//	norm_u();
  std::cout << " _Ax has passed"<<std::endl;

	output();
	  std::cout << " output has passed"<<std::endl;
	
}

void Reconstruct::output()
{
  u_h->writeOpenDXData("u_h.dx");
};

void Reconstruct::initialize()
{
	triangle_template_geometry.readData("triangle.tmp_geo");
	triangle_coord_transform.readData("triangle.crd_trs");
	triangle_template_dof.reinit(triangle_template_geometry);
	triangle_template_dof.readData("triangle.1.tmp_dof");
	triangle_basis_function.reinit(triangle_template_dof);
	triangle_basis_function.readData("triangle.1.bas_fun");

 	template_element.resize(1);
 	template_element[0].reinit(triangle_template_geometry,
 								triangle_template_dof,
 								triangle_coord_transform,
 								triangle_basis_function);

 	fixed_mesh.readData(mesh_file);
}

void Reconstruct::buildFEMSpace()
{
	fem_space = new FEMSpace<double,DIM>(fixed_mesh,template_element);

	int n_element = fixed_mesh.n_geometry(DIM);
	fem_space->element().resize(n_element);

	for(int i = 0;i < n_element;++ i){
		fem_space->element(i).reinit(*fem_space,i,0);
	}

	fem_space->buildElement();
	fem_space->buildDof();
	fem_space->buildDofBoundaryMark();

std::cout << "build  initialize  has passed"<<std::endl;
	/******接下来根据单元构造重构片***/
	buildReconstructPatch();
}

void Reconstruct::build_rhs()
{
	rhs.reinit(fem_space->n_dof());
	
 	FEMSpace<double,DIM>::ElementIterator the_ele = fem_space->beginElement();
  	FEMSpace<double,DIM>::ElementIterator end_ele = fem_space->endElement();
	for(;the_ele != end_ele;++ the_ele){
		double volume = the_ele->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_ele->findQuadratureInfo(3); //3 is accuracy
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<double> jacobian = the_ele->local_to_global_jacobian(quad_info.quadraturePoint());
		std::vector<AFEPack::Point<DIM> > q_point = the_ele->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> > bas_val = the_ele->basis_function_value(q_point);
		
		u_int n_ele_dof = the_ele->dof().size();

		for(u_int l = 0;l < n_quadrature_point;++ l){
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			for (u_int i = 0;i < n_ele_dof;++ i) {
			//  std::cout << 30<<std::endl;
				rhs(the_ele->dof()[i]) +=  Jxw*f(q_point[l])*bas_val[i][l] ;
			//  std::cout << 31<<std::endl;
			  }
		//  std::cout << 33<<std::endl;
		}
		//  std::cout <<"ele"<< 31<<std::endl;
	}

}

void Reconstruct::build_Ax()
{
	u_h = new FEMFunction<double,DIM> (*fem_space);

	int n_total_dof = fem_space->n_dof();

		FullMatrix<double> stiff_A;
		FullMatrix<double> stiff_B;
		FullMatrix<double> stiff_C;
		FullMatrix<double> stiff_D;

		stiff_A.reinit(n_total_dof,n_total_dof);
		stiff_B.reinit(n_total_dof,n_total_dof);
		stiff_C.reinit(n_total_dof,n_total_dof);
		stiff_D.reinit(n_total_dof,n_total_dof);


	FullMatrix<double> A(n_total_dof,n_total_dof);
/*
std::cout << "n_total_dof = "<<n_total_dof<< "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
getchar();
*/


	FullMatrix<double> stiff_P;
	FullMatrix<double> stiff_Q;
	FullMatrix<double> stiff_S;
	FullMatrix<double> stiff_T;
	stiff_P.reinit(n_total_dof,n_total_dof);
	stiff_Q.reinit(n_total_dof,n_total_dof);
	stiff_S.reinit(n_total_dof,n_total_dof);
	stiff_T.reinit(n_total_dof,n_total_dof);

	FullMatrix<double> stiff_P_BC;
	stiff_P_BC.reinit(n_total_dof,n_total_dof);

	for(u_int i = 0;i < n_total_dof;++ i){
	  	for(u_int j = 0;j < n_total_dof;++ j){
			stiff_A[i][j] = 0.;// A.set(1,2,3.141); is entirely equivalent to the operation A(1,2) = 3.141;
			stiff_B[i][j] = 0.;
			stiff_C[i][j] = 0.;
			stiff_D[i][j] = 0.;

			stiff_P[i][j] = 0.;
			stiff_Q[i][j] = 0.;
			stiff_S[i][j] = 0.;
			stiff_T[i][j] = 0.;

			stiff_P_BC[i][j] = 0.;
	  	}
	}

	FEMSpace<double,DIM>::ElementIterator the_ele = fem_space->beginElement();
	FEMSpace<double,DIM>::ElementIterator end_ele = fem_space->endElement();

	for(;the_ele != end_ele;++ the_ele){
		double volume = the_ele->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_ele->findQuadratureInfo(3); //3 is accuracy
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<double> jacobian = the_ele->local_to_global_jacobian(quad_info.quadraturePoint());
		std::vector<AFEPack::Point<DIM> > q_point = the_ele->local_to_global(quad_info.quadraturePoint());

		std::vector<std::vector<std::vector<double> > > basis_grad = the_ele->basis_function_gradient(q_point);
		std::vector<std::vector<double> >  bas_val = the_ele->basis_function_value(q_point);

		const std::vector<int>& ele_dof = the_ele->dof();
		u_int n_ele_dof = ele_dof.size();


		const int the_ele_idx = the_ele->index();
		const AFEPack::Point<DIM>& p0 = fixed_mesh.point(fixed_mesh.geometry(2,the_ele_idx).vertex(0));//regular_mesh.geometry(DIM,ele_idx).vertex(0) 当前网格下DIM维几何体第ele_idx个几何体的第0个点，.point()函数是取指定点的坐标
        const AFEPack::Point<DIM>& p1 = fixed_mesh.point(fixed_mesh.geometry(2,the_ele_idx).vertex(1));
        const AFEPack::Point<DIM>& p2 = fixed_mesh.point(fixed_mesh.geometry(2,the_ele_idx).vertex(2));
		h_ele = getFaceSize(p0,p1,p2);

/*
std::cout << "h_ele = "<<h_ele<< "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
getchar(); 	*/

		for(u_int l = 0;l < n_quadrature_point;++ l){
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;

			for(u_int i = 0;i < n_ele_dof;++ i){
				for(u_int j = 0;j < n_ele_dof;++ j){
				stiff_P[i][j] += Jxw*( basis_grad[i][l][0]*basis_grad[j][l][0] );
				stiff_Q[i][j] += Jxw*( basis_grad[i][l][0]*basis_grad[j][l][1] );
				stiff_S[i][j] += Jxw*( basis_grad[i][l][1]*basis_grad[j][l][0] );
				stiff_T[i][j] += Jxw*( basis_grad[i][l][1]*basis_grad[j][l][1] );

				stiff_P_BC[i][j] += 1.0/h_ele/h_ele*Jxw*( bas_val[i][l]*bas_val[j][l] );//1/h^2,h是单元的尺寸
				}
			}
		}
	}

std::cout << "Press ENTER to continue or CTRL+C to stop ..."<<" h_ele = "<<h_ele<<std::endl;

// 	Mesh<DIM,DIM>& mesh = fixed_mesh;
  	u_int n_ele = fixed_mesh.n_geometry(DIM);

	int n_dof_ = fem_space->n_dof();

	for (u_int i = 0;i < n_ele;i ++) {//对所有的单元做遍历
		GeometryBM& ele = fixed_mesh.geometry(DIM, i);
		const std::vector<int>& vtx = ele.vertex();
		u_int n_vtx = vtx.size();

		std::cout << "2 Press ENTER to continue or CTRL+C to stop ..."<<std::endl;


		for (u_int j = 0;j < n_vtx;j ++) {//对该单元上所有的节点遍历

			const std::vector<u_int>& rp = reconstruct_patch[vtx[j]];
			//the elements are the neighborhood of vtx[j]. vtx[j]与自由度的编号一一对应吗？？？？？？？？
			int N_patch_ele = rp.size();
				
				std::cout << "N_patch_ele ="<<  N_patch_ele <<" j = "<< j <<" vtx[j] = "<< vtx[j]<<std::endl;
			//	std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
			 //   getchar();
			    

			for(int l = 0; l < N_patch_ele; ++ l)
			{
				Element<double, DIM>& ele1 = fem_space->element(reconstruct_patch[vtx[j]][l]);
				const std::vector<int>& ele1_dof = ele1.dof();
				int n_ele1_dof = ele1_dof.size();

				for(u_int k = 0;k < n_ele1_dof;++ k)//想找到自由度所对应的坐标，以便求出基函数的导数值
				{
					const int& dof = ele1_dof[k];
      				const AFEPack::Point<DIM>& interp_point = fem_space->dofInfo(dof).interp_point;/// 自由度dof的坐标
					std::vector<double> basis_grade_ = ele1.basis_function_gradient(k,interp_point); //j=element_dof[k],A(i,j)
				//			std::vector<value_t> basis_function_gradient(int i, const Point<DOW>&) const; < Gradient of certain basis function at a point. i= 0,1,2 线性元的基函数不连续，只有在特定单元上有定义，所以无法识别全局自由度的编号

				stiff_A[vtx[j]][ele1_dof[k]] += weight_patch[vtx[j]]*basis_grade_[0];
				stiff_B[vtx[j]][ele1_dof[k]] += weight_patch[vtx[j]]*basis_grade_[1];
					
				/*	std::cout << "vtx[j] patch l ele k_dof of l"<<" k = "<< k<< "l n_ele1_dof = " << n_ele1_dof <<"vtx[j] = "<<vtx[j]<<"N_patch_ele[i]= "<< N_patch_ele<< "l =  "<< l <<std::endl;

					std::cout << "ele1_dof[k] = "<< ele1_dof[k]<<" k = "<< k <<std::endl;
					std::cout << "stiff_A[vtx[j]][ele1_dof[k]] = "<< stiff_A[vtx[j]][ele1_dof[k]] <<std::endl;
					*/	

				}
			}
		}
	}

	for(u_int i = 0;i < n_total_dof;++ i){
	  	for(u_int j = 0;j < n_total_dof;++ j){
			stiff_C[i][j] = stiff_A[j][i];
			stiff_D[i][j] = stiff_B[j][i];
	  	}
	}

	stiff_C.mmult(stiff_P,stiff_P);//stiff_P_1=stiff_C*stiff_P_2
	stiff_P.mmult(stiff_P,stiff_A);//stiff_P_()=stiff_P*stiff_A=stiff_C*stiff_P_2*stiff_A

	stiff_C.mmult(stiff_Q,stiff_Q);//stiff_Q_1=stiff_C*stiff_Q_2
	stiff_Q.mmult(stiff_Q,stiff_B);//stiff_Q_()=stiff_Q*stiff_B=stiff_C*stiff_Q_2*stiff_B

	stiff_D.mmult(stiff_S,stiff_S);//stiff_S_1=stiff_D*stiff_S_2
	stiff_S.mmult(stiff_S,stiff_A);//stiff_S_()=stiff_S*stiff_A=stiff_D*stiff_S_2*stiff_A

	stiff_D.mmult(stiff_T,stiff_T);//stiff_T_1=stiff_D*stiff_T_2
	stiff_T.mmult(stiff_T,stiff_B);//stiff_T_()=stiff_T*stiff_B=stiff_D*stiff_T_2*stiff_B

	A.equ(1,stiff_P,1,stiff_Q,1,stiff_S);//A=1*stiff_P+1*stiff_Q+1*stiff_S=CPA+CQB+DSA
	A.equ(1,A,1,stiff_T);//A=1*A+1*stiff_T=CPA+CQB+DSA+DTB

//////////////////////////////////////////////////////
/// 下面手工处理边界条件。NBC
//////////////////////////////////////////////////////////
	FullMatrix<double> stiff_A_2(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_A_4(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_B_1(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_B_3(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_C_2(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_C_4(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_D_1(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_D_3(n_total_dof,n_total_dof);

	FullMatrix<double> stiff_I_1(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_I_2(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_I_3(n_total_dof,n_total_dof);
	FullMatrix<double> stiff_I_4(n_total_dof,n_total_dof);


	for(u_int i = 0;i < n_total_dof;++ i){
	  	for(u_int j = 0;j < n_total_dof;++ j){
			stiff_A_2[i][j] = 0.;// A.set(1,2,3.141); is entirely equivalent to the operation A(1,2) = 3.141;
			stiff_A_4[i][j] = 0.;
			stiff_B_1[i][j] = 0.;
			stiff_B_3[i][j] = 0.;

			stiff_C_2[i][j] = 0.;
			stiff_C_4[i][j] = 0.;
			stiff_D_1[i][j] = 0.;
			stiff_D_3[i][j] = 0.;

			stiff_I_1[i][j] = 0.;
			stiff_I_2[i][j] = 0.;
			stiff_I_3[i][j] = 0.;
			stiff_I_4[i][j] = 0.;
	  	}
	}

	for(u_int i = 0;i < n_total_dof;++ i)
	{
		const std::vector<u_int>& rp = reconstruct_patch[i];
			//the elements are the neighborhood of i-th dof 
		int N_patch_ele = rp.size();

		int bm = fem_space->dofInfo(i).boundary_mark;
	/// 若该自由度标识为 1，这个是在生成网格时手工指定的。
    	if (bm == 1)
    	{
			for(int l = 0; l < N_patch_ele; ++ l)
			{
				Element<double, DIM>& theEle = fem_space->element(reconstruct_patch[i][l]);
				
				const std::vector<int>& element_dof = theEle.dof();
				int n_element_dof = theEle.n_dof();	
				//int n_element_dof = element_dof.size();

				for(u_int k = 0;k < n_element_dof;++ k)
				{
					stiff_B_1[i][element_dof[k]] = stiff_B[i][element_dof[k]];
					stiff_D_1[i][element_dof[k]] = stiff_D[i][element_dof[k]];

					stiff_I_1[i][element_dof[k]] = 1.;
				}

			}
		}
		else if(bm == 2)
    	{
			for(int l = 0; l < N_patch_ele; ++ l)
			{
				Element<double, DIM>& theEle = fem_space->element(reconstruct_patch[i][l]);
				
				const std::vector<int>& element_dof = theEle.dof();
				int n_element_dof = theEle.n_dof();	

				for(u_int k = 0;k < n_element_dof;++ k)
				{
					stiff_A_2[i][element_dof[k]] = stiff_A[i][element_dof[k]];
					stiff_C_2[i][element_dof[k]] = stiff_D[i][element_dof[k]];

					stiff_I_2[i][element_dof[k]] = 1.;
				}

			}
		}
		else if (bm == 3)
    	{
			for(int l = 0; l < N_patch_ele; ++ l)
			{
				Element<double, DIM>& theEle = fem_space->element(reconstruct_patch[i][l]);
				
				const std::vector<int>& element_dof = theEle.dof();
				int n_element_dof = theEle.n_dof();	

				for(u_int k = 0;k < n_element_dof;++ k)
				{
					stiff_B_3[i][element_dof[k]] = stiff_B[i][element_dof[k]];
					stiff_D_3[i][element_dof[k]] = stiff_D[i][element_dof[k]];

					stiff_I_3[i][element_dof[k]] = 1.;
				}

			}
		}
		else if(bm == 4)
    	{
			for(int l = 0; l < N_patch_ele; ++ l)
			{
				Element<double, DIM>& theEle = fem_space->element(reconstruct_patch[i][l]);
				
				const std::vector<int>& element_dof = theEle.dof();
				int n_element_dof = theEle.n_dof();	

				for(u_int k = 0;k < n_element_dof;++ k)
				{
					stiff_A_4[i][element_dof[k]] = stiff_A[i][element_dof[k]];
					stiff_C_4[i][element_dof[k]] = stiff_D[i][element_dof[k]];

					stiff_I_4[i][element_dof[k]] = 1.;
				}
			}
		}
	}


	stiff_D_1.mmult(stiff_D_1,stiff_P_BC);//stiff_D_1_1=stiff_D_1*stiff_P_BC
	stiff_D_1.mmult(stiff_D_1,stiff_B_1);//stiff_D_1_()=stiff_D_1*stiff_B_1=stiff_D_1*stiff_P_BC*stiff_B_1

	stiff_D_3.mmult(stiff_D_3,stiff_P_BC);//stiff_D_3_1=stiff_D_3*stiff_P_BC
	stiff_D_3.mmult(stiff_D_3,stiff_B_3);//stiff_D_3_()=stiff_D_3*stiff_B_3=stiff_D_3*stiff_P_BC*stiff_B_3

	stiff_C_2.mmult(stiff_C_2,stiff_P_BC);//stiff_C_2_1=stiff_C_2*stiff_P_BC
	stiff_C_2.mmult(stiff_C_2,stiff_A_2);//stiff_C_2_()=stiff_C_2*stiff_A_2=stiff_C_2*stiff_P_BC*stiff_A_2

	stiff_C_4.mmult(stiff_C_4,stiff_P_BC);//stiff_C_4_1=stiff_C_4*stiff_P_BC
	stiff_C_4.mmult(stiff_C_4,stiff_A_4);//stiff_C_4_()=stiff_C_4*stiff_A_4=stiff_C_4*stiff_P_BC*stiff_A_4

	A.equ(1.,stiff_D_1,1.,stiff_C_2,1.,stiff_D_3);//A=1*stiff_D_1+1*stiff_C_2+1*stiff_D_3=D_1*P_BC*B_1+C_2*P_BC*A_2+D_3*P_BC*B_3
	A.equ(1.,A,1.,stiff_C_4);//A=1*A+1*stiff_C_4=D_1*P_BC*B_1+C_2*P_BC*A_2+D_3*P_BC*B_3+C_4*P_BC*A_4

//////////////////////////////////////////////////////
/// 下面手工处理边界条件。NBC
//////////////////////////////////////////////////////////

	stiff_I_1.mmult(stiff_I_1,stiff_P_BC);//stiff_I_1_1=stiff_I_1*stiff_P_BC
	stiff_I_3.mmult(stiff_I_3,stiff_P_BC);//stiff_I_1_1=stiff_I_1*stiff_P_BC
	stiff_I_2.mmult(stiff_I_2,stiff_P_BC);//stiff_I_1_1=stiff_I_1*stiff_P_BC
	stiff_I_4.mmult(stiff_I_4,stiff_P_BC);//stiff_I_1_1=stiff_I_1*stiff_P_BC

	A.equ(1./h_ele/h_ele,stiff_I_1,1/h_ele/h_ele,stiff_I_2,1/h_ele/h_ele,stiff_I_3);//A=1*stiff_D_1+1*stiff_C_2+1*stiff_D_3=D_1*P_BC*B_1+C_2*P_BC*A_2+D_3*P_BC*B_3
	A.equ(1.,A,1.,stiff_I_4);//A=1.0/h_ele/h_ele*I_1*P_BC+1.0/h_ele/h_ele*I_2*P_BC+1.0/h_ele/h_ele*I_3*P_BC+1.0/h_ele/h_ele*I_4*P_BC


	/****
	***A to SparseMatrix<double> stiff_matrix 
	******/
	
	SparsityPattern sp_stiff_matrix;
	sp_stiff_matrix.copy_from(A);
	//void SparsityPattern::copy_from	(	const FullMatrix< number > & 	matrix	)
	//	Take a full matrix and use its nonzero entries to generate a sparse matrix entry pattern for this object.
	//Previous content of this object is lost, and the sparsity pattern is in compressed mode afterwards.

	/// 矩阵模板压缩. 创建矩阵.
    sp_stiff_matrix.compress();
    SparseMatrix<double> stiff_matrix(sp_stiff_matrix);
/*//给出SparseMatrix的尺寸m \times n
std::cout << "m = "<<stiff_matrix.m()<<std::endl;
std::cout << "n = "<<stiff_matrix.n()<<std::endl;
std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
getchar();
*/

	/****
	***A to SparseMatrix<double> stiff_matrix 
	******/
	solver.reinit(stiff_matrix);
	solver.solve(*u_h, rhs, 1e-08, 50);	

}

void Reconstruct::norm_u()
{
	norm =  0.0;
	FEMSpace<double,DIM>::ElementIterator the_element = fem_space->beginElement();
	FEMSpace<double,DIM>::ElementIterator end_element = fem_space->endElement();
	for (;the_element != end_element;the_element ++) {
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<AFEPack::Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());

		std::vector<double> u_h_val = u_h->value(q_point, *the_element);

		for (int l = 0;l < n_quadrature_point;l ++) {
			double Jxw = quad_info.weight(l)*jacobian[l]*volume;
			norm += Jxw*(u_h_val[l]-u(q_point[l]))*(u_h_val[l]-u(q_point[l]));
		}
	}
	err = sqrt(fabs(norm));

}

void Reconstruct::norm_gradient_u()
{


}

void Reconstruct::norm_gradient_u_R()
{

}

void Reconstruct::norm_laplace_u()
{

}






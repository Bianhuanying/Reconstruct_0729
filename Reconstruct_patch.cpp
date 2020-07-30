/*================================================================
*   Copyright (C) 2020 BHY. All rights reserved.
*   
*   文件名称：Reconstruct_patch.cpp
*   创 建 者：HuanyingBIAN
*   创建日期：2020年07月13日
*   描    述：
*
================================================================*/


#include "Reconstruct.h"

void centroid_point(AFEPack::Point<DIM>& bc, GeometryBM& geo, Mesh<DIM, DIM>& mesh)
{
  const int& n_vtx = geo.n_vertex();
  for(int i = 0;i < n_vtx;++ i){
    const int& vtx_idx = geo.vertex(i);
    AFEPack::Point<DIM>& vtx = mesh.point(vtx_idx);
    bc[0] += vtx[0];     bc[1] += vtx[1];
  }
  bc[0] /= n_vtx; bc[1] /= n_vtx;
}


void Reconstruct::buildReconstructPatch()
{
  reconstruct_patch.clear();

 /// 使用有公共顶点的所有单元作为重构片
  /**
   * 这里我们使用了缺省的假设“几何体单元的编号和有限元单元的编号一致”，
   * 潜在地会出现当这个情况不对的时候的问题。
   * 
   */
  //Mesh<DIM,DIM>& mesh = irregular_mesh->regularMesh();
  Mesh<DIM,DIM>& mesh = fixed_mesh;
  u_int n_pnt = mesh.n_geometry(0);
  u_int n_ele = mesh.n_geometry(DIM);

 // std::vector<std::list<u_int> > ele_of_pnt(n_pnt); // 记录与每个点相交的所在的单元集合
 // std::vector<std::vector<u_int> > reconstruct_patch;//defined in ".h"
	int n_dof_ = fem_space->n_dof();
	reconstruct_patch.resize(n_dof_);//n_dof_ = n_pnt

  for (u_int i = 0;i < n_ele;i ++) {
    GeometryBM& ele = mesh.geometry(DIM, i);
    const std::vector<int>& vtx = ele.vertex();
    u_int n_vtx = vtx.size();

    for (u_int j = 0;j < n_vtx;j ++) {
      reconstruct_patch[vtx[j]].push_back(i);

    /*
    int k_ = reconstruct_patch[vtx[j]].size();
    for(int i_=0;i_<k_;++i_){
    std::cout << "reconstruct_patch[vtx[j]][i_] = "<<reconstruct_patch[vtx[j]][i_]<<" k_ = "<< k_<<" vtx[j] = "<< vtx[j] << " j = "<<j<<std::endl;}
    std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
    getchar();
    */ /*
      int k_ = reconstruct_patch[vtx[j]].size();
      for(int i_=0;i_<k_;++i_){
        if(vtx[j]==0)
          std::cout << "reconstruct_patch[vtx[j]][i_] = "<<reconstruct_patch[vtx[j]][i_]<<" k_ = "<< k_<<" vtx[j] = "<< vtx[j] << " j = "<<j<<std::endl;
      }
      */
     
    }

  }


  ///*******weights for linear FEM dof = n_pnt**********/
	weight_patch.resize(n_pnt);
  //所有的vector在用之前先def 尺寸，否则Segmentation fault (core dumped)

  for(u_int j_=0; j_<n_pnt; j_++){
    int n_patch_ele = reconstruct_patch[j_].size();
    weight_patch[j_] = 1.0/n_patch_ele;
  }

}


/**
 * end of file
 * 
 */

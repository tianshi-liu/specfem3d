#include "mesh_constants_gpu.h"

extern EXTERN_LANG
void FC_FUNC_(wavefield_discontinuity_add_traction_cuda,
              WAVEFIELD_DISCONTINUITY_ADD_TRACTION_CUDA)(int* size_points, 
                                                         int* size_faces,
                                                         long* Mesh_pointer){
  TRACE("wavefield_discontinuity_add_traction_cuda");
  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  if (mp->is_wavefield_discontinuity) {
    int size = (*size_points);
    int blocksize = BLOCKSIZE_KERNEL1;
    int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);
    
#ifdef USE_CUDA
    if (run_cuda) {
      add_acceleration_discontinuity_kernel
       <<<grid, threads, 0, mp->compute_stream>>>(mp->d_accel_wd,
                                                  mp->d_mass_in_wd,
                                                  mp->d_boundary_to_iglob_wd,
                                                  size, mp->d_accel);
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
      hipLaunchKernelGGL(add_acceleration_discontinuity_kernel, 
       dim3(grid), dim3(threads), 0, mp->compute_stream, mp->d_accel_wd,
                                                  mp->d_mass_in_wd,
                                                  mp->d_boundary_to_iglob_wd,           
                                                  size, mp->d_accel);
    }
#endif
    
    size = (*size_faces);
    blocksize = NGLL2;
    
    get_blocks_xy(size,&num_blocks_x,&num_blocks_y);

    dim3 grid2(num_blocks_x,num_blocks_y);
    dim3 threads2(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda) {
       add_traction_discontinuity_kernel
        <<<grid2, threads2, 0, mp->compute_stream>>>(mp->d_traction_wd,
                                                   mp->d_face_ispec_wd,
                                                   mp->d_face_ijk_wd,
                                                   mp->d_face_jacobian2Dw_wd,
                                                   mp->d_ibool,
                                                   size, mp->d_accel);
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
      hipLaunchKernelGGL(add_traction_discontinuity_kernel,
        dim3(grid2), dim3(threads2), 0, mp->compute_stream, mp->d_traction_wd,
                                                   mp->d_face_ispec_wd,
                                                   mp->d_face_ijk_wd,
                                                   mp->d_face_jacobian2Dw_wd,
                                                   mp->d_ibool,
                                                   size, mp->d_accel);
    }
#endif
  }
}

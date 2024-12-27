__global__ void add_acceleration_discontinuity_kernel(
                                  realw_const_p accel_wd, 
                                  realw_const_p mass_in_wd,
                                  const int* boundary_to_iglob_wd,
                                  const int size, realw* accel
                                  ) {
  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;
  int iglob = boundary_to_iglob_wd[id] - 1;
  realw mass_in = mass_in_wd[id];
  if (id < size) {
    accel[iglob*3] = accel[iglob*3] - accel_wd[id*3] * mass_in;
    accel[iglob*3 + 1] = accel[iglob*3 + 1] - accel_wd[id*3 + 1] * mass_in;
    accel[iglob*3 + 2] = accel[iglob*3 + 2] - accel_wd[id*3 + 2] * mass_in;
  }
}

__global__ void add_traction_discontinuity_kernel(
                                  realw_const_p traction_wd,
                                  const int* face_ispec_wd,
                                  const int* face_ijk_wd,
                                  realw_const_p face_jacobian2Dw_wd,
                                  const int* d_ibool,
                                  const int size, realw* accel) {
  int igll = threadIdx.x;
  int iface_wd = blockIdx.x + gridDim.x*blockIdx.y;
  int i, j, k, ispec, iglob;
  realw jacobianw;
  if (iface_wd < size) {
    ispec = face_ispec_wd[iface_wd] - 1;
    i = face_ijk_wd[INDEX3(NDIM,NGLL2,0,igll,iface_wd)]-1;
    j = face_ijk_wd[INDEX3(NDIM,NGLL2,1,igll,iface_wd)]-1;
    k = face_ijk_wd[INDEX3(NDIM,NGLL2,2,igll,iface_wd)]-1;

    iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    jacobianw = face_jacobian2Dw_wd[INDEX2(NGLL2,igll,iface_wd)];
    atomicAdd(&accel[iglob*3],  traction_wd[INDEX3(NDIM,NGLL2,0,igll,iface_wd)] * jacobianw);
    atomicAdd(&accel[iglob*3+1],  traction_wd[INDEX3(NDIM,NGLL2,1,igll,iface_wd)] * jacobianw);
    atomicAdd(&accel[iglob*3+2],  traction_wd[INDEX3(NDIM,NGLL2,2,igll,iface_wd)] * jacobianw);
  }
}

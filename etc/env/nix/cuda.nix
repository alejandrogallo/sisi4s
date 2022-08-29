{pkgs}:
{

  shellHook = ''
    export NVCC=${pkgs.cudatoolkit}/bin/nvcc

    export CUDA_ROOT_PATH=${pkgs.cudatoolkit}
    export CUDA_CXXFLAGS="-I${pkgs.openmpi.out}/include -I$CUDA_ROOT_PATH/include"

    # cudart
    export CUDA_LIB_PATH=${pkgs.cudatoolkit.lib}
    export CUDA_OUT_PATH=${pkgs.cudatoolkit.out}

    # in here is libcuda
    export CUDA_X11LIB="${pkgs.linuxPackages.nvidia_x11}/lib"
    export CUDA_LDFLAGS="-L$CUDA_X11LIB -L$CUDA_ROOT_PATH/lib -L$CUDA_LIB_PATH/lib -lcuda -lcudart -lcublas"
  '';

}

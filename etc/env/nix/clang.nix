{ pkgs, ...}:

{
  buildInputs = with pkgs; [
    clang
    llvmPackages.openmp
  ];

  shellHook = ''
    export CC=clang
    export CXX=clang++
    export OMPI_CC=clang
    export OMPI_CXX=clang++
  '';
}


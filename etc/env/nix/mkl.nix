{ pkgs, ...}:

{
  buildInputs = with pkgs; [
    mkl
  ];

  shellHook = ''
    export MKL_PATH=${pkgs.mkl}
    export LD_LIBRARY_PATH=$MKL_PATH:$LD_LIBRARY_PATH
  '';
}


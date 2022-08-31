{ pkgs, ... }:

let

  /*
  myopenblas = pkgs.openblas.overrideDerivation (old:
    {
      dontStrip = true;
      doCheck = false;
      CFLAGS = "-g";
      FCLAGS = "-g";
    }
  );
  */

  # myopenblas = pkgs.enableDebugging pkgs.openblas;
  myopenblas = pkgs.openblas;

in

{

  buildInputs = with pkgs; [
    myopenblas
    scalapack
  ];

  shellHook = ''
    export OPENBLAS_PATH=${myopenblas}
    export SCALAPACK_PATH=${pkgs.scalapack}
    export LD_LIBRARY_PATH=${pkgs.scalapack}/lib:$LD_LIBRARY_PATH
  '';

}

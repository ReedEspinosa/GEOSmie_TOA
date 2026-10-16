# macOS with Homebrew gfortran: the linker needs the SDK path
if [ "$(uname)" = "Darwin" ] && [ -z "$SDKROOT" ]; then export SDKROOT=$(xcrun --show-sdk-path); fi

# compile gsf
gfortran -o gsf/a.out gsf/spher_expan.f

# compile the GRASP sphere-kernel generator (graspkernels/; see graspkernels/README.md)
gfortran -O2 -o graspkernels/spher_tri graspkernels/spher_tri.f graspkernels/spher.f
gfortran -O2 -o graspkernels/spher_main graspkernels/spher_main.f graspkernels/spher.f

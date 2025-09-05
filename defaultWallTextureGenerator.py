import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Default Parameters: 
# ---------------------------
totalNumberOfAtomTypes = 4
totallength = 7.0
textureAtomtype = 2
textureRadius = 0.00159
textureDensity = 2650
boxxHalfSize = 0.03
boxyHalfSize = boxxHalfSize
boxzlo = -2e-15
boxzhi = totallength
skipLength = 0.1
output_dir = "./"

def defaultWallTextureGenerator(
                            totalNumberOfAtomTypes, 
                            totallength,
                            textureAtomtype,
                            textureRadius,
                            textureDensity,
                            boxxHalfSize,
                            boxyHalfSize,
                            boxzlo,
                            boxzhi,
                            output_dir):
    # ---------------------------
    # Export 1: walltexture.dat
    # ---------------------------
    target_path = output_dir + "/walltexture.dat"
    with open(target_path, "w") as f:
        f.write("LIGGGHTS data file\n\n")
        f.write(f"1 atoms\n\n")  #movable atom
        f.write(f"{totalNumberOfAtomTypes} atom types\n\n")
        f.write(f"{-boxxHalfSize:.9f} {boxxHalfSize:.9f} xlo xhi\n")
        f.write(f"{-boxyHalfSize:.9f} {boxyHalfSize:.9f} ylo yhi\n")
        f.write(f"{boxzlo} {boxzhi:.9f} zlo zhi\n\n")
        f.write("Atoms\n\n")
        f.write(f"{1} {textureAtomtype} {textureRadius*2:.9f} " f"{textureDensity:.9f} 0 0 0\n")
        print("write file: {}".format(target_path))
    # ---------------------------
    # Export 2: walltexture.in
    # ---------------------------
    target_path = output_dir + "/walltexture.in"
    with open(target_path, "w") as f:
        print("write file: {}".format(target_path))

    return

if __name__ == "__main__":
    defaultWallTextureGenerator( totalNumberOfAtomTypes, 
                            totallength,
                            textureAtomtype,
                            textureRadius,
                            textureDensity,
                            boxxHalfSize,
                            boxyHalfSize,
                            boxzlo,
                            boxzhi,
                            output_dir)

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Default Parameters: 
# ---------------------------
totalNumberOfAtomTypes = 4
amplitude = 0.02815
period = 0.636
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

def helixWallTextureGenerator(
                            totalNumberOfAtomTypes, 
                            amplitude,
                            period,
                            totallength,
                            textureAtomtype,
                            textureRadius,
                            textureDensity,
                            boxxHalfSize,
                            boxyHalfSize,
                            boxzlo,
                            boxzhi,
                            skipLength,
                            output_dir):
     # ---------------------------
    # Derived values
    # ---------------------------
    distanceBetweenTextureParticles = (textureRadius * 2) * 1.05
    amplitude = amplitude - (textureRadius * 1.05)

    # Helix geometry
    arcDistance = 2 * np.pi * amplitude
    diagonalDistance = np.sqrt(period**2 + arcDistance**2)
    numPointsPerPeriod = int(diagonalDistance / distanceBetweenTextureParticles)

    u = np.linspace(0, 2*np.pi, numPointsPerPeriod)
    temp_x = amplitude * np.cos(u)
    temp_y = amplitude * np.sin(u)
    temp_z = u * period / (2*np.pi)

    x, y, z = temp_x.copy(), temp_y.copy(), temp_z.copy()

    # Repeat helix for total length
    for i in range(1, int((totallength - period) / period) + 1):
        x = np.concatenate((x, temp_x[1:]))
        y = np.concatenate((y, temp_y[1:]))
        z = np.concatenate((z, temp_z[1:] + i * period))

    # Skip particles below skipLength
    indices = np.where(z > skipLength)[0]
    if indices.size > 0:
        firstIndex = indices[0]
        print(f"The first index is: {firstIndex+1}")  # +1 to match MATLAB indexing
        x, y, z = x[firstIndex+1:], y[firstIndex+1:], z[firstIndex+1:]

    # ---------------------------
    # Plot helix
    # ---------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(x, y, z, 'r')
    ax.set_box_aspect([1, 1, 1])  # equal aspect ratio
    plt.show()

    # ---------------------------
    # Export 1: walltexture.dat
    # ---------------------------
    target_path = output_dir + "/walltexture.dat"
    with open(target_path, "w") as f:
        f.write("LIGGGHTS data file\n\n")
        f.write(f"{len(x)} atoms\n\n")
        f.write(f"{totalNumberOfAtomTypes} atom types\n\n")
        f.write(f"{-boxxHalfSize:.9f} {boxxHalfSize:.9f} xlo xhi\n")
        f.write(f"{-boxyHalfSize:.9f} {boxyHalfSize:.9f} ylo yhi\n")
        f.write(f"{boxzlo} {boxzhi:.9f} zlo zhi\n\n")
        f.write("Atoms\n\n")
        for i in range(len(x)):
            f.write(f"{i+1} {textureAtomtype} {textureRadius*2:.9f} "
                    f"{textureDensity:.9f} {x[i]:.9f} {y[i]:.9f} {z[i]:.9f}\n")
        print("write file: {}".format(target_path))
    # ---------------------------
    # Export 2: walltexture.in
    # ---------------------------
    target_path = output_dir + "/walltexture.in"
    with open(target_path, "w") as f:
        for i in range(len(x)):
            f.write(f"create_atoms {textureAtomtype} single "
                    f"{x[i]:.9f} {y[i]:.9f} {z[i]:.9f}\n")
        f.write(f"set type {textureAtomtype} diameter {textureRadius*2:.9f}\n")
        f.write(f"set type {textureAtomtype} density {textureDensity:.9f}\n")
        print("write file: {}".format(target_path))
    
    return

if __name__ == "__main__":
    helixWallTextureGenerator( totalNumberOfAtomTypes, 
                            amplitude,
                            period,
                            totallength,
                            textureAtomtype,
                            textureRadius,
                            textureDensity,
                            boxxHalfSize,
                            boxyHalfSize,
                            boxzlo,
                            boxzhi,
                            skipLength,
                            output_dir)
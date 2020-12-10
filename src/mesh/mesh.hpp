#ifndef DG_MESH_HPP
#define DG_MESH_HPP

namespace DGHydro {

  class ConfigFile;

  class Mesh {
  public:
    Mesh(ConfigFile *cf);
    ~Mesh(void);

    void Decompose(int rank, int num_proc);

    double *x, *y, *z;
    double dx, dy, dz;

    double minX, maxX;
    double minY, maxY;
    double minZ, maxZ;

    int Nx, Ny, Nz, nGhost;

    int startX, endX;
    int startY, endY;
    int startZ, endZ;

  private:
  };

} // namespace DGHydro

#endif  // DG_MESH_HPP

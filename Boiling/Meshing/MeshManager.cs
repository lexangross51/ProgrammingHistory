using Boiling.Meshing.Interfaces;

namespace Boiling.Meshing;

public class MeshManager(IMeshBuilder meshBuilder)
{
    public Mesh CreateMesh()
    {
        meshBuilder.CreatePoints();
        meshBuilder.CreateElements();
        meshBuilder.CreateBoundaries();
        return meshBuilder.GetMesh();
    }
}
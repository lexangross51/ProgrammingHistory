using Boiling.Meshing.Interfaces;

namespace Boiling.Meshing;

public class MeshManager
{
    public IMeshBuilder? MeshBuilder { get; set; }

    public MeshManager() { }
    
    public MeshManager(IMeshBuilder meshBuilder)
        => MeshBuilder = meshBuilder;

    public Mesh CreateMesh()
    {
        MeshBuilder!.CreatePoints();
        MeshBuilder.CreateElements();
        MeshBuilder.CreateBoundaries();
        return MeshBuilder.GetMesh();
    }
}
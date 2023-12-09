using Boiling.Meshing.Geometry;

namespace Boiling.Meshing.Interfaces;

public interface IMeshBuilder
{
    IEnumerable<Point> Points { get; set; }
    void CreatePoints();
    void CreateElements();
    void CreateBoundaries();
    Mesh GetMesh();
}
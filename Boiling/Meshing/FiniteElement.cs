using Boiling.Meshing.Geometry;

namespace Boiling.Meshing;

public class FiniteElement(int[] nodes, int areaNumber = 0)
{
    public IList<int> Nodes { get; } = nodes;
    public int AreaNumber { get; set; } = areaNumber;
    public IList<Edge> Edges { get; } = new Edge[]
    {
        new(nodes[0], nodes[1]),
        new(nodes[0], nodes[2]),
        new(nodes[1], nodes[3]),
        new(nodes[2], nodes[3])
    };
}
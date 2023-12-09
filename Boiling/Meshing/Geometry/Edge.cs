namespace Boiling.Meshing.Geometry;

public struct Edge
{
    public int Node1 { get; set; }
    public int Node2 { get; set; }

    public Edge(int node1, int node2)
    {
        Node1 = node1;
        Node2 = node2;
    }

    public bool Equal(Edge other)
        => Node1 == other.Node1 && Node2 == other.Node2 ||
           Node1 == other.Node2 && Node2 == other.Node1;
}
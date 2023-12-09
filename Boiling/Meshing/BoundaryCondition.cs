using Boiling.Meshing.Geometry;

namespace Boiling.Meshing;

public struct Dirichlet(int node, Func<double, double, double> value)
{
    public int Node { get; set; } = node;
    public Func<double, double, double> Value { get; set; } = value;

    public void Deconstruct(out int nodeIndex, out Func<double, double, double> func)
    {
        nodeIndex = Node;
        func = Value;
    }
}

public struct Neumann(Edge border, Func<double, double, double> theta)
{
    public Edge Border { get; set; } = border;
    public Func<double, double, double> Theta { get; set; } = theta;
}

public struct Newton(Edge border, Func<double, double, double> ubeta, double beta)
{
    public Edge Border { get; set; } = border;
    public Func<double, double, double> Ubeta { get; set; } = ubeta;
    public double Beta { get; set; } = beta;
}
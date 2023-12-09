using Boiling.Meshing.Geometry;

namespace Boiling.Meshing;

public class Mesh(
    IEnumerable<Point> points,
    IEnumerable<FiniteElement> elements,
    IEnumerable<AreaProperty> materials,
    IEnumerable<Dirichlet> dirichlet,
    IEnumerable<Neumann>? neumann = null,
    IEnumerable<Newton>? newton = null)
{
    public Point[] Points { get; } = points.ToArray();
    public FiniteElement[] Elements { get; } = elements.ToArray();
    public AreaProperty[] Materials { get; } = materials.ToArray();
    public Dirichlet[] Dirichlet { get; } = dirichlet.ToArray();
    public Neumann[]? Neumann { get; } = neumann?.ToArray();
    public Newton[]? Newton { get; } = newton?.ToArray();
}
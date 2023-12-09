namespace Boiling.Meshing;

public struct AreaProperty(double lambda, double gamma, 
    Func<double, double, double> source)
{
    public double Lambda { get; set; } = lambda;
    public double Gamma { get; set; } = gamma;
    public Func<double, double, double> Source { get; set; } = source;
}
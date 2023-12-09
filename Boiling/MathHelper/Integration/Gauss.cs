using Boiling.Meshing.Geometry;

namespace Boiling.MathHelper.Integration;

public class Integration(IEnumerable<QuadratureNode<double>> quadratures)
{
    public double Integrate1D(Func<double, double> f, Interval primitives)
    {
        double a = primitives.LeftBorder;
        double b = primitives.RightBorder;
        double h = primitives.Length;

        double sum = 0.0;

        foreach (var quad in quadratures)
        {
            double qi = quad.Weight;
            double pi = (a + b + quad.Node * h) / 2.0;

            sum += qi * f(pi);
        }

        return sum * h / 2.0;
    }

    public double Integrate2D(Func<double, double, double> f, Rectangle rectangle)
    {
        var leftBottom = rectangle.LeftBottom;
        var rightTop = rectangle.RightTop;

        double hx = rightTop.X - leftBottom.X;
        double hy = rightTop.Y - leftBottom.Y;

        double sum = 0.0;

        foreach (var iquad in quadratures)
        {
            double qi = iquad.Weight;
            double pi = (leftBottom.X + rightTop.X + iquad.Node * hx) / 2.0;

            foreach (var jquad in quadratures)
            {
                double qj = jquad.Weight;
                double pj = (leftBottom.Y + rightTop.Y + jquad.Node * hy) / 2.0;

                sum += qi * qj * f(pi, pj);
            }
        }

        return sum * hx * hy / 4.0;
    }
}
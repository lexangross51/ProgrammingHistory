using System.Globalization;
using Boiling.Meshing.Geometry;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace Boiling.Meshing.JsonConverter;

public class MeshParametersJsonConverter : Newtonsoft.Json.JsonConverter
{
    public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer) { }

    public override object? ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer)
    {
        if (reader.TokenType is JsonToken.Null) return null;

        var tokens = JObject.Load(reader);

        if (tokens is null) return null;

        // Abscissa points count
        if (!int.TryParse(tokens["AbscissaPointsCount"].Value<string>(), out int abscissaPointsCount))
            return null;
        
        // Ordinate points count
        if (!int.TryParse(tokens["OrdinatePointsCount"].Value<string>(), out int ordinatePointsCount))
            return null;
        
        // Control points
        var jValue = tokens["ControlPoints"];
        var controlPoints = new Point[jValue.Count()];
        int index = 0;

        foreach (var pointJValue in jValue)
        {
            if (!Point.TryParse(pointJValue.Value<string>(), out var point))
                return null;

            controlPoints[index++] = point!.Value;
        }
        
        // Areas
        index = 0;
        jValue = tokens["Areas"];
        var areas = new Area[jValue.Count()];

        foreach (var areaJValue in jValue)
        {
            if (!int.TryParse(areaJValue["ParameterNumber"].Value<string>(), out var paramNumber))
                return null;

            var borderIndices = JsonConvert.DeserializeObject<int[]?>(areaJValue["BordersIndices"].ToString());

            if (borderIndices is null or { Length: < 1 } or { Length: > 4 }) return null;

            areas[index].ParameterNumber = paramNumber;
            areas[index].LeftBorderNumber = borderIndices[0];
            areas[index].RightBorderNumber = borderIndices[1];
            areas[index].BottomBorderNumber = borderIndices[2];
            areas[index].TopBorderNumber = borderIndices[3];
            index++;
        }
        
        // Borders
        index = 0;
        jValue = tokens["Borders"];
        var borders = new Border[jValue.Count()];

        foreach (var borderJValue in jValue)
        {
            borders[index++] = JsonConvert.DeserializeObject<Border>(borderJValue.ToString() ?? "");
        }
        
        // Abscissa splits
        var abscissaSplits = JsonConvert.DeserializeObject<int[]>(tokens["AbscissaSplits"].ToString() ?? "");
        if (abscissaSplits is null) return null;
        
        // Ordinate splits
        var ordinateSplits = JsonConvert.DeserializeObject<int[]>(tokens["OrdinateSplits"].ToString() ?? "");
        if (ordinateSplits is null) return null;
        
        // Abscissa coefficients
        var abscissaK = JsonConvert.DeserializeObject<double[]>(tokens["AbscissaK"].ToString() ?? "");
        if (abscissaK is null) return null;
        
        // Ordinate coefficients
        var ordinateK = JsonConvert.DeserializeObject<double[]>(tokens["OrdinateK"].ToString() ?? "");
        if (ordinateK is null) return null;
        
        // Mesh refinement
        if (!int.TryParse(tokens["Refinement"].ToString() ?? "", out int refinement))
            return null;
        
        // Boundary formulas
        index = 0;
        int betaIndex = 0;
        jValue = tokens["BoundaryFormulas"];
        var boundaryFormulas = new Func<double, double, double>[jValue.Count()];
        var betas = new double[2];

        foreach (var formula in jValue)
        {
            string expression = formula.Value<string>().Replace(" ", "");

            if (expression.Contains("beta"))
            {
                var commaIndex = expression.IndexOf(';');
                var funcExpression = expression.Substring(0, commaIndex).Replace("Ubeta(x,y)=", "");
                var betaExpression = expression.Substring(commaIndex + 1).Replace("beta=", "");
                expression = funcExpression;
                
                betas[betaIndex++] = double.Parse(betaExpression, CultureInfo.InvariantCulture);
            }
            else
            {
                expression = expression.Replace("f(x,y)=", "");
            }
            
            boundaryFormulas[index++] = Utilities.GetFunctionFromString(expression);   
        }
        
        // Area properties
        index = 0;
        jValue = tokens["AreaProperties"];
        var areaProperties = new AreaProperty[jValue.Count()];

        foreach (var propertyJValue in jValue)
        {
            if (!double.TryParse(propertyJValue["Lambda"].ToString() ?? "", out var lambda))
                return null;
            
            if (!double.TryParse(propertyJValue["Gamma"].ToString() ?? "", out var gamma))
                return null;

            var expression = propertyJValue["F"].Value<string>()
                .Replace(" ", "")
                .Replace("f(x,y)=", "");

            areaProperties[index].Lambda = lambda;
            areaProperties[index].Gamma = gamma;
            areaProperties[index].Source = Utilities.GetFunctionFromString(expression);
            index++;
        }

        return new MeshParameters
        {
            AbscissaPointsCount = abscissaPointsCount,
            OrdinatePointsCount = ordinatePointsCount,
            ControlPoints = controlPoints,
            Areas = areas,
            Borders = borders,
            BoundaryFormulas = boundaryFormulas,
            AreaProperties = areaProperties,
            AbscissaSplits = abscissaSplits,
            OrdinateSplits = ordinateSplits,
            AbscissaK = abscissaK,
            OrdinateK = ordinateK,
            Refinement = refinement
        };
    }

    public override bool CanConvert(Type objectType)
        => objectType == typeof(MeshParameters);
}
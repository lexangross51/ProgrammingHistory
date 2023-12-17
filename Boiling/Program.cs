using Boiling.FemContext;
using Boiling.FemContext.BasisInfo;
using Boiling.Meshing;

var meshParameters = MeshParameters.ReadJson("Input/Area.json");
var meshManager = new MeshManager(new MeshBuilder(meshParameters));
var mesh = meshManager.CreateMesh();
// Utilities.SaveMesh(mesh, @"C:\Users\lexan\source\repos\Python\meshes");
// Utilities.ElementsPerAreas(mesh, @"C:\Users\lexan\source\repos\Python\meshes");

var femSolver = new FemSolver(mesh, new BiLinearBasis(), new TimeMesh(0, 1200, 2400));
femSolver.Solve();
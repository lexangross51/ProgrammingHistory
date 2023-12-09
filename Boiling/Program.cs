using Boiling;
using Boiling.Meshing;

var meshParameters = MeshParameters.ReadJson("Input/Area.json");
var meshManager = new MeshManager(new MeshBuilder(meshParameters));
var mesh = meshManager.CreateMesh();
Utilities.SaveMesh(mesh, @"C:\Users\lexan\source\repos\Python\meshes");

// var femSolver = new FemSolver(mesh, new BiQuadraticBasis());
// femSolver.Solve();
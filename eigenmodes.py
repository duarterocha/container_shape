from pyoomph import *
from pyoomph.meshes.simplemeshes import CircularMesh
from pyoomph.output.plotting import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *


class DiffusionEquation(Equations):
    def __init__(self, g="g"):
        super(DiffusionEquation, self).__init__()
        self.g = g

    def define_fields(self):
        self.define_scalar_field(self.g, "C2")

    def define_residuals(self):
        g, gtest = var_and_test(self.g)
        self.add_residual(weak(partial_t(g), gtest) + weak(grad(g), grad(gtest)))


class Plotter(MatplotlibPlotter):

    def __init__(self, problem, circular_mesh=True, aspect_ratio=1, filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
        super(Plotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)
        self.circular_mesh = circular_mesh
        self.aspect_ratio

    def define_plot(self):
        self.set_view(xmin=-1, xmax=1, ymax=1, ymin=-1)
        cb_g = self.add_colorbar("g", position="lower left")
        self.add_plot("domain/g", colorbar=cb_g)


class LaplaceEigenProblem(Problem):

    def __init__(self, circular_mesh=False, aspect_ratio=1, name = "g"):
        super(LaplaceEigenProblem, self).__init__()
        self.circular_mesh = circular_mesh
        self.aspect_ratio = aspect_ratio
        self.name = name

    def define_problem(self):
        if self.circular_mesh:
            mesh = CircularMesh(radius=1, domain_name="domain", outer_interface="outer")
        else:
            mesh = RectangularQuadMesh(size=[self.aspect_ratio, 1], lower_left=[0, 0], N=10)
        self.add_mesh_template(mesh)

        g=self.name
        eqs = DiffusionEquation(g=g)
        eqs += MeshFileOutput()
        if self.circular_mesh:
            eqs += DirichletBC(g=0) @ "outer"
        else:
            for boundary in ["top", "bottom", "left", "right"]:
                eqs += DirichletBC(g=0) @ boundary
        eqs += RefineToLevel(level=3)

        self.plotter = []
        for n in range(10):
            self.plotter.append(Plotter(self, "real_" + str(n) + "_{:05d}", eigenmode="real", eigenvector=n))
            self.plotter.append(Plotter(self, "imag_" + str(n) + "_{:05d}", eigenmode="imag", eigenvector=n))

        self.add_equations(eqs @ "domain")

if __name__ == "__main__":
    with LaplaceEigenProblem() as problem:
        problem.circular_mesh = False
        problem.name = 'g'
        problem.set_c_compiler("tcc")
        problem.solve()
        problem.solve_eigenproblem(10)
        problem.output()

from pyoomph import *
from pyoomph.meshes.simplemeshes import CircularMesh
from pyoomph.output.plotting import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *

class Plotter(MatplotlibPlotter):

    def __init__(self, problem, filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
        super(Plotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)

    def define_plot(self):
        p = self.get_problem()
        self.background_color = "darkgrey"

        # view
        self.set_view(-problem.aspect_ratio - 0.1, -0.2, problem.aspect_ratio + 0.1, 1.1)

        # colorbars
        cb_v = self.add_colorbar("velocity", cmap="seismic", position="bottom right", factor=1)
        cb_T = self.add_colorbar("temperature", cmap="coolwarm", position="bottom left")

        # plots
        self.add_plot("domain/velocity", colorbar=cb_v)
        self.add_plot("domain/temperature", colorbar=cb_T, transform="mirror_x")
        self.add_plot("domain/velocity", mode="arrows", linecolor="black")


class ContainerShapeProblem(Problem):

    def __init__(self, aspect_ratio=1):
        super(ContainerShapeProblem, self).__init__()
        self.aspect_ratio = aspect_ratio
        self.Pr = 1
        self.Ra = self.get_global_parameter("Ra")
        self.thermal_coefficient = 1
        self.plotter = Plotter(self)

    def set_Ra(self, value):
        self.Ra.value = value

    def get_Ra(self, symbolic=False):
        if symbolic:
            return self.Ra.get_symbol()  # Symbolic for expressions
        else:
            return self.Ra.value  # current value

    def define_problem(self):
        mesh = RectangularQuadMesh(size=[self.aspect_ratio, 1], lower_left=[0, 0], N=10)
        self.add_mesh_template(mesh)
        equations = NavierStokesEquations(
            mass_density=self.thermal_coefficient / (self.Pr * self.get_Ra(symbolic=True)),
            bulkforce=(1 - self.thermal_coefficient * var('temperature')) * vector(0, -1),
            dynamic_viscosity=self.thermal_coefficient / self.get_Ra(symbolic=True))
        equations += AdvectionDiffusionEquations(fieldnames='temperature')
        for boundary in ['top', 'bottom']:
            equations += DirichletBC(velocity_x=0, velocity_y=0) @ boundary
        for boundary in ['left', 'right']:
            equations += NeumannBC(temperature=0) @ boundary
            equations += DirichletBC(velocity_x=0, velocity_y=0) @ boundary
        equations += DirichletBC(temperature=0) @ 'top'
        equations += DirichletBC(temperature=1) @ 'bottom'
        equations += InitialCondition(temperature=1 - var('coordinate_y'))
        equations += InitialCondition(velocity_x=0, velocity_y=0)
        equations += SpatialErrorEstimator(velocity=1, temperature=1)
        equations += DirichletBC(pressure=0) @ "bottom/left"
        equations += MeshFileOutput()
        self.add_equations(equations @ "domain")


if __name__ == "__main__":
    with ContainerShapeProblem() as problem:
        problem.set_c_compiler("tcc")
        problem.max_refinement_level = 1
        problem.aspect_ratio = 2
        problem.set_Ra(1200)
        problem.plotter = [problem.plotter]
        problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="real", filetrunk="eigenreal_{:05d}"))
        problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="abs", filetrunk="eigenabs_{:05d}"))
        for p in problem.plotter:
            p.file_ext = ["png"]
        # problem.run(5, 0.1, spatial_adapt=1)
        problem.solve()
        problem.output()

        # Eigenproblem
        for parameter in problem.find_bifurcation_via_eigenvalues("Ra", 100, epsilon=1e-3):
            print(parameter)
        problem.set_Ra(problem.get_Ra() + 100)
        problem.solve()
        problem.solve_eigenproblem(6)
        problem.perturb_dofs(numpy.real(problem.get_last_eigenvectors()[0]))
        problem.run(5, startstep=0.001, spatial_adapt=1)
        problem.output_at_increased_time()

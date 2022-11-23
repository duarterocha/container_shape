from pyoomph import *
from pyoomph.expressions import *
from pyoomph.meshes.simplemeshes import CircularMesh
from pyoomph.output.plotting import MatplotlibPlotter


class DiffusionEquation(Equations):
    def __init__(self, name='u', scalar_field=True):
        super(DiffusionEquation, self).__init__()
        self.name = name
        self.scalar_field = scalar_field

    def define_fields(self):
        if self.scalar_field:
            self.define_scalar_field(self.name, "C2")
        else:
            self.define_vector_field(self.name, "C2")

    def define_residuals(self):
        u, utest = var_and_test(self.name)
        if self.scalar_field:
            self.add_residual(weak(partial_t(u), utest) + weak(grad(u), grad(utest)) + weak((u - 1), utest))
        else:
            self.add_residual(weak(partial_t(u), utest) + weak(grad(u), grad(utest)) + weak((u - vector(0.5,0.5)), utest))


class Plotter(MatplotlibPlotter):

    def __init__(self, problem, name='c', filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
        super(Plotter, self).__init__(problem, filetrunk=filetrunk, fileext=fileext, eigenvector=eigenvector,
                                      eigenmode=eigenmode)
        self.name = name

    def define_plot(self):
        self.set_view(xmin=-1, xmax=1, ymax=1, ymin=-1)
        cb_u = self.add_colorbar(self.name, position="lower left")
        self.add_plot("domain/" + self.name, colorbar=cb_u)


class AxisymmetryBreaking(Problem):
    def __init__(self, size=1, one_dim=True, scalar_field=True, Neigen=6):
        super(AxisymmetryBreaking, self).__init__()
        self.one_dim = one_dim
        self.size = size
        self.scalar_field = scalar_field
        self.Neigen = Neigen
        self.param_m = self.get_global_parameter("m")

    def mesh_2d(self):
        mesh = CircularMesh(radius=self.size, domain_name="domain", outer_interface="outer")
        return mesh

    def mesh_1d(self):
        mesh = LineMesh(size=self.size, name="domain", left_name="axis", right_name="outer")
        return mesh

    def define_problem_2d(self):
        mesh = self.mesh_2d()
        self.add_mesh_template(mesh)

        eqs = DiffusionEquation(name='c', scalar_field=self.scalar_field)
        eqs += MeshFileOutput()
        if self.scalar_field:
            eqs += DirichletBC(c=0) @ "outer"
        else:
            eqs += DirichletBC(c_x=0, c_y=0)
        eqs += RefineToLevel()
        for n in range(self.Neigen):
            eqs += TextFileOutput(eigenvector=n, eigenmode="real", filename="real_" + str(n))

        self.plotter = [Plotter(self, name='c')]
        for n in range(self.Neigen):
            self.plotter.append(Plotter(self, filetrunk="real_" + str(n) + "_{:05d}", eigenmode="real", eigenvector=n))
            self.plotter.append(Plotter(self, filetrunk="imag_" + str(n) + "_{:05d}", eigenmode="imag", eigenvector=n))

        self.add_equations(eqs @ "domain")

    def define_problem_1d(self):
        mesh = self.mesh_1d()
        self.add_mesh_template(mesh)
        self.set_coordinate_system(axisymmetric)

        eqs = DiffusionEquation(name='c', scalar_field=self.scalar_field)

        m_sym = self.param_m.get_symbol()
        r = var("coordinate_x")
        c, ctest = var_and_test("c")

        eqs += WeakContribution(m_sym ** 2 / r ** 2 * c, ctest)
        eqs += TextFileOutput()
        if self.scalar_field:
            eqs += DirichletBC(c=0) @ "outer"
        else:
            eqs += DirichletBC(c_x=0, c_y=0)
        eqs += RefineToLevel()

        for n in range(self.Neigen):
            eqs += TextFileOutput(eigenvector=n, eigenmode="real", filename="real_" + str(n))

        self.add_equations(eqs @ "domain")

    def define_problem(self):
        if not self.one_dim:
            self.define_problem_2d()
        else:
            self.define_problem_1d()


with AxisymmetryBreaking() as problem:
    problem.one_dim = False
    problem.scalar_field = False
    problem.set_c_compiler("tcc")
    problem.param_m.value = 2
    problem.max_refinement_level = 2
    problem.solve()
    problem.solve_eigenproblem(problem.Neigen)
    problem.output()

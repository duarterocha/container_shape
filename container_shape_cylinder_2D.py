from pyoomph import *
from pyoomph.meshes.simplemeshes import CylinderMesh
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *
from plotter import Plotter
import math
import pyoomph.solvers.petsc

# Since curved entities do not work yet, we manually have to map the nodes on the outer mantle on the radius
class MapOnCylinderMantle(Equations):
    def __init__(self,nondim_radius):
        super(MapOnCylinderMantle, self).__init__()
        self.nondim_radius=nondim_radius

    def after_mapping_on_macro_elements(self):
        import math
        for n in self.get_mesh().nodes():
            x,y=n.x(0),n.x(1)
            r=math.sqrt(x**2+y**2)
            xn,yn=x/r,y/r
            n.set_x(0,self.nondim_radius*xn)
            n.set_x(1,self.nondim_radius*yn)

class ContainerShapeProblem(Problem):

    def __init__(self):
        super(ContainerShapeProblem, self).__init__()
        self.aspect_ratio = self.get_global_parameter("GammaF")
        self.Ra = self.get_global_parameter("Ra")
        self.Pr = 1
        self.thermal_coefficient = 1
        # self.plotter = Plotter(self)

    def set_Ra(self, value):
        self.Ra.value = value

    def get_Ra(self, symbolic=False):
        if symbolic:
            return self.Ra.get_symbol()  # Symbolic for expressions
        else:
            return self.Ra.value  # current value

    def define_problem(self):
        mesh = RectangularQuadMesh(size=[self.aspect_ratio.value, 1], lower_left=[0, 0],
                                   N=[math.ceil(self.aspect_ratio.value * 10), 10])
        self.add_mesh_template(mesh)
        equations = NavierStokesEquations(
            mass_density=self.thermal_coefficient / (self.Pr * self.get_Ra(symbolic=True)),
            bulkforce=(1 - self.thermal_coefficient * var('temperature')) * vector(0, -1),
            dynamic_viscosity=self.thermal_coefficient / self.get_Ra(symbolic=True))
        equations += AdvectionDiffusionEquations(fieldnames='temperature')
        for boundary in ['top', 'bottom', 'left', 'right']:
            equations += DirichletBC(velocity_x=0, velocity_y=0) @ boundary
        for boundary in ['left', 'right']:
            equations += NeumannBC(temperature=0) @ boundary
        equations += DirichletBC(temperature=0) @ 'top'
        equations += DirichletBC(temperature=1) @ 'bottom'
        equations += DirichletBC(pressure=0) @ "bottom/left"
        equations += InitialCondition(temperature=1 - var('coordinate_y'))
        equations += InitialCondition(velocity_x=0, velocity_y=0)
        equations += InitialCondition(pressure=- 0.5 * var('coordinate_y') ** 2)
        equations += MeshFileOutput()
        # equations += Scaling(coordinate_x=self.aspect_ratio.get_symbol())
        self.add_equations(equations @ "domain")


if __name__ == "__main__":
    with ContainerShapeProblem() as problem:
        problem.max_refinement_level = 2
        #problem.set_eigensolver("slepc")
        problem.aspect_ratio.value = 1
        problem.set_Ra(1600)

        # Eigenproblem
        for parameter in problem.find_bifurcation_via_eigenvalues("Ra", 200, epsilon=1e-3, do_solve=False, max_ds=5000, neigen=3):
            print(parameter)
        write_outfile(problem)
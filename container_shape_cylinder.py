from pyoomph import *
from pyoomph.meshes.simplemeshes import CylinderMesh
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *
from plotter import Plotter
import math

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
        #self.plotter = Plotter(self)

    def set_Ra(self, value):
        self.Ra.value = value

    def get_Ra(self, symbolic=False):
        if symbolic:
            return self.Ra.get_symbol()  # Symbolic for expressions
        else:
            return self.Ra.value  # current value

    def define_problem(self):
        mesh = CylinderMesh(radius=1, height=1/self.aspect_ratio.value, nsegments_h=4*math.ceil(1/self.aspect_ratio.value))
        self.add_mesh_template(mesh)
        equations = MeshFileOutput()
        equations += RefineToLevel(self.max_refinement_level)
        equations += MapOnCylinderMantle(1) @ "mantle"
        equations += MeshFileOutput() @ "mantle"
        equations += MeshFileOutput() @ "bottom"
        equations += MeshFileOutput() @ "top"
        equations += NavierStokesEquations(
            mass_density=self.thermal_coefficient / (self.Pr * self.get_Ra(symbolic=True)),
            bulkforce=(1 - self.thermal_coefficient * var('temperature')) * vector(0, 0, -1),
            dynamic_viscosity=self.thermal_coefficient / self.get_Ra(symbolic=True)).with_pressure_fixation()
        equations += AdvectionDiffusionEquations(fieldnames='temperature')
        for boundary in ['top', 'bottom', 'mantle']:
            equations += DirichletBC(velocity_x=0, velocity_y=0, velocity_z=0) @ boundary
        equations += NeumannBC(temperature=0) @ 'mantle'
        equations += DirichletBC(temperature=0) @ 'top'
        equations += DirichletBC(temperature=1) @ 'bottom'
        #equations += DirichletBC(pressure=0) @ "bottom/mantle"
        equations += InitialCondition(temperature=1 - var('coordinate_z'))
        equations += InitialCondition(velocity_x=0, velocity_y=0, velocity_z=0)
        equations += InitialCondition(pressure=- 0.5 * var('coordinate_z') ** 2)
        self.add_equations(equations @ "domain")

def write_outfile(problem):
    outfile = open(problem.get_output_directory("bifurcation_file.txt"), "w")
    Gamma = problem.aspect_ratio.value
    RaC = problem.get_Ra()
    outfile.write(str(Gamma) + "\t" + str(RaC) + "\n")
    outfile.flush()


if __name__ == "__main__":
    with ContainerShapeProblem() as problem:
        problem.max_refinement_level = 2
        problem.aspect_ratio.value = 10
        problem.set_Ra(1600)

        # Eigenproblem
        for parameter in problem.find_bifurcation_via_eigenvalues("Ra", 200, epsilon=1e-3, do_solve=False, max_ds=5000, neigen=3):
            print(parameter)
        write_outfile(problem)
        '''problem.solve_eigenproblem(6)
        problem.perturb_dofs(0.001*numpy.real(problem.get_last_eigenvectors()[0]))
        problem.run(3,0.2, do_not_set_IC=True)
        problem.output_at_increased_time()'''

'''problem.plotter = [problem.plotter]
    problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="real", filetrunk="eigenreal_{:05d}"))
    problem.plotter.append(Plotter(problem, eigenvector=0, eigenmode="abs", filetrunk="eigenabs_{:05d}"))
    for p in problem.plotter:
        p.file_ext = ["png"]
    # problem.run(5, 0.1, spatial_adapt=1)
    problem.solve()
    problem.output()'''
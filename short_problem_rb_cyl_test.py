from pyoomph import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.advection_diffusion import *


class RBConvectionProblem(Problem):
    def __init__(self):
        super(RBConvectionProblem, self).__init__()
        self.Ra=self.get_global_parameter("Ra")
        self.Ra.value=100
        self.Pr=self.get_global_parameter("Pr")
        self.Pr.value=1
        self.Gamma=self.get_global_parameter("Gamma")
        self.Gamma.value=8
        self.desired_ndof=10000 # Approx number of degrees of freedom

    def define_problem(self):
        self.set_coordinate_system(axisymmetric)
        
        ndesired=self.desired_ndof/(10) # About 10 dofs per element if arrange in a larger grid
        Nx=math.ceil(math.sqrt(ndesired*(self.Gamma.value/2)))
        Ny=math.ceil(ndesired/Nx)
        mesh=RectangularQuadMesh(size=[self.Gamma.value/2, 1],N=[Nx,Ny])        
        self.add_mesh_template(mesh)
        
        Pr,Ra=self.Pr.get_symbol(), self.Ra.get_symbol()
        # we scale the -grad(p) term in the NS with Pr*Ra. This is just a rescaling of the pressure, but it helps since the initial condition for p is then independent of Ra
        eqs=NavierStokesEquations(mass_density=1,dynamic_viscosity=Pr,bulkforce=Pr*Ra*var("T")*vector(0,1),pressure_factor=Pr*Ra)
        eqs+=AdvectionDiffusionEquations(fieldnames="T",diffusivity=1,space="C1")
        eqs+=MeshFileOutput()
        eqs+=MeshFileOutput(eigenvector=0,eigenmode="real",filetrunk="eigen0real")
        eqs+=DirichletBC(T=0)@"bottom"
        eqs+=DirichletBC(T=-1)@"top"
        z=var("coordinate_y")
        eqs+=InitialCondition(T=-z)
        eqs+=InitialCondition(pressure=-z**2/2)

        for ns_wall in ["top","right","bottom"]:
            eqs+=NoSlipBC()@ns_wall
           
        # This will set velocity_x=0, but it also takes care of everything else when investigating the axial symmetry breaking
        # See below
        eqs+=AxisymmetryBC()@"left"

        self.add_equations(eqs@"domain")
        
        # Activate the axial symmetry investigation
        # it will do the following:
        #	1.	Introduce a global parameter m that holds the current angular mode (default: 0 -> axisymmetric)
        #   2. Switch the coordinate system to a variant of the AxisymmetricCoordinateSystem
        #		a) this will add a vector component velocity_theta (and likewise to all other vector fields)
        #		b)	div and grad will include theta derivatives, which expand to \partial_\theta f = i*m f        
        #	3.	Activate complex-conjugation in the weak forms, i.e. weak(a,b)=integral[a*conjugate(b)]*dx
        #			Complex terms will only arise due to \partial_\theta terms!
        #	4.	Add a modificator for all add_residual statements in the equations
        #		a) it will add real(residual) to the normal residual
        #		b)	and imag(residual) to another residual
        #	5.	Tell the eigensolver to assemble the mass matrix M and the Jacobian J by
        #			M=M_real+I*M_imag
        #			J=J_real+I*J_imag			
        #		where we both the real and imaginary residual are evaluated.
        #		Normal solving of the system (problem.solve/run) just considers the real residual for m=0. 
        #		Hence, all \partial_\theta are zero for problem.solve/run => axisymmetric!
        #	6.	All AxisymmetryBC objects will do the following:
        #		When the normal problem is solved (problem.solve/run):
        #			a)	just pin all radial vector components to zero, i.e. like DirichletBC(velocity_x=0)@"left"
        #			b)	add a DirichletBC(velocity_theta=0)@"left"	# We cannot have a theta component at r=0 for axisymmetric solves
        #			c)	likewise for all other vector fields 
        #			d)	keep all other field untouched for normal solving
        #		When an eigenproblem is solved:
        #			a) deactivate the DirichletBCs for e.g. velocity_x and velocity_theta. Hence, in general, they can be nonzero in the eigensolution
        #			b)	add matrix modificators to change the matrices M and J during eigensolving
        #				- [m=0]:	Enforce velocity_x=0 and velocity_theta=0 of the eigensolution at r=0 (likewise all other vector fields)
        #				- [|m|=1]:	Enforce velocity_y=0 (likewise all other vector fields) and all scalar fields=0 of the eigensolution at r=0  
        #				- [|m|>=2]: Enforce all vector and scalar fields of the eigensolution to be zero at r=0
        #				This enforcing is done as follows:
        #					we solve the generalized eigenproblem:
        #						lambda*M*v=J*v
        #				To set a degree of freedom to zero, we just remove the row corresponding to that degree from M (so it is a zero row)
        #				and we remove the row from J, but set a 1 on the diagonal.
        #				For this degree (call it u here), lambda*M*v=J*v will hence just give 0=u, i.e. the value in the eigensolution must be zero
        #				Thereby, all eigenfunctions have the same number of unknowns, which is required for internal reasons in pyoomph
        self.define_problem_for_axial_symmetry_breaking_investigation()
        # What we must do by hand, however, is setting the no-slip also for velocity_theta
        self.add_equations((DirichletBC(velocity_theta=0)@["bottom","right","top"])@"domain")
        
        

with RBConvectionProblem() as problem:
    problem.set_c_compiler("tcc")
    problem.set_eigensolver("pardiso")
    problem.Ra.value=1000
    problem.initialise()
    # We can pass the angular modes we want to inspect here:
    for param in problem.find_bifurcation_via_eigenvalues("Ra",100,do_solve=False,axial_m=[0,1,2,3]):
        print(param) # Eigenvalues will be sorted again by real part
        problem.output_at_increased_time()
    problem.output()

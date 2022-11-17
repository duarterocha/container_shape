from pyoomph.output.plotting import *


class Plotter(MatplotlibPlotter):

    def __init__(self, problem, filetrunk="plot_{:05d}", fileext="png", eigenvector=None, eigenmode="abs"):
        super(Plotter, self).__init__(problem, filetrunk, fileext, eigenvector, eigenmode)

    def define_plot(self):
        p = self.get_problem()
        self.background_color = "darkgrey"

        # view
        self.set_view(-p.aspect_ratio.value - 0.1, -0.2, p.aspect_ratio.value + 0.1, 1.1)

        # colorbars
        cb_v = self.add_colorbar("velocity", cmap="seismic", position="bottom right", factor=1)
        cb_T = self.add_colorbar("temperature", cmap="coolwarm", position="bottom left")

        # plots
        self.add_plot("domain/velocity", colorbar=cb_v)
        self.add_plot("domain/temperature", colorbar=cb_T, transform="mirror_x")
        self.add_plot("domain/velocity", mode="arrows", linecolor="black")
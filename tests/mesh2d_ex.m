fname = "cycle.mesh2d";

mesh = mesh2d_read(fname)

mesh2d_plot(mesh,2)

exportgraphics(gcf, "cycle.png")
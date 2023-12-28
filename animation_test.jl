using PyCall
using PyPlot
mpla = pyimport("matplotlib.animation")
np = pyimport("numpy")
Cell = pyimport("celluloid")
Camera = Cell.Camera


# Default plot params
function plot_params()
    plt.rc("xtick", labelsize=16)
    plt.rc("ytick", labelsize=16)
    plt.rc("xtick", direction="in")
    plt.rc("ytick", direction="in")
    plt.rc("font", size=15)
    plt.rc("figure", autolayout="True")
    plt.rc("axes", titlesize=16)
    plt.rc("axes", labelsize=17)
    plt.rc("xtick.major", width=1.5)
    plt.rc("xtick.major", size=5)
    plt.rc("ytick.major", width=1.5)
    plt.rc("ytick.major", size=5)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=1.5)
    plt.rc("legend", fontsize=13)
    plt.rc("mathtext", fontset="stix")
    plt.rc("font", family="STIXGeneral")
  
    # Default width for Nature is 7.2 inches, 
    # height can be anything
    #plt.rc("figure", figsize=(7.2, 4.5))
end

plot_params()
fig = PyPlot.figure(figsize=(20, 20))

ax = fig.add_subplot(2,2,1)

t = np.linspace(0, 3, 40)
g = -9.81
v0 = 12
z = g * t.^2 / 2 + v0 * t

v02 = 5
z2 = g * t.^2 / 2 + v02 * t

# generste the first figure
# scat = ax.scatter(t[1], z[1], c="b", s=5,)
# line2 = ax.plot(t[1], z2[1], )[1]
ax.set(xlim=[0, 3], ylim=[-4, 10], xlabel="Time [s]", ylabel="Z [m]")
# ax.legend()
camera = Camera(fig)

for i in eachindex(t)
    ax.scatter(t[1:i], z[1:i], c="b", s=5,)
    plt.pause(0.1)
    camera.snap()
end

animation = camera.animate()
animation.save("animation.gif", writer="pillow", fps=2)


# function update(frame)
#     # for each frame, update the data stored on each artist.
#     x = t[:frame]
#     y = z[:frame]
#     # update the scatter plot:
#     data = np.stack([x, y]).T
#     scat.set_offsets(data)
#     # update the line plot:
#     line2.set_xdata(t[:frame])
#     line2.set_ydata(z2[:frame])
#     return (scat, line2)
# end

# ani = mpla.FuncAnimation(fig=fig, func=update, frames=40, interval=30)
# figname = string("movie.mp4")
# writer = mpla.PillowWriter(fps=15, bitrate=1800)
# ani.save(figname, writer=writer, dpi = 300)

# pylint: disable=no-member
# pylint: disable=import-error
# pylint: disable=no-name-in-module
import random
import sys
import warnings
import os
import copy
import importlib
import inspect
import multiprocessing
import tqdm
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import kmc.variables
from kmc.system import System
import kmc.bimolecular
from kmc import __version__ as kmc_version
import kmc.utils

warnings.filterwarnings("ignore")

print("####################################################################")
print("Xcharge: A Kinetic Monte Carlo Model for Exciton and Charge Dynamics")
print("Version  : " + kmc_version.__version__)
print()
output_header = (
    "# Version  : "
    + kmc_version.__version__
    + "\n"
)
### end interface


# importing param module
spec = importlib.util.spec_from_file_location(
    sys.argv[1].split(".")[0], os.path.join(os.getcwd(), sys.argv[1])
)
param = importlib.util.module_from_spec(spec)
spec.loader.exec_module(param)

argumentos = []
for name, value in vars(param).items():
    if hasattr(value, "assign_to_system") and not inspect.isclass(
        value
    ):  # getting only instancied classes
        argumentos.append(value)
        # print(name,hasattr(value, 'make'))

# getting all essential info from user's input
monomolecular = param.monomolecular
processes = param.processes


def set_variables(variable_name):
    try:
        return getattr(param, variable_name)
    except AttributeError:
        return getattr(kmc.variables, variable_name)


identifier = set_variables("identifier")
animation_mode = set_variables("animation_mode")
save_animation = set_variables("save_animation")
animation_exten = set_variables("animation_exten")
time_limit = set_variables("time_limit")
pause = set_variables("pause")
marker_type = set_variables("marker_type")
rotate = set_variables("rotate")
frozen_lattice = set_variables("frozen_lattice")
bimolec = set_variables("bimolec")
periodic = set_variables("periodic")
n_proc = set_variables("n_proc")
rounds = set_variables("rounds")
cutoff = set_variables("cutoff")
#####


def passar(*args):
    pass


def anni(system, array, local):
    anni_general(system, array, local)


if bimolec:
    bi_func = anni
else:
    bi_func = passar


def regular_distance(system, local, destination=None):
    if destination is not None:
        dx = system.X[destination] - system.X[local]
        dy = system.Y[destination] - system.Y[local]
        dz = system.Z[destination] - system.Z[local]
    else:
        # dx,dy,dz = kmc.utils.distance(system.X,system.Y,system.Z,len(system.X),local)
        dx = system.X - system.X[local]
        dy = system.Y - system.Y[local]
        dz = system.Z - system.Z[local]
    return dx, dy, dz


def periodic_distance(system, local, destination=None):
    if destination is None:
        dx = system.X - system.X[local]
        dy = system.Y - system.Y[local]
        dz = system.Z - system.Z[local]
    else:
        dx = system.X[destination] - system.X[local]
        dy = system.Y[destination] - system.Y[local]
        dz = system.Z[destination] - system.Z[local]
    if system.Lx > 0:
        dx -= system.Lx * np.round(dx / system.Lx)
    if system.Ly > 0:
        dy -= system.Ly * np.round(dy / system.Ly)
    if system.Lz > 0:
        dz -= system.Lz * np.round(dz / system.Lz)
    return dx, dy, dz


if periodic:
    distance = periodic_distance
else:
    distance = regular_distance


def make_system():
    # Create instance of system
    system = System()
    # Sets system properties
    for argumento in argumentos:
        argumento.assign_to_system(system)
    system.set_basic_info(
        monomolecular,
        processes,
        identifier,
        animation_mode,
        time_limit,
        pause,
        bimolec,
        distance,
    )

    return system


# runs the annihilations defined in anni_funcs_array
def anni_general(system, anni_dict, local):
    Ss = system.particles.copy()
    locs = np.array([s.position for s in Ss])
    superpostos = np.where(locs == local)[0]
    if len(superpostos) > 1:
        try:
            anni_dict[tuple(sorted(Ss[i].species for i in superpostos[:2]))](
                Ss, system, superpostos[:2]
            )
        except KeyError:  # in case of the set (particle1,particle2) is not defined
            pass


def decision(s, system):
    kind = s.species
    local = s.position
    dx, dy, dz = distance(system, local)
    r = kmc.utils.distances(dx, dy, dz, len(dx))
    cut = np.where(r < cutoff)[0]
    r = r[cut]
    dx = dx[cut]
    dy = dy[cut]
    dz = dz[cut]
    mats = system.mats[cut]
    hop = system.processes[kind]
    mono = system.monomolecular[kind]
    jump_rate = [
        transfer.rate(
            r=r,
            dx=dx,
            dy=dy,
            dz=dz,
            system=system,
            particle=s,
            mats=mats,
            matlocal=system.mats[local],
            cut=cut,
        )
        for transfer in hop
    ]
    mono_rate = [[m.rate(material=system.mats[local])] for m in mono]
    jump_rate.extend(mono_rate)
    sizes = np.array([len(i) for i in jump_rate])
    jump_rate = np.concatenate(jump_rate)
    labels = hop + mono
    soma, jump = kmc.utils.jump(jump_rate, len(jump_rate), random.uniform(0, 1))
    destino = np.argmax(np.cumsum(sizes) - 1 >= jump)
    s.process = labels[destino]
    if destino < len(hop):
        s.destination = cut[int(jump - np.sum(sizes[:destino]))]
    else:
        s.destination = local
    return soma


########ITERATION FUNCTIONS#######################################################
def step_noanimation(system):
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        R, dests = [], []
        for s in Ss:
            if s in system.particles:
                Rs = decision(s, system)
                if s.destination not in dests:
                    s.process.action(s, system, s.destination)
                    bi_func(system, kmc.bimolecular.bimolec_funcs_array, s.destination)
                    R.append(Rs)
                    dests.append(s.destination)
        R = np.array(R)
        system.time += np.mean((1 / R) * np.log(1 / random.uniform(0, 1)))
        for s in Ss:
            s.stamp_time(system)
    Ss = system.particles.copy()
    for s in Ss:
        s.kill("alive", system, system.s1, "alive")
        s.stamp_time(system)


def step_animation(system):
    while system.count_particles() > 0 and system.time < system.time_limit:
        system.IT += 1
        Ss = system.particles.copy()
        random.shuffle(Ss)
        R, dests = [], []
        for s in Ss:
            if s in system.particles:
                Rs = decision(s, system)
                if s.destination not in dests:
                    s.process.action(s, system, s.destination)
                    bi_func(system, kmc.bimolecular.bimolec_funcs_array, s.destination)
                    R.append(Rs)
                    dests.append(s.destination)
        R = np.array(R)
        system.time += np.mean((1 / R) * np.log(1 / random.uniform(0, 1)))
        for s in Ss:
            s.stamp_time(system)
        return Ss
    Ss = system.particles.copy()
    for s in Ss:
        s.kill("alive", system, system.s1, "alive")
        s.stamp_time(system)


##########################################################################################

if animation_mode:
    step = step_animation
else:
    step = step_noanimation


def open_log():
    filename = f"Simulation_{identifier}.txt"
    if os.path.isfile(filename) == False:
        with open(filename, "w") as f:
            f.write(output_header)
            texto = "Time,DeltaX,DeltaY,DeltaZ,Type,Energy,Location,FinalX,FinalY,FinalZ,CausaMortis,Status"
            f.write(texto + "\n")
    return filename


# Prints Spectra
def spectra(system, f):
    texto = ""
    for s in system.dead:
        texto += s.write()
    f.write(texto + "END\n")


def animate(num, system, ax, marker_option, do_rotate):
    Ss = step(system)
    X, Y, Z = system.X, system.Y, system.Z
    mats = system.mats
    ax.clear()
    # ploting the sites according to mat index
    colors_dic = {0: "black", 1: "blue", 2: "red", 3: "green", 4: "yellow"}
    n_mats = np.unique(mats)
    for mat in n_mats:
        X_mat = X[mats == mat]
        Y_mat = Y[mats == mat]
        Z_mat = Z[mats == mat]
        ax.scatter(X_mat, Y_mat, Z_mat, alpha=0.25, color=colors_dic.get(int(mat)))
    try:
        for s in Ss:
            xs = X[s.position]
            ys = Y[s.position]
            zs = Z[s.position]

            if marker_option == 1:
                ax.scatter(
                    xs,
                    ys,
                    zs,
                    marker=s.marker,
                    color=s.color,
                    s=200,
                    alpha=1,
                    label=s.species,
                )
            if marker_option == 0:
                ax.scatter(xs, ys, zs, color=s.color, s=100, alpha=1, label=s.species)
    except (IndexError,TypeError):
        pass

    if do_rotate:  # rotating the animation by an angle of IT
        ax.view_init(azim=system.IT)

    # removing duplicates on legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    # sort legend by label
    by_label = dict(sorted(by_label.items(), key=lambda x: x[0]))
    plt.legend(by_label.values(), by_label.keys())
    ax.text2D(
        0.03, 0.98, f"Time: {system.time:.2e} ps", transform=ax.transAxes
    )  # time
    ax.text2D(
        0.03, 0.94,f"Particles: {len(system.particles):.0f}", transform=ax.transAxes
    )  # npart

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")

    return (ax,)

# resets particles' initial position for a given system
def reroll_system(system):
    system.reset_particles()
    for argumento in argumentos:
        class_name = argumento.__class__.__name__
        if class_name in ["CreateParticles"]:
            argumento.assign_to_system(system)
    return system


# RUN of a single round
def RUN(
    dynamic,
):  # ROUND DYNAMICS WHERE, FOR EACH INSTANCE, THE LATTICE IS RECALCULATED
    system = make_system()
    step(system)
    return system


def RUN_FREEZE(
    dynamic,
):  # ROUND DYNAMICS WHERE, FOR EACH INSTANCE, THE LATTICE REMAINS INTACT
    syst = dynamic[1]
    system = reroll_system(copy.deepcopy(syst))
    step(system)
    return system


# setting up the animation object and adding responses to events
def run_animation():
    ani_running = True

    def onClick(event):  # if somenone clicks on the ani, this happens
        nonlocal ani_running
        if ani_running:
            ani.event_source.stop()
            ani_running = False
        else:
            ani.event_source.start()
            ani_running = True

    def pause_plot(event, pause):  # if pause = true, this will happen
        nonlocal ani_running
        if pause:
            ani.event_source.stop()
            ani_running = False

    system = make_system()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    fig.canvas.mpl_connect("button_press_event", onClick)  # pausing if clicking
    fig.canvas.mpl_connect(
        "draw_event", lambda event: pause_plot(event, pause)
    )  # pausing if pause = True at the first frame

    ani = animation.FuncAnimation(
        fig,
        animate,
        fargs=[system, ax, marker_type, rotate],
        interval=25,
        blit=False,
        repeat=False,
        cache_frame_data=True,
    )  # ,save_count=1000)

    return ani


def main():
    if animation_mode:
        ani = run_animation()
        path = identifier + "_animation." + animation_exten

        if save_animation:
            # save .gif
            if animation_exten == "gif":
                ani.save(path, writer="imagemagick", fps=10)

            # save .mp4
            if animation_exten == "mp4":
                writervideo = animation.FFMpegWriter(fps=10)
                ani.save(path, writer=writervideo)

        plt.show()
    else:
        p = multiprocessing.Pool(n_proc)
        filename = open_log()
        if not frozen_lattice:  # at every round, the entire lattice is recalculated
            run = RUN
            args = [(i) for i in range(rounds)]
        else:  # at every round, only particle creation is recalculated
            syst = make_system()
            run = RUN_FREEZE
            args = [(i, syst) for i in range(rounds)]
        with open(filename, "a") as f:
            for result in tqdm.tqdm(p.imap(run, args), total=rounds):
                spectra(result, f)


if __name__ == "__main__":
    sys.exit(main())

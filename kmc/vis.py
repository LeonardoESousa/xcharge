# kmc/viewer.py

import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection


class KMCViewer:
    def __init__(
        self,
        system,
        step_function,
        fig,
        ax,
        colors_dic,
        sizes_dic,
        material_label=None,
        distance_function=None,
        marker_option=None,
        rotate=False,
        square_ratio=True,
        clean_vis=False,
        scatter_alpha=0.35,
        image_scatter_alpha=0.12,
        print_site_position=False,
        material_leg=True,
        time_units="s",
        periodic=False,
        lattice_marker_size_scale=1.0,
        image_lattice_marker_size_scale=1.0,
        particle_marker_size=200,
        particle_alpha=1.0,
        fp_linewidth=4.0,
        fp_alpha=0.75,
        fp_image_marker_size=120,
        fp_image_alpha=0.45,
        axis_padding_fraction=0.08,
        default_material_color="gray",
        default_material_size=50,
        default_particle_color="black",
        default_particle_marker="o",
    ):
        self.system = system
        self.step_function = step_function
        self.distance_function = distance_function

        self.fig = fig
        self.ax = ax

        self.colors_dic = colors_dic
        self.sizes_dic = sizes_dic
        self.material_label = material_label or {}

        self.marker_option = marker_option
        self.rotate = rotate
        self.square_ratio = square_ratio
        self.clean_vis = clean_vis
        self.scatter_alpha = scatter_alpha
        self.image_scatter_alpha = image_scatter_alpha
        self.print_site_position = print_site_position
        self.material_leg = material_leg
        self.time_units = time_units
        self.periodic = periodic

        self.lattice_marker_size_scale = lattice_marker_size_scale
        self.image_lattice_marker_size_scale = image_lattice_marker_size_scale
        self.particle_marker_size = particle_marker_size
        self.particle_alpha = particle_alpha
        self.fp_linewidth = fp_linewidth
        self.fp_alpha = fp_alpha
        self.fp_image_marker_size = fp_image_marker_size
        self.fp_image_alpha = fp_image_alpha
        self.axis_padding_fraction = axis_padding_fraction

        self.default_material_color = default_material_color
        self.default_material_size = default_material_size
        self.default_particle_color = default_particle_color
        self.default_particle_marker = default_particle_marker

        self.empty = np.asarray([], dtype=float)

        (
            self.cell_origin,
            self.cell_min,
            self.cell_max,
            self.cell_length,
        ) = self._infer_cell()

        self.material_artists = {}
        self.image_lattice_artists = {}
        self.active_image_shifts = set()
        self.axis_signature = None

        self.particle_artists = {}
        self.fp_image_artists = {}
        self.particle_style = {}

        self.empty_segment = np.asarray(
            [
                [
                    [np.nan, np.nan, np.nan],
                    [np.nan, np.nan, np.nan],
                ]
            ],
            dtype=float,
        )

        self.fp_collection = Line3DCollection(
            self.empty_segment,
            linewidths=self.fp_linewidth,
            alpha=self.fp_alpha,
        )

        self.frame_text = None
        self.npart_text = None

        self.particle_legend = None
        self.material_legend = None
        self.last_legend_species = None

        self._configure_axes()
        self._draw_static_lattice()
        self._draw_site_indices()
        self._draw_material_legend()
        self._draw_text_panel()

        self.ax.add_collection3d(self.fp_collection)

    def _infer_cell(self):
        data_min = np.asarray(
            [
                np.amin(self.system.X),
                np.amin(self.system.Y),
                np.amin(self.system.Z),
            ],
            dtype=float,
        )

        data_max = np.asarray(
            [
                np.amax(self.system.X),
                np.amax(self.system.Y),
                np.amax(self.system.Z),
            ],
            dtype=float,
        )

        length = np.asarray(
            [
                getattr(self.system, "Lx", 0.0),
                getattr(self.system, "Ly", 0.0),
                getattr(self.system, "Lz", 0.0),
            ],
            dtype=float,
        )

        origin = data_min.copy()
        cell_min = data_min.copy()
        cell_max = data_max.copy()

        if self.periodic:
            tol = 1.0e-8

            for axis in range(3):
                if length[axis] <= 0.0:
                    continue

                if (
                    data_min[axis] >= -tol
                    and data_max[axis] <= length[axis] + tol
                ):
                    origin[axis] = 0.0
                else:
                    origin[axis] = data_min[axis]

                cell_min[axis] = origin[axis]
                cell_max[axis] = origin[axis] + length[axis]

        return origin, cell_min, cell_max, length

    def _configure_axes(self):
        self._set_axes_for_shifts(set())

    def _set_axes_for_shifts(self, shifts):
        all_shifts = set(shifts)
        all_shifts.add((0.0, 0.0, 0.0))

        shift_array = np.asarray(list(all_shifts), dtype=float)

        min_shift = np.amin(shift_array, axis=0)
        max_shift = np.amax(shift_array, axis=0)

        lim_min = self.cell_min + min_shift
        lim_max = self.cell_max + max_shift

        span = lim_max - lim_min
        pad = self.axis_padding_fraction * max(
            span[0],
            span[1],
            span[2],
            1.0,
        )

        signature = tuple(np.round(np.r_[lim_min, lim_max], 8))

        if signature == self.axis_signature:
            return

        self.axis_signature = signature

        self.ax.set_xlim(lim_min[0] - pad, lim_max[0] + pad)
        self.ax.set_ylim(lim_min[1] - pad, lim_max[1] + pad)
        self.ax.set_zlim(lim_min[2] - pad, lim_max[2] + pad)

        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")

        self.ax.set_proj_type("persp")

        if self.square_ratio:
            self.ax.set_box_aspect(
                [
                    max(span[0], 1.0e-12),
                    max(span[1], 1.0e-12),
                    max(span[2], 1.0e-12),
                ]
            )

        if self.clean_vis:
            self.ax.set_axis_off()
            self.ax.grid(False)
            self.ax.set_position([0, 0, 1, 1])

    def _draw_static_lattice(self):
        X = self.system.X
        Y = self.system.Y
        Z = self.system.Z
        mats = np.asarray(self.system.mats)

        for mat in [m for m in np.unique(mats) if m != 999]:
            mat = int(mat)
            mask = mats == mat

            color = self.colors_dic.get(mat, self.default_material_color)
            size = self.sizes_dic.get(mat, self.default_material_size)
            size *= self.lattice_marker_size_scale

            artist = self.ax.scatter(
                X[mask],
                Y[mask],
                Z[mask],
                marker="o",
                color=color,
                s=size,
                alpha=self.scatter_alpha,
                depthshade=True,
                edgecolors="none",
            )

            self.material_artists[mat] = artist

    def _draw_image_lattice(self, shift):
        shift = self._canonical_shift(shift)

        if shift == (0.0, 0.0, 0.0):
            return

        if shift in self.image_lattice_artists:
            return

        X = self.system.X + shift[0]
        Y = self.system.Y + shift[1]
        Z = self.system.Z + shift[2]
        mats = np.asarray(self.system.mats)

        artists = []

        for mat in [m for m in np.unique(mats) if m != 999]:
            mat = int(mat)
            mask = mats == mat

            color = self.colors_dic.get(mat, self.default_material_color)
            size = self.sizes_dic.get(mat, self.default_material_size)
            size *= self.lattice_marker_size_scale
            size *= self.image_lattice_marker_size_scale

            artist = self.ax.scatter(
                X[mask],
                Y[mask],
                Z[mask],
                marker="o",
                color=color,
                s=size,
                alpha=self.image_scatter_alpha,
                depthshade=True,
                edgecolors="none",
            )

            artists.append(artist)

        self.image_lattice_artists[shift] = artists

    def _update_periodic_images(self, required_shifts):
        new_shift = False

        for shift in required_shifts:
            shift = self._canonical_shift(shift)

            if shift == (0.0, 0.0, 0.0):
                continue

            if shift not in self.active_image_shifts:
                self.active_image_shifts.add(shift)
                self._draw_image_lattice(shift)
                new_shift = True

        if new_shift:
            self._set_axes_for_shifts(self.active_image_shifts)

    def _draw_site_indices(self):
        if not self.print_site_position:
            return

        X = self.system.X
        Y = self.system.Y
        Z = self.system.Z

        for i in range(len(X)):
            self.ax.text(
                X[i],
                Y[i],
                Z[i],
                str(i),
                fontsize=7,
            )

    def _draw_text_panel(self):
        self.frame_text = self.ax.text2D(
            0.03,
            0.97,
            "",
            transform=self.ax.transAxes,
        )

        self.npart_text = self.ax.text2D(
            0.03,
            0.93,
            "",
            transform=self.ax.transAxes,
        )

    def _draw_material_legend(self):
        if not self.material_leg:
            return

        if len(self.material_label) == 0:
            return

        legend_elements = []

        for name, mat in self.material_label.items():
            mat = int(mat)
            color = self.colors_dic.get(mat, self.default_material_color)

            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="None",
                    label=name,
                    markerfacecolor=color,
                    markeredgecolor=color,
                    markersize=10,
                )
            )

        self.material_legend = self.ax.legend(
            handles=legend_elements,
            title="Materials",
            loc="lower right",
        )

        self.ax.add_artist(self.material_legend)

    def _hashable_color(self, color):
        try:
            hash(color)
            return color
        except TypeError:
            return tuple(np.asarray(color).reshape(-1))

    def _particle_marker(self, particle):
        return getattr(
            particle,
            "marker",
            self.default_particle_marker,
        )

    def _particle_color(self, particle):
        return getattr(
            particle,
            "color",
            self.default_particle_color,
        )

    def _particle_key(self, particle):
        color = self._particle_color(particle)

        return (
            particle.species,
            self._particle_marker(particle),
            self._hashable_color(color),
        )

    def _ensure_particle_artist(self, particle):
        key = self._particle_key(particle)

        if key in self.particle_artists:
            return

        species = particle.species
        marker = self._particle_marker(particle)
        color = self._particle_color(particle)

        artist = self.ax.scatter(
            self.empty,
            self.empty,
            self.empty,
            marker=marker,
            color=color,
            s=self.particle_marker_size,
            alpha=self.particle_alpha,
            depthshade=True,
            label=species,
        )

        self.particle_artists[key] = artist

        self.particle_style[species] = {
            "marker": marker,
            "color": color,
        }

    def _ensure_fp_image_artist(self, particle):
        key = self._particle_key(particle)

        if key in self.fp_image_artists:
            return

        marker = self._particle_marker(particle)
        color = self._particle_color(particle)

        artist = self.ax.scatter(
            self.empty,
            self.empty,
            self.empty,
            marker=marker,
            color=color,
            s=self.fp_image_marker_size,
            alpha=self.fp_image_alpha,
            depthshade=True,
        )

        self.fp_image_artists[key] = artist

    def _set_scatter_xyz(self, artist, coords):
        if len(coords) == 0:
            artist._offsets3d = (
                self.empty,
                self.empty,
                self.empty,
            )
            return

        coords = np.asarray(coords, dtype=float)

        artist._offsets3d = (
            coords[:, 0],
            coords[:, 1],
            coords[:, 2],
        )

    def _update_particle_markers(self, particles):
        X = self.system.X
        Y = self.system.Y
        Z = self.system.Z

        alive = [
            p
            for p in particles
            if p.status == "alive"
        ]

        by_key = {}
        current_species = set()

        for particle in alive:
            self._ensure_particle_artist(particle)

            key = self._particle_key(particle)

            by_key.setdefault(key, []).append(
                [
                    X[particle.position],
                    Y[particle.position],
                    Z[particle.position],
                ]
            )

            current_species.add(particle.species)

        for key, artist in self.particle_artists.items():
            self._set_scatter_xyz(
                artist,
                by_key.get(key, []),
            )

        return current_species

    def _as_scalar(self, value):
        return float(np.asarray(value).reshape(-1)[0])

    def _is_pair_particle(self, particle):
        species = particle.species.lower()

        return (
            "frenkelpair" in species
            or "i2" in species
        )

    def _raw_position(self, site):
        return np.asarray(
            [
                self.system.X[site],
                self.system.Y[site],
                self.system.Z[site],
            ],
            dtype=float,
        )

    def _distance_displacement(self, local, destination):
        if self.distance_function is not None:
            dx, dy, dz = self.distance_function(
                self.system,
                local,
                destination,
            )

            return np.asarray(
                [
                    self._as_scalar(dx),
                    self._as_scalar(dy),
                    self._as_scalar(dz),
                ],
                dtype=float,
            )

        dx = self.system.X[destination] - self.system.X[local]
        dy = self.system.Y[destination] - self.system.Y[local]
        dz = self.system.Z[destination] - self.system.Z[local]

        if self.periodic:
            if self.system.Lx > 0:
                dx -= self.system.Lx * np.floor(
                    dx / self.system.Lx + 0.5
                )

            if self.system.Ly > 0:
                dy -= self.system.Ly * np.floor(
                    dy / self.system.Ly + 0.5
                )

            if self.system.Lz > 0:
                dz -= self.system.Lz * np.floor(
                    dz / self.system.Lz + 0.5
                )

        return np.asarray(
            [
                dx,
                dy,
                dz,
            ],
            dtype=float,
        )

    def _canonical_shift(self, shift):
        shift = np.asarray(shift, dtype=float).copy()

        for axis, length in enumerate(self.cell_length):
            if length <= 0.0:
                shift[axis] = 0.0
                continue

            n = int(np.round(shift[axis] / length))
            shift[axis] = n * length

            if abs(shift[axis]) < 1.0e-12:
                shift[axis] = 0.0

        return tuple(float(x) for x in shift)

    def _nearest_image_endpoint_and_shift(self, particle):
        local = particle.position
        destination = particle.ghost_site

        r0 = self._raw_position(local)
        r_destination = self._raw_position(destination)

        dr = self._distance_displacement(
            local,
            destination,
        )

        r1 = r0 + dr
        shift = self._canonical_shift(r1 - r_destination)
        r1 = r_destination + np.asarray(shift, dtype=float)

        return r0, r1, shift

    def _update_fp_lines(self, particles):
        segments = []
        colors = []
        image_coords_by_key = {}
        required_shifts = set()

        for particle in particles:
            if particle.status != "alive":
                continue

            if not self._is_pair_particle(particle):
                continue

            if getattr(particle, "ghost_site", None) is None:
                continue

            self._ensure_fp_image_artist(particle)

            r0, r1, shift = self._nearest_image_endpoint_and_shift(particle)
            color = self._particle_color(particle)
            key = self._particle_key(particle)

            segments.append(
                np.asarray(
                    [
                        r0,
                        r1,
                    ],
                    dtype=float,
                )
            )

            colors.append(color)

            if shift != (0.0, 0.0, 0.0):
                required_shifts.add(shift)

                image_coords_by_key.setdefault(key, []).append(
                    [
                        r1[0],
                        r1[1],
                        r1[2],
                    ]
                )

        self._update_periodic_images(required_shifts)

        self.fp_collection.set_segments(segments)

        if len(colors) > 0:
            self.fp_collection.set_colors(colors)
        else:
            self.fp_collection.set_colors([])

        for key, artist in self.fp_image_artists.items():
            self._set_scatter_xyz(
                artist,
                image_coords_by_key.get(key, []),
            )

    def _update_particle_legend(self, current_species):
        if current_species == self.last_legend_species:
            return

        self.last_legend_species = set(current_species)

        if self.particle_legend is not None:
            self.particle_legend.remove()
            self.particle_legend = None

        if len(current_species) == 0:
            if self.material_legend is not None:
                self.ax.add_artist(self.material_legend)
            return

        handles = []

        for species in sorted(current_species):
            style = self.particle_style[species]

            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker=style["marker"],
                    linestyle="None",
                    label=species,
                    markerfacecolor=style["color"],
                    markeredgecolor=style["color"],
                    color=style["color"],
                    markersize=10,
                )
            )

        self.particle_legend = self.ax.legend(
            handles=handles,
            title="Particles",
            loc="upper right",
        )

        self.ax.add_artist(self.particle_legend)

        if self.material_legend is not None:
            self.ax.add_artist(self.material_legend)

    def update(self):
        particles = self.step_function(self.system)

        if particles is None:
            particles = []

        current_species = self._update_particle_markers(particles)

        self._update_fp_lines(particles)

        self._update_particle_legend(current_species)

        if self.rotate:
            self.ax.view_init(azim=self.system.IT)

        self.frame_text.set_text(
            f"time = {self.system.time:.2e} {self.time_units}"
        )

        self.npart_text.set_text(
            f"npart = {len(self.system.particles)}"
        )

        return (
            list(self.material_artists.values())
            + [
                artist
                for artists in self.image_lattice_artists.values()
                for artist in artists
            ]
            + list(self.particle_artists.values())
            + list(self.fp_image_artists.values())
            + [
                self.fp_collection,
                self.frame_text,
                self.npart_text,
            ]
        )

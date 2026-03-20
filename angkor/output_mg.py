# output_mg.py
# ANGKOR - Multi-Group Output Module
# Generates professional plots and reports for SolverMG

import numpy as np
import matplotlib.pyplot as plt
import os

class OutputMG:
    """
    Generates plots and text reports from SolverMG results.
    Supports any number of energy groups G.
    """
    def __init__(self, solver, reader):
        """
        Args:
            solver (SolverMG)    : completed multi-group solver
            reader (InputReader) : input file reader
        """
        self.solver  = solver
        self.reader  = reader
        self.engine  = reader.engine
        self.title   = reader.title

        # Compute power map immediately
        self.power_map = self._compute_power()

    def _compute_power(self):
        """Compute normalized power distribution from all groups."""
        nx, ny, G = self.solver.nx, self.solver.ny, self.solver.G
        power = np.zeros((ny, nx))

        for j in range(ny):
            for i in range(nx):
                xs = self.solver._get_xs(i, j)
                # Sum over all groups: nu_sigma_f[g] * flux[g]
                for g in range(G):
                    power[j, i] += xs["nu_sigma_f"][g] * self.solver.flux[g][j, i]

        # Normalize: max power = 1.0
        if np.max(power) > 0:
            power /= np.max(power)
        return power

    def plot_all(self, save_dir=".", show=True):
        """Generate all plots: flux maps, power map, centerlines."""
        self._plot_flux_maps(save_dir, show)
        self._plot_power_map(save_dir, show)
        self._plot_centerline_profiles(save_dir, show)

    def _plot_flux_maps(self, save_dir=".", show=True):
        """Plot flux maps for all groups."""
        G = self.solver.G
        fig, axes = plt.subplots(1, G, figsize=(6*G, 5))
        if G == 1:  # Single group returns single axes object
            axes = [axes]
        fig.suptitle(f"ANGKOR — {self.title}\n"
                     f"k-eff = {self.solver.k_eff:.6f}", fontsize=13, fontweight="bold")

        extent = [0, self.engine.domain_x, 0, self.engine.domain_y]

        for g in range(G):
            im = axes[g].imshow(
                self.solver.flux[g],
                origin="lower",
                extent=extent,
                cmap="hot",
                aspect="equal"
            )
            axes[g].set_title(f"Group {g+1} Flux", fontsize=11)
            axes[g].set_xlabel("x (cm)")
            axes[g].set_ylabel("y (cm)")
            plt.colorbar(im, ax=axes[g], label="Normalized flux")

        # Add material boundaries
        for ax in axes:
            for region in self.engine.regions:
                ax.axvline(x=region.x_min, color="cyan", linewidth=0.8, linestyle="--", alpha=0.6)
                ax.axvline(x=region.x_max, color="cyan", linewidth=0.8, linestyle="--", alpha=0.6)

        plt.tight_layout()
        path = os.path.join(save_dir, "flux_maps.png")
        plt.savefig(path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {path}")
        if show:
            plt.show()
        plt.close()

    def _plot_power_map(self, save_dir=".", show=True):
        """Plot normalized power distribution."""
        fig, ax = plt.subplots(figsize=(8, 7))
        fig.suptitle(f"ANGKOR — {self.title}\n"
                     f"Power Distribution | k-eff = {self.solver.k_eff:.6f}",
                     fontsize=13, fontweight="bold")

        extent = [0, self.engine.domain_x, 0, self.engine.domain_y]

        im = ax.imshow(
            self.power_map,
            origin="lower",
            extent=extent,
            cmap="jet",
            aspect="equal",
            vmin=0
        )
        plt.colorbar(im, ax=ax, label="Normalized power (max=1.0)")
        ax.set_xlabel("x (cm)")
        ax.set_ylabel("y (cm)")
        ax.set_title("Power Peaking Map")

        # Add region boundaries
        for region in self.engine.regions:
            ax.axvline(x=region.x_min, color="white", linewidth=1.0, linestyle="--", alpha=0.7)

        plt.tight_layout()
        path = os.path.join(save_dir, "power_map.png")
        plt.savefig(path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {path}")
        if show:
            plt.show()
        plt.close()

    def _plot_centerline_profiles(self, save_dir=".", show=True):
        """Plot flux and power along x and y centerlines."""
        nx, ny, G = self.solver.nx, self.solver.ny, self.solver.G
        mid_j = ny // 2
        mid_i = nx // 2
        x_centers = self.engine.x_centers
        y_centers = self.engine.y_centers

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f"ANGKOR — {self.title} — Centerline Profiles", fontsize=13, fontweight="bold")

        # X centerline
        for g in range(G):
            axes[0].plot(x_centers, self.solver.flux[g][mid_j, :], linewidth=2, label=f"Group {g+1} Flux")
        axes[0].plot(x_centers, self.power_map[mid_j, :], "k--", linewidth=2, label="Power")
        axes[0].set_xlabel("x (cm)")
        axes[0].set_ylabel("Normalized value")
        axes[0].set_title(f"Horizontal profile at y={y_centers[mid_j]:.1f} cm")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        # Y centerline
        for g in range(G):
            axes[1].plot(y_centers, self.solver.flux[g][:, mid_i], linewidth=2, label=f"Group {g+1} Flux")
        axes[1].plot(y_centers, self.power_map[:, mid_i], "k--", linewidth=2, label="Power")
        axes[1].set_xlabel("y (cm)")
        axes[1].set_ylabel("Normalized value")
        axes[1].set_title(f"Vertical profile at x={x_centers[mid_i]:.1f} cm")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        path = os.path.join(save_dir, "centerline_profiles.png")
        plt.savefig(path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {path}")
        if show:
            plt.show()
        plt.close()

    def print_report(self):
        """Print professional text report to terminal."""
        G = self.solver.G
        print(f"\n{'='*55}")
        print(f"  ANGKOR — Simulation Report (Multi-Group)")
        print(f"{'='*55}")
        print(f"  Title    : {self.title}")
        print(f"  Domain   : {self.engine.domain_x} x {self.engine.domain_y} cm")
        print(f"  Mesh     : {self.solver.nx} x {self.solver.ny} cells")
        print(f"  Groups   : {G}")
        print(f"\n  {'─'*45}")
        print(f"  RESULTS")
        print(f"  {'─'*45}")
        print(f"  k-eff          : {self.solver.k_eff:.6f}")
        for g in range(G):
            print(f"  Max flux G{g+1}   : {np.max(self.solver.flux[g]):.4f}")
        print(f"  Max power peak  : {np.max(self.power_map):.4f}")

        # Peak location
        idx = np.unravel_index(np.argmax(self.power_map), self.power_map.shape)
        j_max, i_max = idx
        x_max = self.engine.x_centers[i_max]
        y_max = self.engine.y_centers[j_max]
        print(f"  Peak location  : ({x_max:.1f}, {y_max:.1f}) cm")
        print(f"{'='*55}\n")

    def save_report(self, save_dir="."):
        """Save text report to file."""
        path = os.path.join(save_dir, "angkor_report.txt")
        G = self.solver.G
        with open(path, "w") as f:
            f.write(f"ANGKOR Multi-Group Simulation Report\n")
            f.write(f"{'='*40}\n")
            f.write(f"Title  : {self.title}\n")
            f.write(f"k-eff  : {self.solver.k_eff:.6f}\n")
            f.write(f"Domain : {self.engine.domain_x} x {self.engine.domain_y} cm\n")
            f.write(f"Mesh   : {self.solver.nx} x {self.solver.ny}\n")
            f.write(f"Groups : {G}\n")
            f.write(f"Max Power location: ({self.engine.x_centers[np.argmax(self.power_map) % self.solver.nx]:.1f} cm)\n")
        print(f"  Saved: {path}")
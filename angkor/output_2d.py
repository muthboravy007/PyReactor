# output_2d.py
# ANGKOR - 2D Output Module
# Generates professional plots and reports

# ============================================================
# IMPORTS
# ============================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

# ============================================================
# Output2D Class
# ============================================================
class Output2D:
    """
    Generates professional plots and text reports
    from ANGKOR 2D solver results.
    """

    def __init__(self, solver, reader):
        """
        Args:
            solver (Solver2D)     : completed solver with results
            reader (InputReader)  : input file reader
        """
        self.solver  = solver          
        self.reader  = reader           
        self.engine  = self.reader.engine           
        self.title   = self.reader.title           

        # Compute power map immediately
        self.power_map = self._compute_power()

    def _compute_power(self):
        """
        Compute normalized power distribution.
        Power = nuSigf1·flux1 + nuSigf2·flux2
        """
        nx = self.solver.nx
        ny = self.solver.ny

        power = np.zeros((ny, nx))

        for j in range(ny):
            for i in range(nx):
                # Get cross sections at this cell
                xs = self.solver._get_xs(i, j)

                # Power = fission source from both groups
                power[j][i] = (
                    xs["nu_sigma_f1"] * self.solver.flux1[j][i] +
                    xs["nu_sigma_f2"] * self.solver.flux2[j][i]
                )

        # Normalize: max power = 1.0
        if np.max(power) > 0:
            power = power / np.max(power)

        return power

    def plot_all(self, save_dir=".", show=True):
        """
        Generate all plots in one call.

        Args:
            save_dir (str) : folder to save plots
            show     (bool): display plots on screen
        """
        self._plot_flux_maps(save_dir, show)
        self._plot_power_map(save_dir, show)
        self._plot_centerline_profiles(save_dir, show)

    def _plot_flux_maps(self, save_dir=".", show=True):
        """Plot group 1 and group 2 flux side by side."""

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f"ANGKOR — {self.title}\n"
                     f"k-eff = {self.solver.k_eff:.6f}",
                     fontsize=13, fontweight="bold")

        domain_x = self.engine.domain_x
        domain_y = self.engine.domain_y
        extent   = [0, domain_x, 0, domain_y]

        # --- Group 1 (Fast) ---
        im1 = axes[0].imshow(
            self.solver.flux1,
            origin="lower",
            extent=extent,
            cmap="hot",
            aspect="equal"
        )
        axes[0].set_title("Group 1 — Fast Flux", fontsize=11)
        axes[0].set_xlabel("x (cm)")
        axes[0].set_ylabel("y (cm)")
        plt.colorbar(im1, ax=axes[0], label="Normalized flux")

        # --- Group 2 (Thermal) ---
        im2 = axes[1].imshow(
            self.solver.flux2,
            origin="lower",
            extent=extent,
            cmap="hot",
            aspect="equal"
        )
        axes[1].set_title("Group 2 — Thermal Flux", fontsize=11)
        axes[1].set_xlabel("x (cm)")
        axes[1].set_ylabel("y (cm)")
        plt.colorbar(im2, ax=axes[1], label="Normalized flux")

        # Add material boundary lines
        for ax in axes:
            for region in self.engine.regions:
                ax.axvline(x=region.x_min, color="cyan",
                           linewidth=0.8, linestyle="--", alpha=0.6)
                ax.axvline(x=region.x_max, color="cyan",
                           linewidth=0.8, linestyle="--", alpha=0.6)

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
                     f"Power Distribution  |  k-eff = {self.solver.k_eff:.6f}",
                     fontsize=13, fontweight="bold")

        extent = [0, self.engine.domain_x,
                  0, self.engine.domain_y]

        im = ax.imshow(
            self.power_map,
            origin="lower",
            extent=extent,
            cmap="jet",
            aspect="equal",
            vmin=0
        )
        plt.colorbar(im, ax=ax, label="Normalized power (max=1.0)")
        ax.set_xlabel("x (cm)", fontsize=11)
        ax.set_ylabel("y (cm)", fontsize=11)
        ax.set_title("Power Peaking Map", fontsize=11)

        # Add region boundaries
        for region in self.engine.regions:
            ax.axvline(x=region.x_min, color="white",
                       linewidth=1.0, linestyle="--", alpha=0.7)

        plt.tight_layout()
        path = os.path.join(save_dir, "power_map.png")
        plt.savefig(path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {path}")
        if show:
            plt.show()
        plt.close()

    def _plot_centerline_profiles(self, save_dir=".", show=True):
        """Plot flux and power along x and y centerlines."""

        nx  = self.solver.nx
        ny  = self.solver.ny
        mid_j = ny // 2      # y-centerline row
        mid_i = nx // 2      # x-centerline column

        x_centers = self.engine.x_centers
        y_centers = self.engine.y_centers

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f"ANGKOR — {self.title} — Centerline Profiles",
                     fontsize=13, fontweight="bold")

        # --- X centerline (horizontal slice at y=midpoint) ---
        axes[0].plot(x_centers,
                     self.solver.flux1[mid_j, :],
                     "r-", linewidth=2, label="Group 1 (fast)")
        axes[0].plot(x_centers,
                     self.solver.flux2[mid_j, :],
                     "b-", linewidth=2, label="Group 2 (thermal)")
        axes[0].plot(x_centers,
                     self.power_map[mid_j, :],
                     "g--", linewidth=2, label="Power")

        # Add region boundaries
        for region in self.engine.regions:
            axes[0].axvline(x=region.x_min, color="gray",
                            linewidth=0.8, linestyle=":")
        axes[0].set_xlabel("x (cm)")
        axes[0].set_ylabel("Normalized value")
        axes[0].set_title(f"Horizontal profile at y={y_centers[mid_j]:.1f} cm")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        # --- Y centerline (vertical slice at x=midpoint) ---
        axes[1].plot(y_centers,
                     self.solver.flux1[:, mid_i],
                     "r-", linewidth=2, label="Group 1 (fast)")
        axes[1].plot(y_centers,
                     self.solver.flux2[:, mid_i],
                     "b-", linewidth=2, label="Group 2 (thermal)")
        axes[1].plot(y_centers,
                     self.power_map[:, mid_i],
                     "g--", linewidth=2, label="Power")

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
        """Print a professional text report to terminal."""

        print(f"\n{'='*55}")
        print(f"  ANGKOR — Simulation Report")
        print(f"{'='*55}")
        print(f"  Title    : {self.title}")
        print(f"  Domain   : {self.engine.domain_x} x "
              f"{self.engine.domain_y} cm")
        print(f"  Mesh     : {self.solver.nx} x {self.solver.ny} cells")
        print(f"  Groups   : 2 (fast + thermal)")
        print(f"\n  {'─'*45}")
        print(f"  RESULTS")
        print(f"  {'─'*45}")
        print(f"  k-eff          : {self.solver.k_eff:.6f}")
        print(f"  Max fast flux  : {np.max(self.solver.flux1):.4f}")
        print(f"  Max therm flux : {np.max(self.solver.flux2):.4f}")
        print(f"  Max power peak : {np.max(self.power_map):.4f}")

        # Find location of max power
        idx = np.unravel_index(
            np.argmax(self.power_map), self.power_map.shape
        )
        j_max, i_max = idx
        x_max = self.engine.x_centers[i_max]
        y_max = self.engine.y_centers[j_max]
        print(f"  Peak location  : ({x_max:.1f}, {y_max:.1f}) cm")
        print(f"{'='*55}\n")

    def save_report(self, save_dir="."):
        """Save text report to file."""
        path = os.path.join(save_dir, "angkor_report.txt")
        with open(path, "w") as f:
            f.write(f"ANGKOR Simulation Report\n")
            f.write(f"{'='*40}\n")
            f.write(f"Title  : {self.title}\n")
            f.write(f"k-eff  : {self.solver.k_eff:.6f}\n")
            f.write(f"Domain : {self.engine.domain_x} x "
                    f"{self.engine.domain_y} cm\n")
            f.write(f"Mesh   : {self.solver.nx} x {self.solver.ny}\n")
            f.write(f"Peak power location: "
                    f"({self.engine.x_centers[np.argmax(self.power_map) % self.solver.nx]:.1f} cm)\n")
        print(f"  Saved: {path}")


# ============================================================
# TEST
# ============================================================
if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from input_reader import InputReader
    from solver_2d import Solver2D

    print("=" * 55)
    print("  ANGKOR — Full Pipeline Test")
    print("=" * 55)

    # Step 1: Read input
    reader = InputReader("input/pwr_2d.yaml")
    reader.read()

    # Step 2: Solve
    solver = Solver2D(reader.engine,
                      reader.materials,
                      reader.solver)
    solver.solve()

    # Step 3: Output
    out = Output2D(solver, reader)
    out.print_report()
    out.plot_all(save_dir=".", show=True)
    out.save_report(save_dir=".")
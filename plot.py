#!/usr/bin/env python3
"""
N-Body Simulation Visualization

Creates an animated GIF showing orbital trajectories from the simulation output.
Supports both full solar system view and Earth-Moon detail view.

Author: Makar Ulesov
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from pathlib import Path

# Configuration
INPUT_FILE = "simulation_output.csv"
OUTPUT_FILE = "orbit_animation.gif"
FRAME_SKIP = 3  # Use every Nth data point (reduces file size)
FPS = 30
DPI = 100

# Astronomical unit for scaling
AU = 1.496e11  # meters


def load_data(filename: str) -> pd.DataFrame:
    """Load and validate simulation output."""
    if not Path(filename).exists():
        raise FileNotFoundError(
            f"'{filename}' not found. Run the simulation first: ./nbody"
        )
    return pd.read_csv(filename)


def create_orbit_animation(df: pd.DataFrame, output_path: str):
    """
    Create animated GIF showing orbital motion.

    Uses a two-panel layout:
    - Left: Full solar system view (Sun-Earth-Moon)
    - Right: Earth-Moon system detail view
    """
    # Extract position data (convert to AU for readability)
    sun_x = df["Sun_x"].values / AU
    sun_y = df["Sun_y"].values / AU
    earth_x = df["Earth_x"].values / AU
    earth_y = df["Earth_y"].values / AU
    moon_x = df["Moon_x"].values / AU
    moon_y = df["Moon_y"].values / AU

    # Moon position relative to Earth (in km for detail view)
    moon_rel_x = (df["Moon_x"].values - df["Earth_x"].values) / 1e6  # Mm
    moon_rel_y = (df["Moon_y"].values - df["Earth_y"].values) / 1e6  # Mm

    # Subsample for smaller GIF
    indices = range(0, len(df), FRAME_SKIP)
    n_frames = len(list(indices))

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor('#0a0a0a')

    # Style axes
    for ax in [ax1, ax2]:
        ax.set_facecolor('#0a0a0a')
        ax.tick_params(colors='white')
        ax.spines['bottom'].set_color('white')
        ax.spines['top'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.spines['right'].set_color('white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')

    # Left plot: Solar system view
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_xlabel("x [AU]")
    ax1.set_ylabel("y [AU]")
    ax1.set_title("Sun-Earth-Moon System")
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.2, color='white')

    # Right plot: Earth-Moon detail
    ax2.set_xlim(-500, 500)
    ax2.set_ylim(-500, 500)
    ax2.set_xlabel("x [Mm]")
    ax2.set_ylabel("y [Mm]")
    ax2.set_title("Moon Orbit (Earth-centered)")
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.2, color='white')

    # Initialize plot elements
    # Trails (trajectory lines)
    trail_length = 50 * FRAME_SKIP  # Show last ~50 days of trajectory

    earth_trail, = ax1.plot([], [], 'b-', alpha=0.5, linewidth=1, label='Earth')
    moon_trail_solar, = ax1.plot([], [], 'gray', alpha=0.3, linewidth=0.5)
    moon_trail_detail, = ax2.plot([], [], 'gray', alpha=0.5, linewidth=1, label='Moon')

    # Bodies
    sun_point, = ax1.plot([], [], 'yo', markersize=15, label='Sun')
    earth_point, = ax1.plot([], [], 'bo', markersize=8)
    moon_point_solar, = ax1.plot([], [], 'o', color='gray', markersize=4, label='Moon')

    earth_center, = ax2.plot([0], [0], 'bo', markersize=12, label='Earth')
    moon_point_detail, = ax2.plot([], [], 'o', color='silver', markersize=8)

    # Time display
    time_text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes,
                          color='white', fontsize=10, verticalalignment='top')

    # Legends
    ax1.legend(loc='upper right', facecolor='#1a1a1a', edgecolor='white',
               labelcolor='white', fontsize=8)
    ax2.legend(loc='upper right', facecolor='#1a1a1a', edgecolor='white',
               labelcolor='white', fontsize=8)

    def init():
        """Initialize animation."""
        return (earth_trail, moon_trail_solar, moon_trail_detail,
                sun_point, earth_point, moon_point_solar,
                moon_point_detail, time_text)

    def animate(frame_idx):
        """Update animation for each frame."""
        i = frame_idx * FRAME_SKIP

        # Calculate trail indices
        trail_start = max(0, i - trail_length)

        # Update solar system view
        earth_trail.set_data(earth_x[trail_start:i+1], earth_y[trail_start:i+1])
        moon_trail_solar.set_data(moon_x[trail_start:i+1], moon_y[trail_start:i+1])

        sun_point.set_data([sun_x[i]], [sun_y[i]])
        earth_point.set_data([earth_x[i]], [earth_y[i]])
        moon_point_solar.set_data([moon_x[i]], [moon_y[i]])

        # Update Earth-Moon detail view
        moon_trail_detail.set_data(moon_rel_x[trail_start:i+1],
                                    moon_rel_y[trail_start:i+1])
        moon_point_detail.set_data([moon_rel_x[i]], [moon_rel_y[i]])

        # Update time display
        days = df["time"].iloc[i] / (24 * 3600)
        time_text.set_text(f"Day {days:.1f}")

        return (earth_trail, moon_trail_solar, moon_trail_detail,
                sun_point, earth_point, moon_point_solar,
                moon_point_detail, time_text)

    # Create animation
    print(f"Creating animation with {n_frames} frames...")
    anim = FuncAnimation(fig, animate, init_func=init, frames=n_frames,
                         interval=1000/FPS, blit=True)

    # Save as GIF
    print(f"Saving to {output_path}...")
    writer = PillowWriter(fps=FPS)
    anim.save(output_path, writer=writer, dpi=DPI)

    plt.close()
    print(f"Animation saved: {output_path}")


def main():
    print("N-Body Simulation Visualization")
    print("=" * 35)

    # Load data
    print(f"Loading {INPUT_FILE}...")
    df = load_data(INPUT_FILE)
    print(f"Loaded {len(df)} time steps")

    # Create animation
    create_orbit_animation(df, OUTPUT_FILE)

    print("\nDone! Open orbit_animation.gif to view the result.")


if __name__ == "__main__":
    main()

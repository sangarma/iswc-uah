from manim import *
import numpy as np

class MagneticMirrorTrajectory(ThreeDScene):
    def construct(self):
        # 1. Set up the 3D Camera, Axes, and Labels
        axes = ThreeDAxes(
            x_range=[-6, 6, 1],
            y_range=[-6, 6, 1],
            z_range=[-6, 6, 1]
        )
        
        axis_labels = axes.get_axis_labels(
            x_label=Text("x", font_size=24), 
            y_label=Text("y", font_size=24), 
            z_label=Text("z", font_size=24)
        )
        
        # Set the initial view and LEAVE IT FIXED for the entire video
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.add(axes, axis_labels)

        # 2. Create a Title fixed to the screen
        title = Text("Magnetic Mirror Trajectory", color=YELLOW, font_size=28)
        title.to_corner(UL)
        self.add_fixed_in_frame_mobjects(title)

        # 3. Load the data, DOWNSAMPLE, and AUTO-SCALE it
        try:
            print("Loading massive dataset, please wait a moment...")
            # Remember to change this to your actual file name
            raw_points = np.loadtxt("mirror.txt")[:200000, 1:4]
            
            # Keep downsampling so Manim doesn't crash on 1,000,000 lines
            MAX_POINTS = 200000
            step_size = 100
            points = raw_points[::step_size] 
            print(f"Reduced {len(raw_points)} points to {len(points)} points (Step size: {step_size})")
            
        except OSError:
            raise FileNotFoundError("Make sure your data file is in the same folder.")

        # Scale the data so it fits on screen
        max_val = np.max(np.abs(points))
        if max_val > 0:
            scale_factor = 5.0 / max_val
            points *= scale_factor

        # 4. Create the invisible path
        path = VMobject().set_points_as_corners(points)

        # 5. Create the particle
        particle = Dot3D(color=RED)

        # 6. Use a ValueTracker to control movement
        tracker = ValueTracker(0.0)
        particle.add_updater(lambda m: m.move_to(path.point_from_proportion(tracker.get_value())))

        # 7. Create the trail
        trace = TracedPath(particle.get_center, stroke_color=RED, stroke_width=3)
        self.add(particle, trace)

        # 8. Set a fixed total animation time
        total_time = 12.0 

        # ---------------------------------------------------------
        # ANIMATION: Play the entire trajectory in one go (No camera moves)
        # ---------------------------------------------------------
        self.play(
            tracker.animate.set_value(1.0), # Move dot from 0% to 100%
            run_time=total_time,
            rate_func=linear
        )

        self.wait(2)
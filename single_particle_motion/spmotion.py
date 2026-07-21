from manim import *
import numpy as np

class CompareMethods(ThreeDScene):
    def construct(self):
        # 1. Set up the 3D Camera, Axes, and Labels
        # We explicitly define the range so the axes span the scaled data
        axes = ThreeDAxes(
            x_range=[-1, 1, 1],
            y_range=[-1, 1, 1],
            z_range=[-1, 10, 1]
        )
        
        # Add X, Y, Z labels to the tips of the axes
        axis_labels = axes.get_axis_labels(
            x_label=Text("x", font_size=24), 
            y_label=Text("y", font_size=24), 
            z_label=Text("z", font_size=24)
        )
        
        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.add(axes, axis_labels)

        # 2. Create a 2D Legend fixed to the screen
        legend_euler = Text("Euler Method", color=RED, font_size=24)
        legend_heun = Text("Heun Method", color=GREEN, font_size=24)
        legend_rk4 = Text("RK4 Method", color=BLUE, font_size=24)
        
        legend = VGroup(legend_euler, legend_heun, legend_rk4).arrange(DOWN, aligned_edge=LEFT)
        legend.to_corner(UL)
        self.add_fixed_in_frame_mobjects(legend)

        # 3. Load the data and AUTO-SCALE it
        try:
            points_euler = np.loadtxt("euler.txt")[:, 1:4]
            points_heun = np.loadtxt("heun.txt")[:, 1:4]
            points_rk4 = np.loadtxt("rk4.txt")[:, 1:4]
        except OSError:
            raise FileNotFoundError("Make sure the .txt files are in the same folder.")

        max_val = max(
            np.max(np.abs(points_euler)),
            np.max(np.abs(points_heun)),
            np.max(np.abs(points_rk4))
        )

        if max_val > 0:
            scale_factor = 15.0 / max_val
            points_euler *= scale_factor
            points_heun *= scale_factor
            points_rk4 *= scale_factor

        # 4. Create the invisible paths
        path_euler = VMobject().set_points_as_corners(points_euler)
        path_heun = VMobject().set_points_as_corners(points_heun)
        path_rk4 = VMobject().set_points_as_corners(points_rk4)

        # 5. Create dots for the particles
        est_euler = Dot3D(color=RED)
        est_heun = Dot3D(color=GREEN)
        est_rk4 = Dot3D(color=BLUE)

        # 6. Use a ValueTracker to control movement along the path
        # 0.0 means the start of the path, 1.0 means the end
        tracker = ValueTracker(0.0)

        # Add updaters so the dots always position themselves based on the tracker's value
        est_euler.add_updater(lambda m: m.move_to(path_euler.point_from_proportion(tracker.get_value())))
        est_heun.add_updater(lambda m: m.move_to(path_heun.point_from_proportion(tracker.get_value())))
        est_rk4.add_updater(lambda m: m.move_to(path_rk4.point_from_proportion(tracker.get_value())))

        # 7. Create the trails
        trace_euler = TracedPath(est_euler.get_center, stroke_color=RED, stroke_width=3)
        trace_heun = TracedPath(est_heun.get_center, stroke_color=GREEN, stroke_width=3)
        trace_rk4 = TracedPath(est_rk4.get_center, stroke_color=BLUE, stroke_width=3)

        self.add(est_euler, est_heun, est_rk4, trace_euler, trace_heun, trace_rk4)
        #self.add(est_euler, trace_euler)


        # 8. Calculate time (Add a safeguard so it doesn't crash if the dataset is too small)
        total_time = max(len(points_euler) * 0.05, 4.0) 
        half_time = total_time / 3

        # ---------------------------------------------------------
        # ANIMATION PART 1: First half of the trajectory
        # ---------------------------------------------------------
        self.begin_ambient_camera_rotation(rate=0.1)
        
        self.play(
            tracker.animate.set_value(0.5), # Move dots to 50% of the path
            run_time=half_time,
            rate_func=linear
        )

        self.stop_ambient_camera_rotation()

        # ---------------------------------------------------------
        # ANIMATION PART 2: Second half + Camera move to Top-Down
        # phi = 0 looks directly down the Z-axis.
        # theta = -90 aligns the X and Y axes like a standard 2D graph.
        # ---------------------------------------------------------
        self.move_camera(
                phi=0 * DEGREES, 
                theta=-90 * DEGREES,
                added_anims=[tracker.animate.set_value(1.0)], # Animates the dot simultaneously
                run_time=total_time - half_time,
                rate_func=linear
            )

        self.wait(2)